# Copyright (c) Meta Platforms, Inc. and affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

from contextlib import contextmanager
import logging
import numpy as np
import random
from typing import Any, List, Optional, Union, Dict

import torch
import torch.nn.functional as F

from .tensor import (
    assert_logprobs,
    is_1hot_tensor,
)

logger = logging.getLogger(__name__)

"""
Mental model:
    1. This file contains everything that enables + modifies sampling from distributions.
"""

SeedDict = Dict[str, Any]
def get_rng_states() -> SeedDict:
    """ Retrieves a dictionary of seeds for all rngs. """
    rng_states = {
        'torch': torch.get_rng_state(),
        'torch.cuda': torch.cuda.get_rng_state(),
        'random': random.getstate(),
        'numpy': np.random.get_state(),
    }
    return rng_states


def set_rng_states(rng_state: SeedDict):
    """ Sets seeds of all rngs based on a dictionary of rng states. """
    torch.set_rng_state(rng_state['torch'])
    torch.cuda.set_rng_state(rng_state['torch.cuda'])
    random.setstate(rng_state['random'])
    np.random.set_state(rng_state['numpy'])


def set_rng_seeds(seed: Union[SeedDict, int]):
    """ Sets seeds of all rngs based on a single integer """
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    random.seed(seed)
    np.random.seed(seed)


@contextmanager
def set_rng_seeds_context(seeds: Union[SeedDict, int], disable=False):
    """ Context manager for temporarily setting seeds.
    For example, can be used when random sampling to get
    deterministic outputs for test cases. """
    if disable:
        yield
    else:
        rng_states = get_rng_states()
        set_rng_seeds(seeds)
        yield
        set_rng_states(rng_states)


def top_p_filtering(logits, top_p=0.0, filter_value=-float('Inf')):
    """
    Filter a distribution of logits using nucleus (top-p) filtering
    Args:
        logits: prob distribution shape (vocabulary size)
        top_p >0.0: keep the top tokens with cumulative probability >= top_p (nucleus filtering).
            Nucleus filtering is described in Holtzman et al. (http://arxiv.org/abs/1904.09751)
    Returns:
        given logits with filter_value at low probabilities (not in top-p)
    """
    sorted_logits, sorted_indices = torch.sort(logits, dim=-1, descending=True)
    cumulative_probs = torch.cumsum(F.softmax(sorted_logits, dim=-1), dim=-1)
    assert (cumulative_probs[..., -1].round() == 1).all(), "Probabilities don't add up to one"

    # Remove tokens with cumulative probability above the threshold
    sorted_indices_to_remove = cumulative_probs > top_p
    # Shift the indices to the right to keep also the first token above the threshold
    sorted_indices_to_remove[..., 1:] = sorted_indices_to_remove[..., :-1].clone()
    sorted_indices_to_remove[..., 0] = False
    # Reorder the mask to locations before sorting and apply it to logits
    sorted_indices_to_remove = sorted_indices_to_remove.gather(-1, sorted_indices.argsort(-1))
    logits[sorted_indices_to_remove] = filter_value
    return logits


def explore_ratio(probs, x_old, exploration_ratio):
    """
    Description:
        Given:
            1. a tensor with arbitrary leading dimensions
                and valid probabilities in the final dimension
            2. a tensor of same size, but 1-hot in the trailing dimension
        apply "eps-force" to all of the distributions.
        "eps-force" is a heuristic used in 2021 to prevent self-mutations.
        It does two things:
            1. It bounds the probability of positions in `probs` that
                are ==1 in `x_old` to be at most == `exploration_ratio`.
            2. If this caused a reduction in probability mass, redistribute
                that reduction (`overage`) among the other positions via
                __scaling__ the other positions multiplicatively.
                This retains the proportionality of the other positions.
    Args:
        probs (torch.float32): [*, K]
        x_old: 1 hot of previous values, to cap. [*, K]
    Example:
        Given:
            probs: [0.8, 0.15, 0.05]
            x_old: [1, 0, 0]
            exploration_ratio = 0.6
        Return:
            [0.8-0.2, 0.15*2, 0.05*2] == [0.6, 0.3, 0.1]
    """
    assert probs.shape == x_old.shape, (probs.shape, x_old.shape)
    assert is_1hot_tensor(x_old)
    thresh = 1 - exploration_ratio

    # Amount that x_old exceeded thresh, in all batch dims.
    overage = ((probs * x_old).sum(-1, keepdim=True) - thresh).clamp(0) # [*]

    # Push down old positions accoring to overage.
    probs = probs - x_old * overage # [*, K]

    # Scale up all other positions (probs_fresh) so that it sums to overage.
    probs_fresh = probs * (1-x_old) + 1e-8 # [*, K]
    probs = probs + probs_fresh * overage / probs_fresh.sum(-1, keepdim=True) # [*, K]

    return probs


def modify_logits(
    # Input:
    logits,

    # Settings:
    x_old = None,
    mask = None,
    force_propose_new_tokens = None,
    nucleus_sampling_rate = None,
    epsilon_force_mutation_ratio = None,
    epsilon_greedy_bypass = None,
    temperature = None,
):
    """
    General purpose modifier of logits/logprobs.
    Args:
        logits: [*, K] (potentially batched) tensor of distributions in final dim, K.
        x_old: [*, K] optional one-hot tensor representing current values.
            Some subsequent flags use this to promote sampling of new values.
        mask: [K] Boolean tensor for final dimension.
            Probability mass is only allowed where mask == True.
        force_propose_new_tokens (bool):  If true, don't allow any probability mass on returns logprobs where x_old == True.
        nucleus_sampling_rate (float in range [0, 1]): Fraction of nucleus sampling to use.
        epsilon_force_mutation_ratio (float in range [0, 1]): Fraction of eps-force to use.
        epsilon_greedy_bypass (float in range [0, 1]): Fraction of epsilon greedy to use.
        temperature (float): The temperature used for conversion of logits -> probs in Softmax operation.
        
    Returns:
        logprobs [*, K].  This is good for numerical precision issues, and is a bit nicer than logits.
    """
    # Ensure logprobs at the start
    logprobs = logits.log_softmax(-1)
    assert_logprobs(logprobs)
    K = logprobs.size(-1)
    if x_old is not None:
        assert logprobs.shape == x_old.shape

    # Handle optional specification of mask, which will clear some positions.
    if mask is None:
        mask = torch.ones(K, device=logprobs.device).bool()
    logits = logprobs.masked_fill(~mask, -float('inf'))

    # Do logprob-based modifications.
    if force_propose_new_tokens:
        logits = logits.masked_fill(x_old.bool(), -float('inf'))
    if nucleus_sampling_rate is not None:
        logits = top_p_filtering(logits=logits, top_p=nucleus_sampling_rate)

    # Convert to probs for prob-based modifications.
    probs = F.softmax(logits/temperature, dim=-1)

    # Do prob-based modifications.
    if epsilon_force_mutation_ratio is not None:
        probs = explore_ratio(probs, x_old, epsilon_force_mutation_ratio)
    if epsilon_greedy_bypass is not None and torch.rand(1) < epsilon_greedy_bypass: # do only for valid AA's?
        probs = torch.ones_like(probs) / K

    # 4. probs -> logits
    # NOTE: This eps used to be way too high.
    # Keep in mind that eps needs to be lower than 1e-32 for proper test case matching
    # in some cases, when adding the global EPS.
    log_probs =  (probs + 1e-100).log()
    logits = log_probs.masked_fill(~mask, -float('inf'))
    return logits