# Copyright (c) Meta Platforms, Inc. and affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

import logging
import math
logger = logging.getLogger(__name__)
from typing import Optional, Dict, List, Any, Union

import torch
from torch import nn
import torch.nn.functional as F
from abc import ABC, abstractmethod

from .tensor import is_1hot_tensor, add_eos_bos
from .masking import apply_mask, assert_valid_mask


def lm_marginal(mlm, x, mask):
    """ 
    Utility to extract logprobs from an MLM, given an input tensor and a mask.
    NOTE:
     - **Mask each set bit in mask's L-dim separately.** - (as opposed to all at once.)
     - Therefore, performs B*n forward passes.
     - In the future, we might try just using a single forward pass per sequence = B passes.
        - We had some success with something akin to this in the ESM-1V paper.
    Args:
        x: [B, L, K]
        mask: [B, L, 1] # n masks in L dim
            if None, use a mask of all ones.
    Returns:
        logits: [B, n, K]
    """
    B, L, K = x.shape
    if mask is None:
        mask = torch.ones(B, L, 1, dtype=torch.bool, device=x.device)
    n = assert_valid_mask(mask, x=x) # this is gross.

    # Get coords of set bits. [B*n]
    b_coords, l_coords, _ = mask.nonzero(as_tuple=True) # [B*n], [B*n]
    # Double-check of mask assumptions.
    assert torch.equal(b_coords, torch.repeat_interleave(torch.arange(B).to(x.device), n).to(x.device))

    # naming: _m1 = mask1.  (B*n leading dim.)
    x_m1 = x[b_coords]  # [B*n, L, K]
    mask_m1 = F.one_hot(l_coords, L).unsqueeze(-1).bool() # [B*n, L, 1]

    # Apply mask, mlm forward, select logits.
    x_masked_m1 = apply_mask(x_m1, mask_m1, mlm.vocab.mask_idx) # [B*n, L, K]
    lm_logits_m1 = mlm(x_masked_m1)['logits']
    lm_logits_select = lm_logits_m1.masked_select(mask_m1).reshape(B, n, K) # [B, n, K]

    # Mlm outputs 'logits' = not logprobabilities, but pre-softmax values.
    # We must log-softmax here to convert from pre-softmax logits -> logprobs.
    lm_logprobs_select = torch.log_softmax(lm_logits_select, axis=-1)

    return lm_logprobs_select # [B, n, K]


class WrapLM(nn.Module, ABC):
    def __init__(self, LM, vocab):
        super().__init__()
        self.model = LM
        self.vocab = vocab

    @abstractmethod
    def forward(self, seq1h, **kwargs):
        raise NotImplementedError


class WrapLmEsm(WrapLM):
    def forward(self, seq1h):
        B, L, K = seq1h.shape
        seq1h, seq_start_idx, seq_end_idx = self._prepare_seq(seq1h)
        seq = seq1h.argmax(-1)
        
        out = self.model(seq)

        return {
            'logits': out['logits'][:, seq_start_idx:seq_end_idx, :K],
        }
    
    def _prepare_seq(self, seq1h):
        assert is_1hot_tensor(seq1h)
        B, L, K = seq1h.shape
        # Prepend bos/cls and append eos.
        seq1h = add_eos_bos(seq1h, bos_idx=self.vocab.cls_idx, eos_idx=self.vocab.eos_idx)
        seq_start_idx = 1
        seq_end_idx = L + 1

        # In some cases, a vocab padded to 64 positions
        # with dummy character was used.
        # As a workaround, pad K dimension, then remove this portion after
        embed_dim = self.model.embed_tokens.weight.size(0)
        if embed_dim != K:
            seq1h = torch.cat([
                seq1h,
                torch.zeros(B, L+2, embed_dim - K).to(seq1h)
            ], -1)

        return seq1h, seq_start_idx, seq_end_idx
