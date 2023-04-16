# Copyright (c) Meta Platforms, Inc. and affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

from typing import List, Optional

import torch
import torch.nn.functional as F


def get_line_info():
    from inspect import getframeinfo, currentframe

    caller = getframeinfo(currentframe().f_back)
    return f"{caller.filename}:{caller.lineno}"


def assert_shape(tensor, *shape_args):
    # This only returns the unbound shapes
    source = f'{get_line_info()}:'
    assert len(tensor.shape) == len(shape_args), f"{source} should have {len(shape_args)} dimensions, actually has shape {tensor.shape}"
    unbound = []
    for i, ii in enumerate(shape_args):
        if ii == -1:
            unbound.append(tensor.shape[i])
        else:
            assert (
                tensor.shape[i] == ii
            ), f"{source} shape should be {shape_args}, actually {tensor.shape}"
    if len(unbound) == 1:
        return unbound[0]
    return unbound


def assert_probs(tensor: torch.FloatTensor):
    """ Assert that [*,K] input tensor is valid probs in final dim. """
    assert is_prob_tensor(tensor)


def assert_logprobs(tensor: torch.FloatTensor):
    """ Assert that [*,K] input tensor is valid logprobs in final dim. """
    try:
        assert_probs(tensor.exp())
    except AssertionError:
        raise AssertionError('Not logprobs, perhaps you have logits and need to apply F.log_softmax?')


def is_1hot_tensor(x1h):
    return x1h.min() == 0 and x1h.max() == 1 and \
           x1h.shape[-1] > 1 and (x1h.sum(-1) == 1).all()


def is_prob_tensor(p, atol=1e-3):
    return (
        p.min() >= 0
        and p.max() <= 1
        and torch.isclose(
            p.sum(axis=-1), torch.ones_like(p.sum(axis=-1)), atol=atol).all()
    ).item()


def add_eos_bos(
    seq1h,
    bos_idx: Optional[int] = None,
    eos_idx: Optional[int] = None,
):
    """
    Helper for (possibly) prepending bos/cls and appending eos tokens.
    """
    B, L, K = seq1h.shape

    to_concat = []
    if bos_idx is not None:
        to_concat.append(F.one_hot(torch.full([B,1], bos_idx), K).to(seq1h))
    to_concat.append(seq1h)
    if eos_idx is not None:
        to_concat.append(F.one_hot(torch.full([B,1], eos_idx), K).to(seq1h))
    seq1h = torch.cat(to_concat, axis=1)

    return seq1h
