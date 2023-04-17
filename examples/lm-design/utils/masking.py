# Copyright (c) Meta Platforms, Inc. and affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

from typing import Union

import torch
import torch.nn.functional as F
from .tensor import (
    assert_shape,
    is_1hot_tensor,
)


def apply_mask(x, mask, mask_tok_idx):
    """
    Args:
        x (torch.float32): [*, K]
        mask (torch.float32): [*, 1]
        mask_tok_idx (int): vocab index for mask token
    Returns:
        x_masked (torch.float32): [*, K] = x,
            with mask_tok_idx inserted, where mask == 1.
    """
    assert is_1hot_tensor(x)
    K = x.size(-1)
    x_masked = x.detach().clone()
    x_masked[mask.squeeze(-1)] = F.one_hot(torch.tensor(mask_tok_idx), K).to(x)
    return x_masked


def assert_valid_mask(mask, x=None):
    """
    Utility that checks if a mask is [B, L, 1], given
    an optional tensor the mask might apply to, of size [B, L, K]
    This utility also returns the number of set bits in the mask, in the L dim == n.
    n is useful to know, because:
    x_masked = x.masked_select(mask)
        x: [B, L, K]
        mask: [B, L, 1]
        x_masked: [B, n, K]
    """
    if x is not None:
        assert x.dtype == torch.float, x.dtype
        assert x.ndim == 3, x.ndim
        assert mask.shape[:2] == x.shape[:2], (mask.shape[:2], x.shape[:2])
    assert mask.dtype == torch.bool
    num_masks_per_b = mask.sum(1)
    assert (num_masks_per_b[0] == num_masks_per_b).all()
    return num_masks_per_b[0]


def masked_insert_from_tokens(x, insert_mask, insert_toks):
    """ 
    [B, L, K] -> [B, I, L, K]
        x contains n masks per row.
        I is a new inserted dimension made
            by duplicating each [L, K] seq in [B] [I] times.
        Inserted tokens come from [B, I, n] insert_toks tensor.
        If insert_toks is [I, n], will broadcast.
    This is most helpful for gibbs:
        producing proposal seqs over AA_vocab.
        len(AA_vocab)**n for n-site gibbs.
    Args:
        x (torch.float32): [B, L, K]
        insert_mask (torch.bool): [B, L, 1]
            indicating which positions to insert tokens in newly created I dim.
        insert_toks (torch.float32): [(B?), I, mask_count=n]
    Returns:
        x_filled (torch.float32): [B, I, L, K]
    """
    B, L, K = assert_shape(x, -1, -1, -1)
    n = assert_valid_mask(insert_mask)
    if insert_toks.dim() == 2:
        I = assert_shape(insert_toks, -1, n) # [I, n]
        insert_toks = insert_toks.unsqueeze(0).repeat(B, 1, 1) # [B, I, n]
    else:
        I = assert_shape(insert_toks, B, -1, n) # [B, I, n]

    # 1hot -> dense (remove the need for these.)
    x_idx = x.argmax(-1) # [B, L]
    insert_mask = insert_mask.squeeze(-1) # [B, L]

    # Open up insert dim
    x_idx = x_idx.unsqueeze(1).repeat(1, I, 1) # [B, I, L]
    insert_mask = insert_mask.unsqueeze(1).repeat(1, I, 1) # [B, I, L]
    x_idx[insert_mask] = insert_toks.flatten() # [B, I, n]

    # dense -> 1hot
    x = F.one_hot(x_idx, K).to(x) # [B, I, L, K]
    return x
