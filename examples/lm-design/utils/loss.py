# Copyright (c) Meta Platforms, Inc. and affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

import torch
import torch.nn.functional as F


def get_cce_loss(probs, labels, mask=None, eps=1e-8,):
    """
    Calculates the categorical cross entropy and averages result.
    Using optional mask to control which labels are included in the loss.
    Args:
      probs (torch.float32): [B, L, L, num_categories]
      labels (torch.int32): [B, L, L]
      mask (torch.float32): [B, L, L]
    Returns:
      average_cce (torch.float32): [B]
    """
    if mask is None:
        B, L = probs.shape[:2]
        mask = torch.ones(B, L, L).to(probs.device)

    num_categories = probs.shape[-1]
    labels_onehot = F.one_hot(labels, num_categories)
    cce_ij = -torch.sum(labels_onehot * torch.log(probs+eps), axis=-1)
    average_cce = torch.sum(mask * cce_ij, axis=(1, 2)) / torch.sum(mask, axis=(1, 2))
    return average_cce
