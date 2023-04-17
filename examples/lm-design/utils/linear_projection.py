# Copyright (c) Meta Platforms, Inc. and affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.
#
import logging
from collections import defaultdict
from functools import partial
from pathlib import Path
import torch
import torch.nn.functional as F
import torch.nn as nn
from collections import namedtuple

import os
import sys
import torch
from pathlib import Path
from argparse import Namespace
from omegaconf import DictConfig, OmegaConf
from typing import Dict, Optional, List, Callable, Union
import math

BinningScheme = namedtuple(
    'BinningScheme',
    [
        'N_BINS',
        'CUTOFF_BIN',
        'MIN_DIST',
        'MAX_DIST',
        'CONTACT_DIST',
        'THETA_BINS',
        'PHI_BINS',
        'OMEGA_BINS',
        'TORSION_BINS',
    ],
)
LinearProjectionDist_1A = BinningScheme(
    N_BINS=18,
    CUTOFF_BIN=5,  # conservative cutoff for 8 angstroms exact ((20-2.5)/16*5 + 2.5 = ~7.9)
    MIN_DIST=2.5,
    MAX_DIST=20,
    CONTACT_DIST=8,
    THETA_BINS=18,
    PHI_BINS=8,
    OMEGA_BINS=18,
    TORSION_BINS=50,
)
BIN_FLAG_TO_ENUM = {
    'LinearProjectionDist_1A': LinearProjectionDist_1A,
}

logger = logging.getLogger(__name__)


def extract_features(
    model,
    inp: torch.Tensor,
    has_cls: bool = True,
    has_eos: bool = True,
    need_head_weights: Union[bool, Dict] = False,
):
    """
    Inflexible way of dealing with the mess of inputs with cls tokens,
    need embedding without, and LSTM doesnt want cls token in the first place.
    """
    
    out_start_idx = 1 if has_cls else 0
    out_end_idx = -1 if has_eos else None
    inp = inp.argmax(dim=-1) # output [B, T]
    result = model(inp, need_head_weights=need_head_weights)
    attentions = result["attentions"]

    batch, layer, head, seqlen, seqlen2 = attentions.size()
    assert seqlen == seqlen2
    attentions = attentions.reshape(
        batch, layer * head, seqlen, seqlen
    )
    emb = result["logits"]

    return emb[:, out_start_idx:out_end_idx], \
        attentions[:, :, out_start_idx:out_end_idx, out_start_idx:out_end_idx]

    
class LinearProjectionDistogramModel(nn.Module):
    """
    This model regresses angles and distances from the attention maps of the LM
    """

    def __init__(self):
        super().__init__()
        self.base_model = None  # To be initialized later
        self.num_heads = 20 * 33
        self.bin_enum = BIN_FLAG_TO_ENUM['LinearProjectionDist_1A']
        self._s1 = self.bin_enum.N_BINS
        self._s2 = self._s1 + self.bin_enum.THETA_BINS
        self._s3 = self._s2 + self.bin_enum.PHI_BINS
        self._s4 = self._s3 + self.bin_enum.OMEGA_BINS
        
        conv_in_channels = self.num_heads

        self.conv1 = torch.nn.Conv2d(
            in_channels=conv_in_channels,
            out_channels=self.bin_enum.N_BINS + self.bin_enum.OMEGA_BINS,
            kernel_size=1,
            stride=1,
            padding=0,
        )

        self.conv2 = torch.nn.Conv2d(
            in_channels=conv_in_channels,
            out_channels=self.bin_enum.THETA_BINS + self.bin_enum.PHI_BINS,
            kernel_size=1,
            stride=1,
            padding=0,
        )


    def forward(self, src_tokens, **kwargs):
        _, attentions_asym = extract_features(self.base_model, src_tokens, need_head_weights=True)
            
        def symmetrize(contacts, scale=1.0):
            return scale * (contacts + contacts.transpose(-1, -2))
             
        attentions_sym = symmetrize(attentions_asym)
        
        # [B, C, N, N] -> [B, output_channels, N, N]
        out1 = self.conv1(attentions_sym)
        out2 = self.conv2(attentions_asym)
        return {
            'logits':       out1[:, :self.bin_enum.N_BINS:,    :, :],
            'omega_logits': out1[:, self.bin_enum.N_BINS:,     :, :],
            'theta_logits': out2[:, :self.bin_enum.THETA_BINS, :, :],
            'phi_logits':   out2[:, self.bin_enum.THETA_BINS:, :, :],
        }



