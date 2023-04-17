# Copyright (c) Meta Platforms, Inc. and affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

import logging
from typing import Tuple, Dict
import torch
from torch import nn
from esm.data import Alphabet
from utils.tensor import add_eos_bos, is_1hot_tensor
from torch import nn
from torch.nn import functional as F
from abc import ABC, abstractmethod
import hydra
import os 

from utils.linear_projection import LinearProjectionDistogramModel

logger = logging.getLogger(__name__)  # Hydra configured

PDB_LOADER_PARAM_REGISTRY = {
    'LinearProjectionDist-1A': {
        "DMIN": 2.5,
        "DMAX": 20.0,
        "DBINS": 17,
        "ABINS": 17,
        "PHI_BINS": 7,
        "WMIN": 0.8,
        "LMIN": 150,
        "LMAX": 400,
        "EPOCHS": 10,
        "NCPU": 8,
        "SLICE": "CONT",
        "contact_bin_cutoff": (0, 5)  # both inclusive
    },
}


class WrapStruct(nn.Module, ABC):
    """
    Wraps structure models used by Designer.
    Companion of WrapLM.
    """
    def __init__(self, model, vocab, has_eos=False):
        super().__init__()
        self.model = model
        self.vocab = vocab
        self.has_eos = has_eos

    @abstractmethod
    def forward(self, seq1h, pos_idx=None):
        raise NotImplementedError

    @property
    def msa_based(self):
        return False


def load(vocab: Alphabet) -> Tuple[nn.Module, Dict]:
        
    pdb_loader_params = PDB_LOADER_PARAM_REGISTRY['LinearProjectionDist-1A']
    model = load_model()
    has_eos = True
    return WrapStructModel(model, vocab, has_eos), pdb_loader_params


def load_model():
    model_url_path = 'https://dl.fbaipublicfiles.com/fair-esm/examples/lm_design/linear_projection_model.pt'
    local_model_path = hydra.utils.to_absolute_path('./linear_projection_model.pt')
    if not os.path.exists(local_model_path):
        logger.info(f'Downloading linear projection model from {model_url_path} to {local_model_path}')
        os.system(f'wget {model_url_path} -O {local_model_path}')
        
    # load model_path
    state = torch.load(local_model_path, map_location='cpu')
    # chkpt_args = state['cfg']
    model = LinearProjectionDistogramModel()
    model_state = state['model']
    model.load_state_dict(model_state)
    from esm.pretrained import esm2_t33_650M_UR50D
    base_model, alphabet = esm2_t33_650M_UR50D()
    model.base_model = base_model
    return model
   

class WrapStructModel(WrapStruct):
    def forward(self, seq1h, pos_idx=None):
        assert is_1hot_tensor(seq1h)
        B, L, K = seq1h.shape
        # Prepend bos/cls only.
        eos_idx = self.vocab.eos_idx if self.has_eos else None
        seq1h = add_eos_bos(seq1h, bos_idx=self.vocab.cls_idx, eos_idx=eos_idx)
        # In designer, all seqs have same length.
        src_lengths = torch.full((B,), L + 1 + int(self.has_eos)).to(seq1h).long()

        # In some cases, a vocab padded to 64 positions
        # with dummy character was used.
        # As a workaround, pad K dimensions up to 64 (embed_dim),
        # then remove this portion after to get K tokens back.
        embed_dim = self.model.base_model.embed_tokens.weight.size(0)
        if embed_dim != K:
            seq1h = torch.cat([
                seq1h,
                torch.zeros(B, L+1, embed_dim - K).to(seq1h)
            ], -1)

        # Call forward.  bos/eos will be stripped inside.
        out_dict = self.model(seq1h, src_lengths=src_lengths)
        return self._process_output(out_dict)

    @property
    def lm(self):
        return self.model.base_model


    @staticmethod
    def _process_output(out_dict):
        for k in ('logits', 'theta_logits', 'phi_logits', 'omega_logits'):
            out_dict[k] = out_dict[k].permute(0, 2, 3, 1).contiguous()
        
        if 'logits' in out_dict:
            out_dict['dist_logits'] = out_dict['logits']
            del out_dict['logits']
        for k in [k for k in out_dict if k.endswith('_logits')]:
            targetname = k[:-len('_logits')]
            out_dict[f'p_{targetname}'] = F.softmax(out_dict[k], -1)
        return out_dict
