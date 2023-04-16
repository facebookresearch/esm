# Copyright (c) Meta Platforms, Inc. and affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

import logging

import torch
import torch.nn.functional as F
from torch.distributions import Bernoulli, OneHotCategorical
from tqdm import tqdm


logger = logging.getLogger(__name__)


@torch.no_grad()
def stage_fixedbb(self, cfg, disable_tqdm=False):
    """Metropolis-Hastings sampling with uniform proposal and energy-based acceptance."""
    B, L, K = self.x_seqs.shape

    # import ipdb; ipdb.set_trace()

    itr = self.stepper(range(cfg.num_iter), cfg=cfg)
    itr = tqdm(itr, total=cfg.num_iter, disable=disable_tqdm)
    for step, s_cfg in itr:
        x = self.x_seqs
        a_cfg = s_cfg.accept_reject

        ##############################
        # Proposal
        ##############################
        # Decide which position to mutate == {i}.
        # mask 1 place
        mask = torch.zeros((B, L, 1), dtype=torch.bool).to(x)  # [B,L,1]
        mask[:, torch.randint(0, L, (B,))] = True  # [B,L,1]
        mask = mask.bool()

        # Uniform proposal distribution.
        log_p_x_i = torch.full((B, K), fill_value=-float("inf")).to(x)  # [B, K]
        log_p_x_i[..., self.vocab_mask_AA] = 0  # [B, K]
        p_x_i = log_p_x_i.softmax(-1)
        xp_i = OneHotCategorical(probs=p_x_i).sample()
        xp = x.masked_scatter(mask, xp_i)  # [B,L,K]

        ##############################
        # Accept / reject
        ##############################
        # log A(x',x) = log P(x') - log P(x))
        # for current input x, proposal x', target distribution P and symmetric proposal.
        log_P_x = self.calc_total_loss(x, mask, **a_cfg.energy_cfg)[0]  # [B]
        log_P_xp = self.calc_total_loss(xp, mask, **a_cfg.energy_cfg)[0]  # [B]
        log_A_xp_x = (-log_P_xp - -log_P_x) / a_cfg.temperature  # [B]
        A_xp_x = (log_A_xp_x).exp().clamp(0, 1)  # [B]
        A_bools = Bernoulli(A_xp_x).sample().bool()  # [B]
        self.x_seqs = torch.where(A_bools[:, None, None], xp, x)  # [B,L,K]
