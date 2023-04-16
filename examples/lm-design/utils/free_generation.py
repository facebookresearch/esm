# Copyright (c) Meta Platforms, Inc. and affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

from pathlib import Path
from typing import List, Optional, Dict
import torch
import tree
from omegaconf import DictConfig

from utils.scheduler import SchedulerSpec
from utils.tensor import  assert_shape
import torch.nn.functional as F
from torch.distributions.categorical import Categorical
from utils.constants import COORDS4D_NAMES
from tqdm.auto import tqdm

from utils.fixedbb import stage_fixedbb


def stage_free_generation(
    designer,
    num_iter: int,
    resample_y_every: int,
    stage_fixedbb_args: Optional[DictConfig] = None,
    resample_y_temp: SchedulerSpec = 1.0,
):
    assert resample_y_every < num_iter, \
        "resample_y_every  must be smaller than num_iter, {} {}" \
            .format(resample_y_every, num_iter)

    _T_resample_y_temp = designer.init_schedulers(resample_y_temp=resample_y_temp)[0]

    def resample_y():
        x_seq = designer.x_seqs
        struct_preds = designer.struct_model(x_seq)

        def sample_logits(preds, logits_name):
            logits = preds[logits_name]
            logits /= _T_resample_y_temp()

            # Alter preds from logits that were changed by temp (for logging after)
            distangle_key_name = logits_name.replace('_logits', '')
            preds[f'p_{distangle_key_name}'] = logits.softmax(-1)
            preds[logits_name] = logits

            sampled_map = Categorical(logits=logits).sample()
            return sampled_map

        # sample structure (dist and angles) bins from logits
        sampled_dist_and_angles = []
        for coord_name in COORDS4D_NAMES:
            sampled = sample_logits(struct_preds, f'{coord_name}_logits')
            sampled_dist_and_angles.append(sampled)
        sampled_dist_and_angles = torch.stack(sampled_dist_and_angles, dim=1)

        # Logits B x L x L x K; first 3 dims are treated as batch dims:
        set_target_structure(designer, sampled_dist_and_angles)

    curr_step = 0
    pbar = tqdm(total=num_iter, desc='stage_hallucination_joint_mh')
    while curr_step < num_iter:
        resample_y()

        # do resample_y_every steps of gibbs p(seq|struct) sampling
        num_iter_mh = resample_y_every
        stage_fixedbb_args['num_iter'] = num_iter_mh
        stage_fixedbb(designer, stage_fixedbb_args, disable_tqdm=True)
        # Set this flag to true so inner schedulers at stage_fixedbb will keep state between calls
        designer.resuming_stage = True
        curr_step += num_iter_mh
        pbar.update(num_iter_mh)

def set_target_structure(designer, sampled_dist_and_angles):
    """
    Set the given sampled contacts as the target structure in designer. This allows designing 
    a sequence for that structure later.
    """

    assert sampled_dist_and_angles.shape[0] == 1, "Only single-batch supported for now"
    if hasattr(designer, 'coords'):
        assert designer.coords.shape == sampled_dist_and_angles.shape

    # Coords is [B=1, 4, L, L]
    designer.coords = sampled_dist_and_angles

    cutoff_bin_max = designer.pdb_loader_params['contact_bin_cutoff'][1]
    # Use resnet predictons to determine the
    designer.target_contacts = (designer.coords[:, 0] <= cutoff_bin_max).squeeze(0)
    
    designer.target_no_contacts = ~designer.target_contacts
    # TEMP; assume  Batchdim==1 (also for gibbs).
    assert_shape(designer.target_contacts, designer.L, designer.L)
