# Copyright (c) Meta Platforms, Inc. and affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

import pytest
from pathlib import Path
import torch
import esm

# Directly from hubconf.py
model_names = """
    esm1_t6_43M_UR50S,
    esm1_t12_85M_UR50S,
    esm1_t34_670M_UR50S,
    esm1_t34_670M_UR50D,
    esm1_t34_670M_UR100,
    esm1b_t33_650M_UR50S,
    esm_msa1_t12_100M_UR50S,
    esm_msa1b_t12_100M_UR50S,
    esm1v_t33_650M_UR90S,
    esm1v_t33_650M_UR90S_1,
    esm1v_t33_650M_UR90S_2,
    esm1v_t33_650M_UR90S_3,
    esm1v_t33_650M_UR90S_4,
    esm1v_t33_650M_UR90S_5,
    esm_if1_gvp4_t16_142M_UR50,
    esm2_t6_8M_UR50D,
    esm2_t12_35M_UR50D,
    esm2_t30_150M_UR50D,
    esm2_t33_650M_UR50D,
    esm2_t36_3B_UR50D,
    esm2_t48_15B_UR50D
"""
model_names = [mn.strip() for mn in model_names.strip(" ,\n").split(",")]


@pytest.mark.parametrize("model_name", model_names)
def test_load_hub_fwd_model(model_name: str) -> None:
    model, alphabet = getattr(esm.pretrained, model_name)()
    # batch_size = 2, seq_len = 3, tokens within vocab
    dummy_inp = torch.tensor([[0, 1, 2], [3, 4, 5]])
    if "esm_msa" in model_name:
        dummy_inp = dummy_inp.unsqueeze(0)
    output = model(dummy_inp)  # dict
    logits = output["logits"].squeeze(0)
    assert logits.shape == (2, 3, len(alphabet))


@pytest.mark.parametrize("model_name", model_names)
def test_load_local(model_name: str) -> None:
    # Assumes everything has already been loaded & cached.
    local_path = Path.home() / ".cache/torch/hub/checkpoints" / (model_name + ".pt")
    if model_name.endswith("esm1v_t33_650M_UR90S"):
        return  # skip; needs to get rerouted to specific instance
    model, alphabet = esm.pretrained.load_model_and_alphabet_local(local_path)
