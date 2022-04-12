# Copyright (c) Facebook, Inc. and its affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

import sys
import subprocess
import tempfile
import requests
import shutil
from pathlib import Path
import torch
import esm


def test_readme_1():
    import torch

    model, alphabet = torch.hub.load("facebookresearch/esm:main", "esm1b_t33_650M_UR50S")


def test_readme_2():
    import torch
    import esm

    # Load ESM-1b model
    model, alphabet = esm.pretrained.esm1b_t33_650M_UR50S()
    batch_converter = alphabet.get_batch_converter()
    model.eval()  # disables dropout for deterministic results

    # Prepare data (first 2 sequences from ESMStructuralSplitDataset superfamily / 4)
    data = [
        ("protein1", "MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG"),
        ("protein2", "KALTARQQEVFDLIRDHISQTGMPPTRAEIAQRLGFRSPNAAEEHLKALARKGVIEIVSGASRGIRLLQEE"),
        (
            "protein2 with mask",
            "KALTARQQEVFDLIRD<mask>ISQTGMPPTRAEIAQRLGFRSPNAAEEHLKALARKGVIEIVSGASRGIRLLQEE"
        ),
        (
            "protein3",
            "K A <mask> I S Q"
        ),
    ]
    batch_labels, batch_strs, batch_tokens = batch_converter(data)

    # Extract per-residue representations (on CPU)
    with torch.no_grad():
        results = model(batch_tokens, repr_layers=[33], return_contacts=True)
    token_representations = results["representations"][33]

    # Generate per-sequence representations via averaging
    # NOTE: token 0 is always a beginning-of-sequence token, so the first residue is token 1.
    sequence_representations = []
    for i, (_, seq) in enumerate(data):
        sequence_representations.append(token_representations[i, 1 : len(seq) + 1].mean(0))

    # Look at the unsupervised self-attention map contact predictions
    import matplotlib.pyplot as plt

    for (_, seq), attention_contacts in zip(data, results["contacts"]):
        plt.matshow(attention_contacts[: len(seq), : len(seq)])
        plt.title(seq)
        plt.show()


def _run_py_cmd(cmd, **kwargs):
    this_python = sys.executable
    cmd.replace("python", this_python)
    subprocess.run(cmd, shell=True, check=True, **kwargs)


def test_readme_3():
    # NOTE modification on copy paste from README for speed:
    # * some_proteins -> few_proteins (subset)
    # * I computed reference values a while ago for: esm1b -> esm1 and layers 33 -> 34
    cmd = """
python scripts/extract.py esm1_t34_670M_UR50S examples/data/few_proteins.fasta examples/data/few_proteins_emb_esm1/ \
    --repr_layers 0 33 34 --include mean per_tok
"""
    _run_py_cmd(cmd)
    confirm_all_tensors_equal(
        "examples/few_proteins_emb_esm1/",
        "https://dl.fbaipublicfiles.com/fair-esm/tests/some_proteins_emb_esm1_t34_670M_UR50S_ref",
    )


def assert_pt_file_equal(f, fref):
    a = torch.load(f)
    b = torch.load(fref)
    # set intersection of dict keys:
    which_layers = a["representations"].keys() & b["representations"].keys()
    assert which_layers, "Expected at least one layer appearing in both dumps"
    for layer in which_layers:
        assert torch.allclose(a["representations"][layer], b["representations"][layer], atol=1e-3)


def confirm_all_tensors_equal(local_dir: str, ref_dir: str) -> None:
    # TODO use pytest built-in fixtures for tmp_path https://docs.pytest.org/en/6.2.x/fixture.html#fixtures
    for fn in Path(local_dir).glob("*.pt"):
        with tempfile.NamedTemporaryFile(mode="w+b", prefix=fn.name) as f:
            ref_url = f"{ref_dir}/{fn.name}"
            with requests.get(ref_url, stream=True) as r:
                shutil.copyfileobj(r.raw, f)
            f.seek(0)
            assert_pt_file_equal(fn, f)


def test_msa_transformers():
    _test_msa_transformer(*esm.pretrained.esm_msa1_t12_100M_UR50S())
    _test_msa_transformer(*esm.pretrained.esm_msa1b_t12_100M_UR50S())


def _test_msa_transformer(model, alphabet):
    batch_converter = alphabet.get_batch_converter()
    # Make an "MSA" of size 3
    data = [
        ("protein1", "MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG"),
        ("protein2", "MHTVRQSRLKSIVRILEMSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG"),
        ("protein3", "MHTVRQSRLKSIVRILEMSKEPVSGAQL---LSVSRQVIVQDIAYLRSLGYNIVAT----VLAGG"),
    ]
    batch_labels, batch_strs, batch_tokens = batch_converter(data)

    with torch.no_grad():
        results = model(batch_tokens, repr_layers=[12], return_contacts=True)
    token_representations = results["representations"][12]
    assert token_representations.shape == (1, 3, 66, 768)


def test_variant_readme_1():
    cmd = """
python predict.py \
    --model-location esm1v_t33_650M_UR90S_1 esm1v_t33_650M_UR90S_2 esm1v_t33_650M_UR90S_3 esm1v_t33_650M_UR90S_4 esm1v_t33_650M_UR90S_5 \
    --sequence HPETLVKVKDAEDQLGARVGYIELDLNSGKILESFRPEERFPMMSTFKVLLCGAVLSRVDAGQEQLGRRIHYSQNDLVEYSPVTEKHLTDGMTVRELCSAAITMSDNTAANLLLTTIGGPKELTAFLHNMGDHVTRLDRWEPELNEAIPNDERDTTMPAAMATTLRKLLTGELLTLASRQQLIDWMEADKVAGPLLRSALPAGWFIADKSGAGERGSRGIIAALGPDGKPSRIVVIYTTGSQATMDERNRQIAEIGASLIKHW \
    --dms-input ./data/BLAT_ECOLX_Ranganathan2015.csv \
    --mutation-col mutant \
    --dms-output ./data/BLAT_ECOLX_Ranganathan2015_labeled.csv \
    --offset-idx 24 \
    --scoring-strategy wt-marginals
    """
    _run_py_cmd(cmd, cwd="examples/variant-prediction/")


def test_variant_readme_2():
    cmd = """
python predict.py \
    --model-location esm_msa1b_t12_100M_UR50S \
    --sequence HPETLVKVKDAEDQLGARVGYIELDLNSGKILESFRPEERFPMMSTFKVLLCGAVLSRVDAGQEQLGRRIHYSQNDLVEYSPVTEKHLTDGMTVRELCSAAITMSDNTAANLLLTTIGGPKELTAFLHNMGDHVTRLDRWEPELNEAIPNDERDTTMPAAMATTLRKLLTGELLTLASRQQLIDWMEADKVAGPLLRSALPAGWFIADKSGAGERGSRGIIAALGPDGKPSRIVVIYTTGSQATMDERNRQIAEIGASLIKHW \
    --dms-input ./data/BLAT_ECOLX_Ranganathan2015.csv \
    --mutation-col mutant \
    --dms-output ./data/BLAT_ECOLX_Ranganathan2015_labeled.csv \
    --offset-idx 24 \
    --scoring-strategy masked-marginals \
    --msa-path ./data/BLAT_ECOLX_1_b0.5.a3m
    """
    _run_py_cmd(cmd, cwd="examples/variant-prediction/")


if __name__ == "__main__":
    confirm_all_tensors_equal(
        "examples/few_proteins_emb_esm1/",
        "https://dl.fbaipublicfiles.com/fair-esm/tests/some_proteins_emb_esm1_t34_670M_UR50S_ref/",
    )
