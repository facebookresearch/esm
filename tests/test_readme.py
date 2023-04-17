# Copyright (c) Meta Platforms, Inc. and affiliates.
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

    model, alphabet = torch.hub.load("facebookresearch/esm:main", "esm2_t33_650M_UR50D")


def test_readme_2():
    import torch
    import esm

    # Load ESM-2 model
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    batch_converter = alphabet.get_batch_converter()
    model.eval()  # disables dropout for deterministic results

    # Prepare data (first 2 sequences from ESMStructuralSplitDataset superfamily / 4)
    data = [
        ("protein1", "MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG"),
        ("protein2", "KALTARQQEVFDLIRDHISQTGMPPTRAEIAQRLGFRSPNAAEEHLKALARKGVIEIVSGASRGIRLLQEE"),
        ("protein2 with mask","KALTARQQEVFDLIRD<mask>ISQTGMPPTRAEIAQRLGFRSPNAAEEHLKALARKGVIEIVSGASRGIRLLQEE"),
        ("protein3",  "K A <mask> I S Q"),
    ]
    batch_labels, batch_strs, batch_tokens = batch_converter(data)
    batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)

    # Extract per-residue representations (on CPU)
    with torch.no_grad():
        results = model(batch_tokens, repr_layers=[33], return_contacts=True)
    token_representations = results["representations"][33]

    # Generate per-sequence representations via averaging
    # NOTE: token 0 is always a beginning-of-sequence token, so the first residue is token 1.
    sequence_representations = []
    for i, tokens_len in enumerate(batch_lens):
        sequence_representations.append(token_representations[i, 1 : tokens_len - 1].mean(0))

    # Look at the unsupervised self-attention map contact predictions
    try:
        import matplotlib.pyplot as plt
        for (_, seq), tokens_len, attention_contacts in zip(data, batch_lens, results["contacts"]):
            plt.matshow(attention_contacts[: tokens_len, : tokens_len])
            plt.title(seq)
            plt.show()
    except ImportError:
        pass # dont need mpl to run test


def _run_py_cmd(cmd, **kwargs):
    this_python = sys.executable
    cmd.replace("python", this_python)
    subprocess.run(cmd, shell=True, check=True, **kwargs)


def test_readme_esmfold():
    import torch
    import esm

    model = esm.pretrained.esmfold_v1()
    model = model.eval().cuda()

    sequence = "MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG"
    # Multimer prediction can be done with chains separated by ':'

    with torch.no_grad():
        output = model.infer_pdb(sequence)

    with open("result.pdb", "w") as f:
        f.write(output)

    #import biotite.structure.io as bsio
    #struct = bsio.load_structure("result.pdb", extra_fields=["b_factor"])
    #print(struct.b_factor.mean())  # this will be the pLDDT
    with open("result.pdb") as f:
        lines = [line for line in f.readlines() if line.startswith('ATOM')]
    bfactors = [float(line[60:66]) for line in lines]
    assert torch.allclose(torch.Tensor(bfactors).mean(), torch.Tensor([88.3]), atol=1e-1)


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
