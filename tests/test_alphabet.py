# Copyright (c) Meta Platforms, Inc. and affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

def _test_esm1b(alphabet):
    import torch

    batch_converter = alphabet.get_batch_converter()

    data = [
        ("protein1", "MKTVRQG"),
        ("protein2 with mask", "KALTA<mask>ISQP"),
        ("protein3", "K A <mask> I S Q"),
    ]
    _, _, batch_tokens = batch_converter(data)
    expected_tokens = torch.tensor(
        [
            [0, 20, 15, 11, 7, 10, 16, 6, 2, 1, 1, 1],
            [0, 15, 5, 4, 11, 5, 32, 12, 8, 16, 14, 2],
            [0, 15, 5, 32, 12, 8, 16, 2, 1, 1, 1, 1],
        ]
    )
    assert torch.allclose(batch_tokens, expected_tokens)


def _test_esm1b_truncation(alphabet):
    import torch

    batch_converter = alphabet.get_batch_converter(truncation_seq_length=10)

    data = [
        ("protein1", "MKTVRQGMKTVRQG"),
        ("protein2 with mask", "KALTA<mask>ISQPISQP"),
        ("protein3", "K A <mask> I S Q"),
    ]
    _, _, batch_tokens = batch_converter(data)
    expected_tokens = torch.tensor(
        [
            [0, 20, 15, 11, 7, 10, 16, 6, 20, 15, 11, 2],
            [0, 15, 5, 4, 11, 5, 32, 12, 8, 16, 14, 2],
            [0, 15, 5, 32, 12, 8, 16, 2, 1, 1, 1, 1],
        ]
    )
    assert torch.allclose(batch_tokens, expected_tokens)


def test_esm1b_alphabet():
    import esm

    _, alphabet = esm.pretrained.esm1b_t33_650M_UR50S()
    _test_esm1b(alphabet)
    _test_esm1b_truncation(alphabet)


def test_esm1v_alphabet():
    import esm

    _, alphabet = esm.pretrained.esm1v_t33_650M_UR90S_1()
    _test_esm1b(alphabet)
    _test_esm1b_truncation(alphabet)


def test_esm1_msa1b_alphabet():
    import torch
    import esm

    # Load ESM-1b model
    _, alphabet = esm.pretrained.esm_msa1b_t12_100M_UR50S()
    batch_converter = alphabet.get_batch_converter()

    data = [
        ("protein1", "MKTVRQG"),
        ("protein2", "KALTRAI"),
        ("protein3", "KAAISQQ"),
    ]
    _, _, batch_tokens = batch_converter(data)
    expected_tokens = torch.tensor(
        [
            [
                [0, 20, 15, 11, 7, 10, 16, 6],
                [0, 15, 5, 4, 11, 10, 5, 12],
                [0, 15, 5, 5, 12, 8, 16, 16],
            ]
        ]
    )
    assert torch.allclose(batch_tokens, expected_tokens)
