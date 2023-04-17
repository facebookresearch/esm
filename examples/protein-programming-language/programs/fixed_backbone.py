# Copyright (c) Meta Platforms, Inc. and affiliates.

# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

from copy import deepcopy

from biotite.database.rcsb import fetch
from biotite.structure import AtomArray

from language import (
    FixedLengthSequenceSegment,
    MaximizePLDDT,
    MaximizePTM,
    MinimizeCRmsd,
    MinimizeDRmsd,
    MinimizeSurfaceHydrophobics,
    ProgramNode,
    pdb_file_to_atomarray,
    sequence_from_atomarray,
)


def fixed_backbone_6mrs() -> ProgramNode:
    template_atoms: AtomArray = pdb_file_to_atomarray(fetch("6mrs", format="pdb"))

    sequence_length = len(sequence_from_atomarray(template_atoms))
    sequence = FixedLengthSequenceSegment(sequence_length)

    return ProgramNode(
        sequence_segment=sequence,
        energy_function_terms=[
            MaximizePTM(),
            MaximizePLDDT(),
            MinimizeSurfaceHydrophobics(),
            MinimizeCRmsd(template=template_atoms, backbone_only=True),
            MinimizeDRmsd(template=template_atoms, backbone_only=True),
        ],
    )
