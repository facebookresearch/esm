# Copyright (c) Meta Platforms, Inc. and affiliates.

# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

from copy import deepcopy

from biotite.database.rcsb import fetch
from biotite.structure import AtomArray

from language import (
    ConstantSequenceSegment,
    MaximizePLDDT,
    MaximizePTM,
    MaximizeSurfaceExposure,
    MinimizeCRmsd,
    MinimizeDRmsd,
    MinimizeSurfaceHydrophobics,
    ProgramNode,
    VariableLengthSequenceSegment,
    get_atomarray_in_residue_range,
    pdb_file_to_atomarray,
    sequence_from_atomarray,
)


def scaffolding_ace2() -> ProgramNode:
    binding_site_atoms: AtomArray = pdb_file_to_atomarray(fetch("6m0j", format="pdb"))
    binding_site_atoms = get_atomarray_in_residue_range(
        binding_site_atoms, start=23, end=42
    )
    binding_site_sequence: str = sequence_from_atomarray(binding_site_atoms)
    print(binding_site_sequence)

    leader_amino_acid_sequence = VariableLengthSequenceSegment(50)
    binding_site_sequence = ConstantSequenceSegment(binding_site_sequence)
    follower_amino_acid_sequence = VariableLengthSequenceSegment(50)

    return ProgramNode(
        energy_function_terms=[
            MaximizePTM(),
            MaximizePLDDT(),
            MinimizeSurfaceHydrophobics(),
        ],
        children=[
            ProgramNode(sequence_segment=leader_amino_acid_sequence),
            ProgramNode(
                sequence_segment=binding_site_sequence,
                energy_function_terms=[
                    MaximizeSurfaceExposure(),
                    MinimizeCRmsd(template=binding_site_atoms),
                    MinimizeDRmsd(template=binding_site_atoms),
                ],
                energy_function_weights=[1.0, 10.0, 10.0],
            ),
            ProgramNode(sequence_segment=follower_amino_acid_sequence),
        ]
    )
