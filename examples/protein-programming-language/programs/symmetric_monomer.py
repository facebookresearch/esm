# Copyright (c) Meta Platforms, Inc. and affiliates.

# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

from language import (
    FixedLengthSequenceSegment,
    MaximizePLDDT,
    MaximizePTM,
    MinimizeSurfaceHydrophobics,
    ProgramNode,
    SymmetryRing,
)


def symmetric_monomer(num_protomers: int) -> ProgramNode:
    protomer_sequence = FixedLengthSequenceSegment(50)
    def _make_protomer_node():
        # A new ProgramNode must be made for each new protomer,
        # but the sequence can (and should) be shared.
        return ProgramNode(
            sequence_segment=protomer_sequence
        )

    return ProgramNode(
        energy_function_terms=[
            MaximizePTM(),
            MaximizePLDDT(),
            SymmetryRing(),
            MinimizeSurfaceHydrophobics(),
        ],
        children=[
            _make_protomer_node()
            for _ in range(num_protomers)
        ],
    )
