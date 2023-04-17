# Copyright (c) Meta Platforms, Inc. and affiliates.

# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

from language import (
    FixedLengthSequenceSegment,
    MaximizePLDDT,
    MaximizePTM,
    MinimizeSurfaceHydrophobics,
    ProgramNode,
)


def free_hallucination(sequence_length: int) -> ProgramNode:
    sequence = FixedLengthSequenceSegment(sequence_length)
    return ProgramNode(
        sequence_segment=sequence,
        energy_function_terms=[
            MaximizePTM(),
            MaximizePLDDT(),
            MinimizeSurfaceHydrophobics(),
        ],
    )
