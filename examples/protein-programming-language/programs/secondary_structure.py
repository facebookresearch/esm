# Copyright (c) Meta Platforms, Inc. and affiliates.

# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

from language import (
    FixedLengthSequenceSegment,
    MaximizePLDDT,
    MaximizePTM,
    MatchSecondaryStructure,
    MinimizeSurfaceHydrophobics,
    ProgramNode,
)


def secondary_structure(
    node1_sse: str = 'a',
    node2_sse: str = 'b',
) -> ProgramNode:
    """
    Free hallucinates a protein while controlling the secondary structure
    corresponding to different segments of the sequence.
    Specify `'a'` for alpha helix, `'b'` for beta sheet, and `'c'` for coils.
    """
    node1 = ProgramNode(
        sequence_segment=FixedLengthSequenceSegment(50),
        energy_function_terms=[
            MatchSecondaryStructure(node1_sse),
            #TODO(brianhie): Add globularity here.
        ]
    )
    node2 = ProgramNode(
        sequence_segment=FixedLengthSequenceSegment(50),
        energy_function_terms=[
            MatchSecondaryStructure(node2_sse),
            #TODO(brianhie): Add globularity here.
        ]
    )

    return ProgramNode(
        energy_function_terms=[
            MaximizePTM(),
            MaximizePLDDT(),
            MinimizeSurfaceHydrophobics(),
        ],
        children=[node1, node2,],
    )
