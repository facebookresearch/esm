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
    MaximizeGlobularity,
)


def symmetric_two_level_multimer(
    num_chains: int,
    num_protomers_per_chain: int,
    protomer_sequence_length: int = 50,
) -> ProgramNode:
    """
    Programs a homo-oligomeric protein with two-level symmetry.
    The number of chains in the multimer is specified by `num_chains`.
    A protomer sequence, with length `protomer_sequence_length`, is
    repeated `num_protomers_per_chain` times.
    For example, a two-chain protein with three protomers per chain
    would repeat the protomer six times.
    """

    # The basic repeated unit.
    protomer_sequence = FixedLengthSequenceSegment(protomer_sequence_length)
    def _make_protomer_node():
        return ProgramNode(sequence_segment=protomer_sequence)

    # Protomers are symmetrically combined into a chain.
    def _make_chain_node():
        return ProgramNode(
            energy_function_terms=[
                SymmetryRing(),
                MaximizeGlobularity(),
            ],
            energy_function_weights=[1., 0.05,],
            children=[
                _make_protomer_node()
                for _ in range(num_protomers_per_chain)
            ],
        )

    # Chains are symmetrically combined into a multimer.
    return ProgramNode(
        energy_function_terms=[
            MaximizePTM(),
            MaximizePLDDT(),
            SymmetryRing(),
            MinimizeSurfaceHydrophobics(),
        ],
        children=[
            _make_chain_node()
            for _ in range(num_chains)
        ],
        children_are_different_chains=True,
    )


