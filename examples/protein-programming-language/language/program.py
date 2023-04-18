# Copyright (c) Meta Platforms, Inc. and affiliates.

# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

from functools import partial
from typing import Callable, List, Optional, Tuple

import numpy as np

from language.energy import EnergyTerm
from language.folding_callbacks import FoldingResult
from language.sequence import SequenceSegmentFactory

MULTIMER_RESIDUE_INDEX_SKIP_LENGTH: int = 1000


class ProgramNode:
    def __init__(
        self,
        children: List["ProgramNode"] = None,
        sequence_segment: SequenceSegmentFactory = None,
        children_are_different_chains: bool = False,
        energy_function_terms: List[EnergyTerm] = [],
        energy_function_weights: Optional[List[float]] = None,
    ) -> None:
        self.children: Optional[List["ProgramNode"]] = children
        self.sequence_segment: SequenceSegmentFactory = sequence_segment
        self.children_are_different_chains: bool = children_are_different_chains
        self.energy_function_terms: List[energy_function_terms] = energy_function_terms
        self.energy_function_weights: List[
            float
        ] = energy_function_weights if energy_function_weights else [
            1.0 for _ in self.energy_function_terms
        ]
        if self.energy_function_weights:
            assert len(self.energy_function_terms) == len(
                self.energy_function_weights
            ), "One must have the same number of energy function terms and weights on a node."

        self.residue_index_range: Optional[Tuple[int, int]] = None

    def get_sequence_and_set_residue_index_ranges(
        self, residue_index_offset: int = 1
    ) -> Tuple[str, List[int]]:
        if self.is_leaf_node():
            sequence = self.sequence_segment.get()
            self.residue_index_range = (
                residue_index_offset,
                residue_index_offset + len(sequence),
            )
            return sequence, list(range(*self.residue_index_range))

        offset: int = residue_index_offset
        sequence = ""
        residue_indices = []
        for child in self.children:
            (
                sequence_segment,
                residue_indices_segment,
            ) = child.get_sequence_and_set_residue_index_ranges(
                residue_index_offset=offset
            )
            sequence += sequence_segment
            residue_indices += residue_indices_segment
            offset = residue_indices[-1] + 1
            if self.children_are_different_chains:
                offset += MULTIMER_RESIDUE_INDEX_SKIP_LENGTH
        self.residue_index_range = (residue_indices[0], residue_indices[-1] + 1)
        return sequence, residue_indices

    def get_residue_index_range(self) -> Tuple[int, int]:
        assert (
            self.residue_index_range
        ), "Must call get_sequence_and_set_residue_index_ranges() first."
        return self.residue_index_range

    def get_children(self) -> List["ProgramNode"]:
        return self.children

    def is_leaf_node(self) -> bool:
        return self.children is None

    def get_energy_term_functions(
        self, name_prefix: str = ""
    ) -> List[Tuple[str, float, Callable[[FoldingResult], float]]]:
        name_prefix = name_prefix if name_prefix else "root"

        terms = [
            (
                f"{name_prefix}:{type(term).__name__}",
                weight,
                partial(term.compute, self),
            )
            for weight, term in zip(
                self.energy_function_weights, self.energy_function_terms
            )
        ]

        if self.is_leaf_node():
            return terms

        for i, child in enumerate(self.children):
            terms += child.get_energy_term_functions(
                name_prefix=name_prefix + f".n{i+1}"
            )

        return terms

    def mutate(self) -> None:
        if self.is_leaf_node():
            return self.sequence_segment.mutate()

        weights = np.array(
            [float(child.num_mutation_candidates()) for child in self.children]
        )
        assert (
            weights.sum() > 0
        ), "Some mutations should be possible if mutate() was called."
        child_to_mutate = np.random.choice(self.children, p=weights / weights.sum())
        child_to_mutate.mutate()

    def num_mutation_candidates(self) -> int:
        if self.is_leaf_node():
            return self.sequence_segment.num_mutation_candidates()

        return sum([child.num_mutation_candidates() for child in self.children])
