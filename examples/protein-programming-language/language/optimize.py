# Copyright (c) Meta Platforms, Inc. and affiliates.

# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

from copy import deepcopy
from dataclasses import dataclass

import numpy as np
from rich.live import Live
from rich.table import Table

from language.folding_callbacks import FoldingCallback
from language.program import ProgramNode


@dataclass
class MetropolisHastingsState:
    program: ProgramNode
    temperature: float
    annealing_rate: float
    num_steps: int
    energy: float
    best_energy: float
    energy_term_fn_values: list


def metropolis_hastings_step(
    state: MetropolisHastingsState,
    folding_callback: FoldingCallback,
    verbose: bool = False,
) -> MetropolisHastingsState:
    temperature = state.temperature * state.annealing_rate

    candidate: ProgramNode = deepcopy(state.program)
    candidate.mutate()

    sequence, residue_indices = candidate.get_sequence_and_set_residue_index_ranges()
    folding_output = folding_callback.fold(sequence, residue_indices)

    energy_term_fns = candidate.get_energy_term_functions()
    energy_term_fn_values = [
        (name, weight, energy_fn(folding_output)) for name, weight, energy_fn in energy_term_fns
    ]
    # TODO(scandido): Log these.
    energy: float = sum([
        weight * value for _, weight, value in energy_term_fn_values
    ])

    accept_candidate = False
    if state.energy is None:
        accept_candidate = True
    else:
        # NOTE(scandido): We are minimizing the function here so instead of
        # candidate - current we do -1 * (candidate - current) = -candidate + current.
        energy_differential: float = -energy + state.energy
        accept_probability: float = np.clip(
            # NOTE(scandido): We approximate the ratio of transition probabilities from
            # current to candidate vs. candidate to current to be equal, which is
            # approximately correct.
            np.exp(energy_differential / temperature),
            a_min=None,
            a_max=1.0,
        )
        accept_candidate: bool = np.random.uniform() < accept_probability

    if accept_candidate:
        sequence, _ = candidate.get_sequence_and_set_residue_index_ranges()
        if verbose:
            print(f"Accepted {sequence} with energy {energy:.2f}.")

    return MetropolisHastingsState(
        program=candidate if accept_candidate else state.program,
        temperature=temperature,
        annealing_rate=state.annealing_rate,
        num_steps=state.num_steps + 1,
        energy=energy if accept_candidate else state.energy,
        best_energy=min(energy, state.energy) if state.energy else energy,
        energy_term_fn_values=energy_term_fn_values,
    )


def run_simulated_annealing(
    program: ProgramNode,
    initial_temperature: float,
    annealing_rate: float,
    total_num_steps: int,
    folding_callback: FoldingCallback,
    display_progress: bool = True,
    progress_verbose_print: bool = False,
) -> ProgramNode:
    # TODO(scandido): Track accept rate.

    state = MetropolisHastingsState(
        program=program,
        temperature=initial_temperature,
        annealing_rate=annealing_rate,
        num_steps=0,
        energy=None,
        best_energy=None,
        energy_term_fn_values=None,
    )

    def _generate_table(state):
        table = Table()
        table.add_column("Energy name")
        table.add_column("Weight")
        table.add_column("Value")
        if state.energy_term_fn_values is None:
            return table
        for name, weight, value in state.energy_term_fn_values:
            table.add_row(name, f"{weight:.2f}", f"{value:.2f}")
        table.add_row("Energy", "", f"{state.energy:.2f}")
        table.add_row(
            "Iterations", "",
            f"{state.num_steps} / {total_num_steps}"
        )
        return table

    with Live() as live:
        for _ in range(1, total_num_steps + 1):
            state = metropolis_hastings_step(
                state,
                folding_callback,
                verbose=progress_verbose_print,
            )
            if display_progress:
                live.update(_generate_table(state))

    return state.program
