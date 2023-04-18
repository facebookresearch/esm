# Copyright (c) Meta Platforms, Inc. and affiliates.

# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

from abc import ABC, abstractmethod
from typing import List, Optional

import numpy as np
from biotite.structure import annotate_sse, AtomArray, rmsd, sasa, superimpose

from language.folding_callbacks import FoldingResult
from language.utilities import get_atomarray_in_residue_range


class EnergyTerm(ABC):
    def __init__(self) -> None:
        pass

    @abstractmethod
    def compute(self, node, folding_result: FoldingResult) -> float:
        pass


class MaximizePTM(EnergyTerm):
    def __init__(self) -> None:
        super().__init__()

    def compute(self, node, folding_result: FoldingResult) -> float:
        del node
        return 1.0 - folding_result.ptm


class MaximizePLDDT(EnergyTerm):
    def __init__(self) -> None:
        super().__init__()

    def compute(self, node, folding_result: FoldingResult) -> float:
        del node
        return 1.0 - folding_result.plddt


class SymmetryRing(EnergyTerm):
    def __init__(self, all_to_all_protomer_symmetry: bool = False) -> None:
        super().__init__()
        self.all_to_all_protomer_symmetry: bool = all_to_all_protomer_symmetry

    def compute(self, node, folding_result: FoldingResult) -> float:
        protomer_nodes = node.get_children()
        protomer_residue_ranges = [
            protomer_node.get_residue_index_range() for protomer_node in protomer_nodes
        ]

        centers_of_mass = []
        for start, end in protomer_residue_ranges:
            backbone_coordinates = get_backbone_atoms(
                folding_result.atoms[
                    np.logical_and(
                        folding_result.atoms.res_id >= start,
                        folding_result.atoms.res_id < end,
                    )
                ]
            ).coord
            centers_of_mass.append(get_center_of_mass(backbone_coordinates))
        centers_of_mass = np.vstack(centers_of_mass)

        return (
            float(np.std(pairwise_distances(centers_of_mass)))
            if self.all_to_all_protomer_symmetry
            else float(np.std(adjacent_distances(centers_of_mass)))
        )


def get_backbone_atoms(atoms: AtomArray) -> AtomArray:
    return atoms[
        (atoms.atom_name == "CA") | (atoms.atom_name == "N") | (atoms.atom_name == "C")
    ]


def _is_Nx3(array: np.ndarray) -> bool:
    return len(array.shape) == 2 and array.shape[1] == 3


def get_center_of_mass(coordinates: np.ndarray) -> np.ndarray:
    assert _is_Nx3(coordinates), "Coordinates must be Nx3."
    return coordinates.mean(axis=0).reshape(1, 3)


def pairwise_distances(coordinates: np.ndarray) -> np.ndarray:
    assert _is_Nx3(coordinates), "Coordinates must be Nx3."
    m = coordinates[:, np.newaxis, :] - coordinates[np.newaxis, :, :]
    distance_matrix = np.linalg.norm(m, axis=-1)
    return distance_matrix[np.triu_indices(distance_matrix.shape[0], k=1)]


def adjacent_distances(coordinates: np.ndarray) -> np.ndarray:
    assert _is_Nx3(coordinates), "Coordinates must be Nx3."
    m = coordinates - np.roll(coordinates, shift=1, axis=0)
    return np.linalg.norm(m, axis=-1)


class MinimizeSurfaceHydrophobics(EnergyTerm):
    def __init__(self) -> None:
        super().__init__()

    def compute(self, node, folding_result: FoldingResult) -> float:
        start, end = node.get_residue_index_range()

        return hydrophobic_score(folding_result.atoms, start, end)


_HYDROPHOBICS = {"VAL", "ILE", "LEU", "PHE", "MET", "TRP"}


def hydrophobic_score(
    atom_array: AtomArray,
    start_residue_index: Optional[int] = None,
    end_residue_index: Optional[int] = None,
) -> float:
    """
    Computes ratio of hydrophobic atoms in a biotite AtomArray that are also surface
    exposed. Typically, lower is better.
    """

    hydrophobic_mask = np.array([aa in _HYDROPHOBICS for aa in atom_array.res_name])

    if start_residue_index is None and end_residue_index is None:
        selection_mask = np.ones_like(hydrophobic_mask)
    else:
        start_residue_index = 0 if start_residue_index is None else start_residue_index
        end_residue_index = (
            len(hydrophobic_mask) if end_residue_index is None else end_residue_index
        )
        selection_mask = np.array(
            [
                i >= start_residue_index and i < end_residue_index
                for i in range(len(hydrophobic_mask))
            ]
        )

    # TODO(scandido): Resolve the float/bool thing going on here.
    hydrophobic_surf = np.logical_and(
        selection_mask * hydrophobic_mask, sasa(atom_array)
    )
    # TODO(brianhie): Figure out how to handle divide-by-zero.
    return sum(hydrophobic_surf) / sum(selection_mask * hydrophobic_mask)


class MinimizeSurfaceExposure(EnergyTerm):
    def __init__(self) -> None:
        super().__init__()

    def compute(self, node, folding_result: FoldingResult) -> float:
        start, end = node.get_residue_index_range()

        return surface_ratio(folding_result.atoms, list(range(start, end)))


class MaximizeSurfaceExposure(EnergyTerm):
    def __init__(self) -> None:
        super().__init__()

    def compute(self, node, folding_result: FoldingResult) -> float:
        start, end = node.get_residue_index_range()

        return 1.0 - surface_ratio(folding_result.atoms, list(range(start, end)))


def surface_ratio(atom_array: AtomArray, residue_indices: List[int]) -> float:
    """Computes ratio of atoms in specified ratios which are on the protein surface."""

    residue_mask = np.array([res_id in residue_indices for res_id in atom_array.res_id])
    surface = np.logical_and(residue_mask, sasa(atom_array))
    return sum(surface) / sum(residue_mask)


class MinimizeSurfaceExposure(EnergyTerm):
    def __init__(self) -> None:
        super().__init__()

    def compute(self, node, folding_result: FoldingResult) -> float:
        start, end = node.get_residue_index_range()

        return surface_ratio(folding_result.atoms, list(range(start, end)))


class MaximizeSurfaceExposure(EnergyTerm):
    def __init__(self) -> None:
        super().__init__()

    def compute(self, node, folding_result: FoldingResult) -> float:
        start, end = node.get_residue_index_range()

        return 1.0 - surface_ratio(folding_result.atoms, list(range(start, end)))


def surface_ratio(atom_array: AtomArray, residue_indices: List[int]) -> float:
    """Computes ratio of atoms in specified ratios which are on the protein surface."""

    residue_mask = np.array([res_id in residue_indices for res_id in atom_array.res_id])
    surface = np.logical_and(residue_mask, sasa(atom_array))
    return sum(surface) / sum(residue_mask)


class MaximizeGlobularity(EnergyTerm):
    def __init__(self) -> None:
        super().__init__()

    def compute(self, node, folding_result: FoldingResult) -> float:
        start, end = node.get_residue_index_range()

        backbone = get_backbone_atoms(
            folding_result.atoms[
                np.logical_and(
                    folding_result.atoms.res_id >= start,
                    folding_result.atoms.res_id < end,
                )
            ]
        ).coord

        return float(np.std(distances_to_centroid(backbone)))


def distances_to_centroid(coordinates: np.ndarray) -> np.ndarray:
    """
    Computes the distances from each of the coordinates to the
    centroid of all coordinates.
    """
    assert _is_Nx3(coordinates), "Coordinates must be Nx3."
    center_of_mass = get_center_of_mass(coordinates)
    m = coordinates - center_of_mass
    return np.linalg.norm(m, axis=-1)


class MinimizeCRmsd(EnergyTerm):
    def __init__(self, template: AtomArray, backbone_only: bool = False) -> None:
        super().__init__()

        self.template: AtomArray = template
        self.backbone_only: bool = backbone_only
        if self.backbone_only:
            self.template = get_backbone_atoms(template)

    def compute(self, node, folding_result: FoldingResult) -> float:
        start, end = node.get_residue_index_range()

        atoms = get_atomarray_in_residue_range(folding_result.atoms, start, end)

        if self.backbone_only:
            atoms = get_backbone_atoms(atoms)

        return crmsd(self.template, atoms)


def crmsd(atom_array_a: AtomArray, atom_array_b: AtomArray) -> float:
    # TODO(scandido): Add this back.
    # atom_array_a = canonicalize_within_residue_atom_order(atom_array_a)
    # atom_array_b = canonicalize_within_residue_atom_order(atom_array_b)
    superimposed_atom_array_b_onto_a, _ = superimpose(atom_array_a, atom_array_b)
    return float(rmsd(atom_array_a, superimposed_atom_array_b_onto_a).mean())


class MinimizeDRmsd(EnergyTerm):
    def __init__(self, template: AtomArray, backbone_only: bool = False) -> None:
        super().__init__()

        self.template: AtomArray = template
        self.backbone_only: bool = backbone_only
        if self.backbone_only:
            self.template = get_backbone_atoms(template)

    def compute(self, node, folding_result: FoldingResult) -> float:
        start, end = node.get_residue_index_range()

        atoms = get_atomarray_in_residue_range(folding_result.atoms, start, end)

        if self.backbone_only:
            atoms = get_backbone_atoms(atoms)

        return drmsd(self.template, atoms)


def drmsd(atom_array_a: AtomArray, atom_array_b: AtomArray) -> float:
    # TODO(scandido): Add this back.
    # atom_array_a = canonicalize_within_residue_atom_order(atom_array_a)
    # atom_array_b = canonicalize_within_residue_atom_order(atom_array_b)

    dp = pairwise_distances(atom_array_a.coord)
    dq = pairwise_distances(atom_array_b.coord)

    return float(np.sqrt(((dp - dq) ** 2).mean()))


def pairwise_distances(coordinates: np.ndarray) -> np.ndarray:
    assert _is_Nx3(coordinates), "Coordinates must be Nx3."
    m = coordinates[:, np.newaxis, :] - coordinates[np.newaxis, :, :]
    distance_matrix = np.linalg.norm(m, axis=-1)
    return distance_matrix[np.triu_indices(distance_matrix.shape[0], k=1)]


class MatchSecondaryStructure(EnergyTerm):
    def __init__(self, secondary_structure_element: str) -> None:
        super().__init__()
        self.secondary_structure_element = secondary_structure_element

    def compute(self, node, folding_result: FoldingResult) -> float:
        start, end = node.get_residue_index_range()

        subprotein = folding_result.atoms[
            np.logical_and(
                folding_result.atoms.res_id >= start,
                folding_result.atoms.res_id < end,
            )
        ]
        sse = annotate_sse(subprotein)

        return np.mean(sse != self.secondary_structure_element)
