# Copyright (c) Meta Platforms, Inc. and affiliates.

# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

from abc import ABC, abstractmethod
from dataclasses import dataclass
from io import StringIO
from typing import List

import esm
import torch
from biotite.structure import AtomArray
import numpy as np
from openfold.np.residue_constants import atom_order
from torch.utils._pytree import tree_map

from language.utilities import pdb_file_to_atomarray


@dataclass
class FoldingResult:
    atoms: AtomArray
    ptm: float
    plddt: float


class FoldingCallback(ABC):
    "Interface for running ESMFold and other folding methods."

    def __init__(self) -> None:
        pass

    @abstractmethod
    def load(self, device: str) -> None:
        pass

    @abstractmethod
    def fold(self, sequence: str, residue_indices: List[int]) -> FoldingResult:
        pass


class EsmFoldv1(FoldingCallback):
    "Runs ESMFold v1.0."

    def __init__(self) -> None:
        super().__init__()

        self.model = None

    def load(self, device: str) -> None:
        self.model = esm.pretrained.esmfold_v1().eval()
        self.model = self.model.to(device)

    def fold(self, sequence: str, residue_indices: List[int]) -> FoldingResult:
        assert self.model is not None, "Must call load() before fold()."

        # TODO: Current `esm.esmfold.v1.misc.output_to_pdb()` adds 1 to the `residx`
        # mistakenly, just subtract 1 for now but fix in a later version.
        residue_indices = np.array(residue_indices) - 1

        raw_output = self.model.infer(
            sequence, residx=torch.Tensor(residue_indices).long().reshape(1, -1),
        )
        raw_output = tree_map(lambda x: x.to("cpu"), raw_output)

        pdb_string = esm.esmfold.v1.misc.output_to_pdb(raw_output)[0]
        atoms: AtomArray = pdb_file_to_atomarray(StringIO(pdb_string))

        plddt = raw_output["plddt"]
        plddt = plddt[0, ...].numpy()
        plddt = plddt.transpose()
        plddt = plddt[atom_order["CA"], :]
        plddt = float(plddt.mean()) / 100.0

        ptm = float(raw_output["ptm"])

        return FoldingResult(atoms=atoms, ptm=ptm, plddt=plddt)
