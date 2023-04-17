# Copyright (c) Meta Platforms, Inc. and affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

from typing import Optional, NamedTuple
from enum import Enum


# Note ordering is important, correlates with Designer.py::coords[:, angle-index]
COORDS_ANGLE_NAMES = ['omega', 'theta', 'phi']
COORDS4D_NAMES = ['dist'] + COORDS_ANGLE_NAMES
COORDS6D_NAMES = COORDS4D_NAMES + ['torsion_phi', 'torsion_psi']
