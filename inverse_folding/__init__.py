# Copyright (c) Facebook, Inc. and its affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

from .gvp_transformer import (
    CoordBatchConverter,
    load_model_and_alphabet,
)
from .util import (
    extract_coords_from_structure,
    get_encoder_output,
    load_coords,
    load_structure,
    score_sequence,
)


def load_examples():
    import json
    from pathlib import Path
    example_file = Path(__file__).absolute().parent / "example.json"
    with open(example_file) as f:
        examples = json.load(f)
    return examples
