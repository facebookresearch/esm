# Copyright (c) Meta Platforms, Inc. and affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

import os
import subprocess
import pytest
from pathlib import Path

notebook_dir = Path(__file__).parents[1] / "examples"
notebook_fns = notebook_dir.glob("*.ipynb")


def convert_notebook_to_py(nb_fn: Path, py_fn: Path) -> None:
    """
    From https://stackoverflow.com/questions/17077494/how-do-i-convert-a-ipython-notebook-into-a-python-file-via-commandline
    """
    import nbformat
    from nbconvert import PythonExporter

    with open(nb_fn) as fh:
        nb = nbformat.reads(fh.read(), nbformat.NO_CONVERT)

    exporter = PythonExporter()
    source, meta = exporter.from_notebook_node(nb)

    # Skip the magic, which gets converted to `get_ipython()`
    source.replace("get_ipython", "# get_ipython")

    with open(py_fn, "w+") as fh:
        fh.writelines(source)


def run_multiple(cmds):
    for cmd in cmds.strip().split("\n"):
        print(cmd)
        subprocess.run(cmd.strip(), shell=True, check=True)


def do_setup(nb_name):
    """ Do any setup work; see intro of the notebook """
    if nb_name == "sup_variant_prediction":
        cmds = """
        curl -O https://dl.fbaipublicfiles.com/fair-esm/examples/P62593_reprs.tar.gz
        tar -xzf P62593_reprs.tar.gz
        curl -O https://dl.fbaipublicfiles.com/fair-esm/examples/P62593.fasta
        """
        run_multiple(cmds)
    else:
        print(f"No setup work for {nb_name}")


@pytest.mark.parametrize("nb_fn", list(notebook_fns))
def test_run_notebook(nb_fn: Path, tmp_path: Path):
    """ Simply make sure the notebooks run from a-z """
    py_fn = tmp_path / (nb_fn.stem + ".py")
    print(py_fn)
    convert_notebook_to_py(nb_fn, py_fn)
    os.chdir(notebook_dir)
    do_setup(nb_fn.stem)
    _globals = {}
    exec(py_fn.read_text(), _globals)
    # No asserts, just running is enough for now
