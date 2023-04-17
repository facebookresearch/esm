# Copyright (c) Meta Platforms, Inc. and affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

from setuptools import setup


with open("esm/version.py") as infile:
    exec(infile.read())

with open("README.md") as f:
    readme = f.read()

extras = {
    "esmfold": [ # OpenFold does not automatically pip install requirements, so we add them here.
        "biopython",
        "deepspeed==0.5.9",
        "dm-tree",
        "pytorch-lightning",
        "omegaconf",
        "ml-collections",
        "einops",
        "scipy",
    ]
}

sources = {
    "esm": "esm",
    "esm.model": "esm/model",
    "esm.inverse_folding": "esm/inverse_folding",
    "esm.esmfold.v1": "esm/esmfold/v1",
    "esm.scripts": "scripts"
}

setup(
    name="fair-esm",
    version=version,
    description="Evolutionary Scale Modeling (esm): Pretrained language models for proteins. From Facebook AI Research.",
    long_description=readme,
    long_description_content_type="text/markdown",
    author="Facebook AI Research",
    url="https://github.com/facebookresearch/esm",
    license="MIT",
    packages=sources.keys(),
    package_dir=sources,
    extras_require=extras,
    data_files=[("source_docs/esm", ["LICENSE", "README.md", "CODE_OF_CONDUCT.rst"])],
    zip_safe=True,
    entry_points={
        "console_scripts": [
            "esm-extract=esm.scripts.extract:main",
            "esm-fold=esm.scripts.fold:main",
        ]
    },
)
