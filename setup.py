# Copyright (c) Facebook, Inc. and its affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

from setuptools import setup


with open("esm/version.py") as infile:
    exec(infile.read())


setup(
    name="esm",
    version=version,
    description="Pre-trained evolutionary scale models for proteins, from Facebook AI Research.",
    author="Facebook AI Research",
    url="https://github.com/facebookresearch/esm",
    license="MIT",
    packages=["esm"],
    data_files=[("source_docs/qhoptim", ["LICENSE", "README.rst", "CODE_OF_CONDUCT.rst"])],
    zip_safe=True,
)
