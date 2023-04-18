# Protein programming language

This repository contains an implementation of the programming language described in the paper [A high-level programming language for generative protein design](https://www.biorxiv.org/content/10.1101/2022.12.21.521526v1).

## Tutorial

We have provided a [tutorial notebook](tutorial.ipynb) describing the basics of writing programs and running optimization loops.
The notebook can be run in Colab:

[<img src="https://colab.research.google.com/assets/colab-badge.svg">](https://colab.research.google.com/github/facebookresearch/esm/blob/main/examples/protein-programming-language/tutorial.ipynb)

## Design programs

We have also provided some example programs as described in our original paper.

| Design task                                 | Figure in paper | Program file                                                                |
|:--------------------------------------------|:----------------|:----------------------------------------------------------------------------|
| Free hallucination                          | Figure 2A       | [free_hallucination.py](programs/free_hallucination.py)                     |
| Fixed backbone design                       | Figure 2D       | [fixed_backbone.py](programs/fixed_backbone.py)                             |
| Secondary structure design                  | Figure 2G       | [secondary_structure.py](programs/secondary_structure.py)                   |
| Functional site scaffolding                 | Figure 2H       | [functional_site_scaffolding.py](programs/functional_site_scaffolding.py)   |
| Symmetric monomer design                    | Figure 3A       | [symmetric_monomer.py](programs/symmetric_monomer.py)                       |
| Two-level symmetric<br>homo-oligomer design | Figure 4A       | [symmetric_two_level_multimer.py](programs/symmetric_two_level_multimer.py) |
| Symmetric binding site<br>scaffolding       | Figure 5A       | [symmetric_binding.py](programs/symmetric_binding.py)                       |

