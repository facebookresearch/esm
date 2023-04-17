# LM design examples

This folder contains code for demonstration of protein design using a language model. The code was used to perform the two design tasks specified at the paper [Language models generalize beyond natural proteins
](https://www.biorxiv.org/content/10.1101/2022.12.21.521521v1).


## Notebook examples

Refer to the two notebooks at this folder to run the fixed backbone and free generation design tasks.


## Shell examples

To run the two design tasks from shell, do the following:

1. First, install additional requirements: ```pip install -r additional_requirements.txt```
2. Running Fixed backbone design: ```python -m lm_design task=fixedbb pdb_fn=$PWD/2N2U.pdb```
3. Running Free generation design: ```python -m lm_design task=free_generation```

Notes:
Use the ```seed=<number>``` flag to generate different designs, e.g:
```python -m lm_design task=free_generation seed=42```

Control generated length in free generation using ```free_generation_length=<number>```, e.g:
```python -m lm_design task=free_generation free_generation_length=68```

Other, more advanced configurations can be observed at [config.yaml](conf/config.yaml)


## Paper data
The data from the preprint is available under [paper-data/](paper-data).
This includes designed sequences, their predicted structures, experimental validation results, linear projection for pairwise distance prediction, and details on dataset construction for model training.
