# Evolutionary Scale Modeling


This repository contains code and pre-trained weights for **Transformer protein language models** from Facebook AI Research, including our state-of-the-art **ESM-1b protein language model**.
The models are described in detail in our paper, ["Biological structure and function emerge from scaling unsupervised learning to 250 million protein sequences" (Rives et al., 2019)](https://doi.org/10.1101/622803),
which first proposed protein language modeling with Transformers.

**ESM-1b outperforms all tested protein language models across a range of structure prediction tasks.**

<details><summary>Citation</summary>

```bibtex
@article{rives2019biological,
  author={Rives, Alexander and Meier, Joshua and Sercu, Tom and Goyal, Siddharth and Lin, Zeming and Liu, Jason and Guo, Demi and Ott, Myle and Zitnick, C. Lawrence and Ma, Jerry and Fergus, Rob},
  title={Biological Structure and Function Emerge from Scaling Unsupervised Learning to 250 Million Protein Sequences},
  year={2019},
  doi={10.1101/622803},
  url={https://www.biorxiv.org/content/10.1101/622803v4},
  journal={bioRxiv}
}
```
</details>

<details><summary>Table of contents</summary>
  
- [Comparison to related works](#perf-related)
- [Usage](#usage)
  - [Quick Start](#quickstart)
  - [Compute embeddings in bulk from FASTA](#bulk-fasta)
  - [Notebooks](#notebooks)
- [Benchmarks](#perf)
  - [Comparison on several tasks](#perf-related)
- [Available Models and Datasets](#available)
  - [Pre-trained Models](#available-models)
  - [ESM Structural Split Dataset](#available-esmssd)
- [Citations](#citations)
- [License](#license)
</details>

<details><summary>What's New</summary>
  
- Dec 2020: [Self-Attention Contacts](#quickstart) for all pre-trained models (see [Rao et al. 2020](https://www.biorxiv.org/content/10.1101/2020.12.15.422761v1))
- Dec 2020: Added new pre-trained model [ESM-1b](#perf-related) (see [Rives et al. 2019](https://doi.org/10.1101/622803) Appendix B)
- Dec 2020: [ESM Structural Split Dataset](#available-esmssd) (see [Rives et al. 2019](https://doi.org/10.1101/622803) Appendix A.10)
  
</details>

## Comparison to related works <a name="perf-related"></a>

### Supervised downstreams

| Model                                                       | Pre-training | Params | SSP  | Contact |
|-------------------------------------------------------------|--------------|--------|------|---------|
| [UniRep](https://www.nature.com/articles/s41592-019-0598-1) | UR50\*       | 18M    |  58.4 | 21.9    |
| [SeqVec](https://github.com/rostlab/SeqVec)                 | UR50\*       | 93M    |  62.1 | 29.0    |
| [TAPE](https://github.com/songlab-cal/tape)                 | PFAM\*       | 38M    |  58.0 | 23.2    |
| [ProtBert-BFD](https://github.com/agemagician/ProtTrans)    | BFD\*        | 420M   |  70.0 | 50.3    | 
| [Prot-T5-XL-BFD](https://github.com/agemagician/ProtTrans)  | BFD\*        | 3B     |  71.4 | 55.9    |
| LSTM biLM (S)                                               | UR50/S       | 28M    |  60.4 | 24.1    |
| LSTM biLM (L)                                               | UR50/S       | 113M   |  62.4 | 27.8    |
| Transformer-6                                               | UR50/S       | 43M    |  62.0 | 30.2    |
| Transformer-12                                              | UR50/S       | 85M    |  65.4 | 37.7    |
| Transformer-34                                              | UR100        | 670M   |  64.3 | 32.7    |
| Transformer-34                                              | UR50/S       | 670M   |  69.2 | 50.2    |
| **ESM-1b**                                                  | **UR50/S**    | **650M** | **71.6**|**56.9** |

Comparison to related protein language models.
(SSP) Secondary structure Q8 accuracy on CB513.
(Contact) Top-L long range contact precision on RaptorX test set.

\* Pre-training datasets from related works have differences from ours.

### Unsupervised downstreams
| Model                                                       | Pre-training | Params | L     | L/5     |
|-------------------------------------------------------------|--------------|--------|-------|---------|
| [mfDCA](https://www.pnas.org/content/108/49/E1293)&dagger;         |              |        |  33.0 | 54.2    |     
| [Psicov](https://academic.oup.com/bioinformatics/article/28/2/184/198108)&dagger;| |        |  32.6 | 58.1    |     
| [Gremlin](https://github.com/nickbhat/mogwai)&dagger;               |              |        |  39.3 | 62.8    |     
| [TAPE](https://github.com/songlab-cal/tape)                 | PFAM\*       | 38M    |  11.2 | 17.9    |
| [ProtBert-BFD](https://github.com/agemagician/ProtTrans)    | BFD\*        | 420M   |  34.1 | 57.4    | 
| [Prot-T5-XL-BFD](https://github.com/agemagician/ProtTrans)  | BFD\*        | 3B     |  35.6 | 57.8    |
| Transformer-6                                               | UR50/S       | 43M    |  13.2 | 21.5    |
| Transformer-12                                              | UR50/S       | 85M    |  23.7 | 39.3    |
| Transformer-34                                              | UR50/S       | 670M   |  34.7 | 56.0    |
| **ESM-1b**                                                  | **UR50/S**   | **650M**   | **41.1**|**66.1** |

Comparison to related protein language models. Average Top-L and Top-L/5 long range contact precision on 14842 test structures for Transformer models trained on 20 structures. &dagger; Direct coupling analysis methods (Gremlin, mfDCA, Psicov) use [trRosetta MSAs](https://yanglab.nankai.edu.cn/trRosetta/benchmark/) which include sequences from metagenomics. 

## Usage <a name="usage"></a>

### Quick Start <a name="quickstart"></a>

As a prerequisite, you must have PyTorch 1.5 or later installed to use this repository.

You can either work in the root of this repository, or use this one-liner for installation:

```bash
$ pip install git+https://github.com/facebookresearch/esm.git
```

We also support PyTorch Hub, which removes the need to clone and/or install this repository yourself:

```python
import torch
model, alphabet = torch.hub.load("facebookresearch/esm", "esm1b_t33_650M_UR50S")
```

Then, you can load and use a pretrained model as follows:

```python
import torch
import esm

# Load ESM-1b model
model, alphabet = esm.pretrained.esm1b_t33_650M_UR50S()
batch_converter = alphabet.get_batch_converter()

# Prepare data (first 2 sequences from ESMStructuralSplitDataset superfamily / 4)
data = [
    ("protein1", "MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG"),
    ("protein2", "KALTARQQEVFDLIRDHISQTGMPPTRAEIAQRLGFRSPNAAEEHLKALARKGVIEIVSGASRGIRLLQEE"),
]
batch_labels, batch_strs, batch_tokens = batch_converter(data)

# Extract per-residue representations (on CPU)
with torch.no_grad():
    results = model(batch_tokens, repr_layers=[33], return_contacts=True)
token_representations = results["representations"][33]

# Generate per-sequence representations via averaging
# NOTE: token 0 is always a beginning-of-sequence token, so the first residue is token 1.
sequence_representations = []
for i, (_, seq) in enumerate(data):
    sequence_representations.append(token_representations[i, 1 : len(seq) + 1].mean(0))

# Look at the unsupervised self-attention map contact predictions
import matplotlib.pyplot as plt
for (_, seq), attention_contacts in zip(data, results["contacts"]):
    plt.matshow(attention_contacts[: len(seq), : len(seq)])
    plt.title(seq)
    plt.show()
```

### Compute embeddings in bulk from FASTA <a name="bulk_fasta"></a>

We provide a script that efficiently extracts embeddings in bulk from a FASTA file.
A cuda device is optional and will be auto-detected.
The following command extracts the final-layer embedding for a FASTA file from the ESM-1b model:

```bash
$ python extract.py esm1b_t33_650M_UR50S examples/some_proteins.fasta my_reprs/ \
    --repr_layers 0 32 33 --include mean per_tok
```

Directory `my_reprs/` now contains one `.pt` file per FASTA sequence; use `torch.load()` to load them.
`extract.py` has flags that determine what's included in the `.pt` file:
* `--repr-layers` (default: final only) selects which layers to include embeddings from.
* `--include` specifies what embeddings to save. You can use the following:
  * `per_tok` includes the full sequence, with an embedding per amino acid (seq_len x hidden_dim).
  * `mean` includes the embeddings averaged over the full sequence, per layer.
  * `bos` includes the embeddings from the beginning-of-sequence token. 
  (NOTE: Don't use with the pre-trained models - we trained without bos-token supervision)

### Notebooks <a name="notebooks"></a> 

#### Variant prediction - using the embeddings

[<img src="https://colab.research.google.com/assets/colab-badge.svg">](https://colab.research.google.com/github/facebookresearch/esm/blob/master/examples/variant_prediction.ipynb)


To help you get started with using the embeddings, this [jupyter notebook tutorial](examples/variant_prediction.ipynb) shows how to train a variant predictor using embeddings from ESM-1.
You can adopt a similar protocol to train a model for any downstream task, even with limited data.
First you can obtain the embeddings for ``examples/P62593.fasta`` either by [downloading the precomputed](https://dl.fbaipublicfiles.com/fair-esm/examples/P62593_reprs.tar.gz) embeddings
as instructed in the notebook or by running the following:

```bash
# Obtain the embeddings
$ python extract.py esm1_t34_670M_UR50S examples/P62593.fasta examples/P62593_reprs/ \
    --repr_layers 34 --include mean
```

Then, follow the remaining instructions in the tutorial. You can also run the tutorial in a [colab notebook](https://colab.research.google.com/github/facebookresearch/esm/blob/master/examples/variant_prediction.ipynb).

#### ESMStructuralSplitDataset and self-attention contact prediction
[<img src="https://colab.research.google.com/assets/colab-badge.svg">](https://colab.research.google.com/github/facebookresearch/esm/blob/master/examples/esm_structural_dataset.ipynb)

And this [jupyter notebook tutorial](examples/esm_structural_dataset.ipynb) shows how to load and index the `ESMStructuralSplitDataset`,
and computes the self-attention map contact predictions as described in our paper "Transformer protein language models are unsupervised structure learners".



## Available Models and Datasets <a name="available"></a>

### Pre-trained Models <a name="available-models"></a>

| Shorthand | Full Name           | #layers | #params | Dataset | Embedding Dim |  Model URL                                                             |
|-----------|---------------------|---------|---------|---------|---------------|-----------------------------------------------------------------------|
| ESM-1b    | esm1b_t33_650M_UR50S | 33     | 650M    | UR50/S  | 1280          | https://dl.fbaipublicfiles.com/fair-esm/models/esm1b_t33_650M_UR50S.pt   |
| ESM1-main | esm1_t34_670M_UR50S | 34      | 670M    | UR50/S  | 1280          |  https://dl.fbaipublicfiles.com/fair-esm/models/esm1_t34_670M_UR50S.pt |
|           | esm1_t34_670M_UR50D | 34      | 670M    | UR50/D  | 1280          |  https://dl.fbaipublicfiles.com/fair-esm/models/esm1_t34_670M_UR50D.pt |
|           | esm1_t34_670M_UR100 | 34      | 670M    | UR100   | 1280          |  https://dl.fbaipublicfiles.com/fair-esm/models/esm1_t34_670M_UR100.pt |
|           | esm1_t12_85M_UR50S  | 12      | 85M     | UR50/S  | 768           |  https://dl.fbaipublicfiles.com/fair-esm/models/esm1_t12_85M_UR50S.pt  |
|           | esm1_t6_43M_UR50S   | 6       | 43M     | UR50/S  | 768           |  https://dl.fbaipublicfiles.com/fair-esm/models/esm1_t6_43M_UR50S.pt   |


### ESM Structural Split Dataset <a name="available-esmssd"></a>
This is a five-fold cross validation dataset of protein domain structures that can be used to measure generalization of representations
across different levels of structural dissimilarity.
The dataset implements structural holdouts at the family, superfamily, and fold
level. The SCOPe database is used to classify domains. Independently for each level of structural hold-out,
the domains are split into 5 equal sets, i.e. five sets of folds, superfamilies, or families. This ensures
that for each of the five partitions, structures having the same classification do not appear in both the
train and test sets. For a given classification level each structure appears in a test set once, so that
in the cross validation experiment each of the structures will be evaluated exactly once.

The dataset provides 3d coordinates, distance maps, and secondary structure labels.
For further details on the construction of the dataset 
see [Rives et al. 2019](https://doi.org/10.1101/622803) Appendix A.10.

This [jupyter notebook tutorial](examples/esm_structural_dataset.ipynb) shows how to load and index the `ESMStructuralSplitDataset`.

`ESMStructuralSplitDataset`, upon initializing, will download `splits` and `pkl`. 
We also provide `msas` for each of the domains. The data can be directly downloaded below.

| Name   | Description                                                                   | URL                                                                   |
|--------|-------------------------------------------------------------------------------|-----------------------------------------------------------------------|
| splits | train/valid splits                                                            | https://dl.fbaipublicfiles.com/fair-esm/structural-data/splits.tar.gz |
| pkl    | pkl objects containing sequence, SSP labels, distance map, and 3d coordinates | https://dl.fbaipublicfiles.com/fair-esm/structural-data/pkl.tar.gz    |
| msas   | a3m files containing MSA for each domain                                      | https://dl.fbaipublicfiles.com/fair-esm/structural-data/msas.tar.gz   | 

## Citations <a name="citations"></a>

If you find the models useful in your research, we ask that you cite the
following paper:

```bibtex
@article{rives2019biological,
  author={Rives, Alexander and Meier, Joshua and Sercu, Tom and Goyal, Siddharth and Lin, Zeming and Liu, Jason and Guo, Demi and Ott, Myle and Zitnick, C. Lawrence and Ma, Jerry and Fergus, Rob},
  title={Biological Structure and Function Emerge from Scaling Unsupervised Learning to 250 Million Protein Sequences},
  year={2019},
  doi={10.1101/622803},
  url={https://www.biorxiv.org/content/10.1101/622803v4},
  journal={bioRxiv}
}
```
 
For the self-attention contact prediction, see [the following paper (biorxiv preprint)](https://www.biorxiv.org/content/10.1101/2020.12.15.422761v1):

```bibtex
@article{rao2020transformer,
  author = {Rao, Roshan M and Meier, Joshua and Sercu, Tom and Ovchinnikov, Sergey and Rives, Alexander},
  title={Transformer protein language models are unsupervised structure learners},
  year={2020},
  doi={10.1101/2020.12.15.422761},
  url={https://www.biorxiv.org/content/10.1101/2020.12.15.422761v1},
  journal={bioRxiv}
}
```


Much of this code builds on the [fairseq](https://github.com/pytorch/fairseq) sequence modeling framework. We use fairseq internally for our protein language modeling research. We highly recommend trying it out if you'd like to pre-train protein language models from scratch.

## License <a name="license"></a>

This source code is licensed under the MIT license found in the `LICENSE` file
in the root directory of this source tree.
