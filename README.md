# Evolutionary Scale Modeling


This repository contains code and pre-trained weights for **Transformer protein language models** from Facebook AI Research, including our state-of-the-art **ESM-2** and **MSA Transformer**, as well as **ESM-1v** for predicting variant effects and **ESM-IF1** for inverse folding.
Transformer protein language models were introduced in our paper, ["Biological structure and function emerge from scaling unsupervised learning to 250 million protein sequences" (Rives et al., 2019)](https://doi.org/10.1101/622803).

**ESM-2 outperforms all tested single-sequence protein language models across a range of structure prediction tasks.**
The MSA Transformer (ESM-MSA-1) can improve performance on some proteins by leveraging MSA information.

<details><summary>Citation</summary>

```bibtex
@article{lin2022language,
  title={Language models of protein sequences at the scale of evolution enable accurate structure prediction},
  author={Lin, Zeming and Akin, Halil and Rao, Roshan and Hie, Brian and Zhu, Zhongkai and Lu, Wenting and dos Santos Costa, Allan and Fazel-Zarandi, Maryam and Sercu, Tom and Candido, Sal and others},
  journal={bioRxiv},
  year={2022},
  publisher={Cold Spring Harbor Laboratory}
}
@article{rives2021biological,
  title={Biological structure and function emerge from scaling unsupervised learning to 250 million protein sequences},
  author={Rives, Alexander and Meier, Joshua and Sercu, Tom and Goyal, Siddharth and Lin, Zeming and Liu, Jason and Guo, Demi and Ott, Myle and Zitnick, C Lawrence and Ma, Jerry and others},
  journal={Proceedings of the National Academy of Sciences},
  volume={118},
  number={15},
  pages={e2016239118},
  year={2021},
  publisher={National Acad Sciences}
}
```
</details>

<details><summary>Table of contents</summary>

- [Main models you should use](#main-models)
- [Comparison to related works](#perf_related)
- [Usage](#usage)
  - [Quick Start](#quickstart)
  - [Compute embeddings in bulk from FASTA](#bulk_fasta)
  - [CPU offloading for inference with large models](#fsdp)
  - [Zero-shot variant prediction](#zs_variant)
  - [Inverse folding](#invf)
- [Notebooks](#notebooks)
- [Available Models and Datasets](#available)
  - [Pre-trained Models](#available-models)
  - [ESM Structural Split Dataset](#available-esmssd)
  - [Pre-training Dataset Split](#available-pretraining-split)
- [Citations](#citations)
- [License](#license)
</details>

<details><summary>What's New</summary>

- August 2022: ESM-2 - new SOTA Language Models released (see [Lin et al. 2022](https://www.biorxiv.org/content/10.1101/2022.07.20.500902v1))
- April 2022: New inverse folding model ESM-IF1 released, trained on CATH and UniRef50 predicted structures.
- August 2021: Added flexibility to tokenizer to allow for spaces and special tokens (like `<mask>`) in sequence.
- July 2021: New pre-trained model ESM-1v released, trained on UniRef90 (see [Meier et al. 2021](https://doi.org/10.1101/2021.07.09.450648)).
- July 2021: New MSA Transformer released, with a minor fix in the row positional embeddings (`ESM-MSA-1b`).
- Feb 2021: MSA Transformer added (see [Rao et al. 2021](https://www.biorxiv.org/content/10.1101/2021.02.12.430858v1)). Example usage in [notebook](#notebooks).
- Dec 2020: [Self-Attention Contacts](#notebooks) for all pre-trained models (see [Rao et al. 2020](https://doi.org/10.1101/2020.12.15.422761))
- Dec 2020: Added new pre-trained model [ESM-1b](#perf_related) (see [Rives et al. 2019](https://doi.org/10.1101/622803) Appendix B)
- Dec 2020: [ESM Structural Split Dataset](#available-esmssd) (see [Rives et al. 2019](https://doi.org/10.1101/622803) Appendix A.10)

</details>

## Main models you should use <a name="main-models"></a>

| Shorthand | `esm.pretrained.`           | Dataset | Description  |
|-----------|-----------------------------|---------|--------------|
| ESM-2    | `esm2_t36_3B_UR50D()` `esm2_t48_15B_UR50D()`       | UR50 (sample UR90)  | SOTA general-purpose protein language model. Can be used to predict structure, function and other protein properties directly from individual sequences. Released with [Lin et al. 2022](https://www.biorxiv.org/content/10.1101/2022.07.20.500902v1) (Aug 2022 update). |
| ESM-1b    | `esm1b_t33_650M_UR50S()`       | UR50  | SOTA general-purpose protein language model. Can be used to predict structure, function and other protein properties directly from individual sequences. Released with [Rives et al. 2019](https://doi.org/10.1101/622803) (Dec 2020 update). |
| ESM-MSA-1b| `esm_msa1b_t12_100M_UR50S()` |  UR50 + MSA  | MSA Transformer language model. Can be used to extract embeddings from an MSA. Enables SOTA inference of structure. Released with [Rao et al. 2021](https://www.biorxiv.org/content/10.1101/2021.02.12.430858v2) (ICML'21 version, June 2021).  |
| ESM-1v    | `esm1v_t33_650M_UR90S_1()` ... `esm1v_t33_650M_UR90S_5()`| UR90  | Language model specialized for prediction of variant effects. Enables SOTA zero-shot prediction of the functional effects of sequence variations. Same architecture as ESM-1b, but trained on UniRef90. Released with [Meier et al. 2021](https://doi.org/10.1101/2021.07.09.450648). |
| ESM-IF1  | `esm_if1_gvp4_t16_142M_UR50()` | CATH + UR50 | Inverse folding model. Can be used to design sequences for given structures, or to predict functional effects of sequence variation for given structures. Enables SOTA fixed backbone sequence design. Released with [Hsu et al. 2022](https://doi.org/10.1101/2022.04.10.487779). |

For a complete list of available models, with details and release notes, see [Pre-trained Models](#available-models).

## Comparison to related works <a name="perf_related"></a>
<!--
DO NOT EDIT THIS TABLE! This is the source of truth:
https://docs.google.com/spreadsheets/d/1RPvWF47rIMEr-Jg-SRCoGElHcwCl5d7RyEeSyPgp59A/edit#gid=0
exported via https://www.tablesgenerator.com/html_tables
-->

<table class="tg">
<thead>
  <tr>
    <th class="tg-0thz"><span style="font-weight:bold">Task</span></th>
    <th class="tg-j6zm" colspan="3"><span style="font-weight:bold">Unsupervised contact prediction</span></th>
    <th class="tg-j6zm" colspan="2"><span style="font-weight:bold">Structure Prediction</span></th>
  </tr>
</thead>
<tbody>
  <tr>
    <td class="tg-j6zm"><span style="font-weight:bold">Test set</span></td>
    <td class="tg-j6zm"><span style="font-weight:bold">Large valid</span></td>
    <td class="tg-j6zm"><span style="font-weight:bold">CASP14</span></td>
    <td class="tg-j6zm"><span style="font-weight:bold">CAMEO (Apr-Jun 2022)</span></td>
    <td class="tg-j6zm"><span style="font-weight:bold">CASP14</span></td>
    <td class="tg-j6zm"><span style="font-weight:bold">CAMEO (Apr-Jun 2022)</span></td>
  </tr>
  <tr>
    <td class="tg-7zrl">Gremlin (Potts)</td>
    <td class="tg-7zrl">39.3</td>
    <td class="tg-7zrl"></td>
    <td class="tg-7zrl"></td>
    <td class="tg-7zrl"></td>
    <td class="tg-7zrl"></td>
  </tr>
  <tr>
    <td class="tg-7zrl">TAPE</td>
    <td class="tg-7zrl">11.2</td>
    <td class="tg-7zrl"></td>
    <td class="tg-7zrl"></td>
    <td class="tg-7zrl"></td>
    <td class="tg-7zrl"></td>
  </tr>
  <tr>
    <td class="tg-7zrl">ProtBert-BFD</td>
    <td class="tg-7zrl">34.1</td>
    <td class="tg-7zrl"></td>
    <td class="tg-7zrl"></td>
    <td class="tg-7zrl"></td>
    <td class="tg-7zrl"></td>
  </tr>
  <tr>
    <td class="tg-7zrl">Prot-T5-XL-BFD</td>
    <td class="tg-7zrl">35.6</td>
    <td class="tg-7zrl"></td>
    <td class="tg-7zrl"></td>
    <td class="tg-2b7s">46.1</td>
    <td class="tg-2b7s">62.6</td>
  </tr>
  <tr>
    <td class="tg-7zrl">Prot-T5-XL-Ur50 (3B)</td>
    <td class="tg-7zrl">47.9</td>
    <td class="tg-7zrl"></td>
    <td class="tg-7zrl"></td>
    <td class="tg-2b7s">49.8</td>
    <td class="tg-2b7s">69.4</td>
  </tr>
  <tr>
    <td class="tg-7zrl">ESM-1</td>
    <td class="tg-7zrl">33.7</td>
    <td class="tg-7zrl"></td>
    <td class="tg-7zrl"></td>
    <td class="tg-7zrl"></td>
    <td class="tg-7zrl"></td>
  </tr>
  <tr>
    <td class="tg-7zrl">ESM-1b</td>
    <td class="tg-7zrl">41.1</td>
    <td class="tg-7zrl">24.4</td>
    <td class="tg-7zrl">39</td>
    <td class="tg-2b7s">41.6</td>
    <td class="tg-2b7s">64.5</td>
  </tr>
  <tr>
    <td class="tg-7zrl">ESM-1v</td>
    <td class="tg-7zrl">35.3</td>
    <td class="tg-7zrl"></td>
    <td class="tg-7zrl"></td>
    <td class="tg-7zrl"></td>
    <td class="tg-7zrl"></td>
  </tr>
  <tr>
    <td class="tg-7zrl">ESM-MSA-1b</td>
    <td class="tg-7zrl">57.4</td>
    <td class="tg-7zrl"></td>
    <td class="tg-7zrl"></td>
    <td class="tg-7zrl"></td>
    <td class="tg-7zrl"></td>
  </tr>
  <tr>
    <td class="tg-7zrl">ESM-2 (8M)</td>
    <td class="tg-7zrl">15.9</td>
    <td class="tg-7zrl">9.8</td>
    <td class="tg-7zrl">15.7</td>
    <td class="tg-2b7s">36.7</td>
    <td class="tg-2b7s">48.1</td>
  </tr>
  <tr>
    <td class="tg-7zrl">ESM-2 (35M)</td>
    <td class="tg-7zrl">28.8</td>
    <td class="tg-7zrl">16.4</td>
    <td class="tg-7zrl">28.4</td>
    <td class="tg-2b7s">41.4</td>
    <td class="tg-2b7s">56.4</td>
  </tr>
  <tr>
    <td class="tg-7zrl">ESM-2 (150M)</td>
    <td class="tg-7zrl">42.2</td>
    <td class="tg-7zrl">26.8</td>
    <td class="tg-7zrl">40.1</td>
    <td class="tg-2b7s">49.0</td>
    <td class="tg-2b7s">64.9</td>
  </tr>
  <tr>
    <td class="tg-7zrl">ESM-2 (700M)</td>
    <td class="tg-7zrl">50.1</td>
    <td class="tg-7zrl">32.5</td>
    <td class="tg-7zrl">47.6</td>
    <td class="tg-2b7s">51.3</td>
    <td class="tg-2b7s">70.1</td>
  </tr>
  <tr>
    <td class="tg-7zrl">ESM-2 (3B)</td>
    <td class="tg-7zrl">52.7</td>
    <td class="tg-7zrl">34.0</td>
    <td class="tg-7zrl">49.9</td>
    <td class="tg-2b7s">52.5</td>
    <td class="tg-2b7s">71.8</td>
  </tr>
  <tr>
    <td class="tg-7zrl">ESM-2 (15B)</td>
    <td class="tg-7zrl">54.5</td>
    <td class="tg-7zrl">37.0</td>
    <td class="tg-7zrl">51.7</td>
    <td class="tg-2b7s">55.4</td>
    <td class="tg-2b7s">72.1</td>
  </tr>
</tbody>
</table>

Comparison to related protein language models on structure prediction tasks.

* All contact numbers are the top-L,LR precision metric, where long range means sequence separation of at least 24 residues
* For unsupervised contact prediction, a sparse linear combination of the attention heads is used to directly predict protein contacts,
fitted with logistic regression on 20 structures.
For more details on the method, see [Rao et al. 2020](https://doi.org/10.1101/2020.12.15.422761).
* For structure prediction, a AlphaFold2 structure module is trained directly from the frozen language model embeddings.
For more details on the method, see [Lin et al. 2022](https://www.biorxiv.org/content/10.1101/2022.07.20.500902v1).
* Direct coupling analysis methods (Gremlin, mfDCA, Psicov) and ESM-MSA-1 use the [trRosetta MSAs](https://yanglab.nankai.edu.cn/trRosetta/benchmark/), while other methods predict from single sequence.


## Usage <a name="usage"></a>

### Quick Start <a name="quickstart"></a>

As a prerequisite, you must have PyTorch installed to use this repository.

You can use this one-liner for installation, using the latest release of esm:

```bash
pip install fair-esm  # latest release, OR:
pip install git+https://github.com/facebookresearch/esm.git  # bleeding edge, current repo main branch
```


We also support PyTorch Hub, which removes the need to clone and/or install this repository yourself:

```python
import torch
model, alphabet = torch.hub.load("facebookresearch/esm:main", "esm2_t33_650M_UR50D")
```

After pip install, you can load and use a pretrained model as follows:

```python
import torch
import esm

# Load ESM-2 model
model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
batch_converter = alphabet.get_batch_converter()
model.eval()  # disables dropout for deterministic results

# Prepare data (first 2 sequences from ESMStructuralSplitDataset superfamily / 4)
data = [
    ("protein1", "MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG"),
    ("protein2", "KALTARQQEVFDLIRDHISQTGMPPTRAEIAQRLGFRSPNAAEEHLKALARKGVIEIVSGASRGIRLLQEE"),
    ("protein2 with mask","KALTARQQEVFDLIRD<mask>ISQTGMPPTRAEIAQRLGFRSPNAAEEHLKALARKGVIEIVSGASRGIRLLQEE"),
    ("protein3",  "K A <mask> I S Q"),
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
The following command extracts the final-layer embedding for a FASTA file from the ESM-2 model:

```bash
python scripts/extract.py esm2_t33_650M_UR50D examples/data/some_proteins.fasta \
  examples/data/some_proteins_emb_esm2 --repr_layers 0 32 33 --include mean per_tok
```

Directory `some_proteins_emb_esm2/` now contains one `.pt` file per FASTA sequence; use `torch.load()` to load them.
`scripts/extract.py` has flags that determine what's included in the `.pt` file:
* `--repr-layers` (default: final only) selects which layers to include embeddings from.
* `--include` specifies what embeddings to save. You can use the following:
  * `per_tok` includes the full sequence, with an embedding per amino acid (seq_len x hidden_dim).
  * `mean` includes the embeddings averaged over the full sequence, per layer.
  * `bos` includes the embeddings from the beginning-of-sequence token.
  (NOTE: Don't use with the pre-trained models - we trained without bos-token supervision)

### CPU offloading for inference with large models <a name="fsdp"></a>
If you want to load very large models like 15B and/or do inference on long sequences on your machine, regular GPU inference may lead to OOM errors.
We show how to load the model with Fairscale's [Fully Sharded Data Parallel (FSDP)](https://fairscale.readthedocs.io/en/stable/api/nn/fsdp.html) and
use its CPU offloading feature.
This allows to do inference of large models on a single GPU.
Please check out `examples/esm2_infer_fairscale_fsdp_cpu_offloading.py` for more details.

### Zero-shot variant prediction <a name="zs_variant"></a>
See "[examples/variant-prediction/](examples/variant-prediction/)" for code and pre-trained weights for the ESM-1v models described in
[Language models enable zero-shot prediction of the effects of mutations on protein function. (Meier et al. 2021)](https://doi.org/10.1101/2021.07.09.450648).

Note that ESM-2 could be used for variant prediction as well, and is expected to have similar performance to ESM-1v.

### Inverse folding <a name="invf"></a>
See "[examples/inverse_folding/](examples/inverse_folding/)" for detailed user guide. The ESM-IF1 model is described as `GVPTransformer` in [Learning inverse folding from millions of predicted structures. (Hsu et al. 2022)](https://doi.org/10.1101/2022.04.10.487779).

We also provide a colab notebook for the sequence design and sequence scoring functionalities.

[<img src="https://colab.research.google.com/assets/colab-badge.svg">](https://colab.research.google.com/github/facebookresearch/esm/blob/main/examples/inverse_folding/notebook.ipynb)

The ESM-IF1 inverse folding model is built for predicting protein sequences
from their backbone atom coordinates. We provide scripts here 1) to sample sequence
designs for a given structure and 2) to score sequences for a given structure.

Trained with 12M protein structures predicted by AlphaFold2, the ESM-IF1
model consists of invariant geometric input processing layers followed by a
sequence-to-sequence transformer, and achieves 51% native sequence recovery on
structurally held-out backbones with 72% recovery for buried residues.
The model is also trained with span masking to tolerate missing backbone
coordinates and therefore can predict sequences for partially masked structures.

#### Sample sequence designs for a given structure
The environment setup is described in [this subsection of examples/inverse_folding](examples/inverse_folding#recommended-environment).

To sample sequences for a given structure in PDB or mmCIF format, use the
`sample_sequences.py` script. The input file can have either `.pdb` or
`.cif` as suffix.

For example, to sample 3 sequence designs for the golgi casein kinase structure
(PDB [5YH2](https://www.rcsb.org/structure/5yh2); [PDB Molecule of the Month
from January 2022](https://pdb101.rcsb.org/motm/265)), we can run the following
command from the esm root directory:
```
python examples/inverse_folding/sample_sequences.py examples/inverse_folding/data/5YH2.pdb \
  --chain C --temperature 1 --num-samples 3 --outpath examples/inverse_folding/output/sampled_sequences.fasta
```

The sampled sequences will be saved in a fasta format to the specified output file.

The temperature parameter controls the sharpness of the probability
distribution for sequence sampling. Higher sampling temperatures yield more
diverse sequences but likely with lower native sequence recovery.
The default sampling temperature is 1. To optimize for native sequence
recovery, we recommend sampling with low temperature such as 1e-6.

#### Scoring sequences
To score the conditional log-likelihoods for sequences conditioned on a given
structure, use the `score_log_likelihoods.py` script.

For example, to score the sequences in `examples/inverse_folding/data/5YH2_mutated_seqs.fasta`
according to the structure in `examples/inverse_folding/data/5YH2.pdb`, we can run
the following command from the esm root directory:
```
python examples/inverse_folding/score_log_likelihoods.py examples/inverse_folding/data/5YH2.pdb \ 
  examples/inverse_folding/data/5YH2_mutated_seqs.fasta --chain C \
  --outpath examples/inverse_folding/output/5YH2_mutated_seqs_scores.csv
```

The conditional log-likelihoods are saved in a csv format in the specified output path.
The output values are the average log-likelihoods averaged over all amino acids in a sequence.

For more information, see "[./examples/inverse_folding/](examples/inverse_folding/)" for detailed user guide.

## Notebooks <a name="notebooks"></a>

### Inverse folding - predicting or scoring sequences based on backbone structures

[<img src="https://colab.research.google.com/assets/colab-badge.svg">](https://colab.research.google.com/github/facebookresearch/esm/blob/main/examples/inverse_folding/notebook.ipynb)

The ESM-IF1 inverse folding model predicts protein sequences from their backbone atom coordinates, trained with 12M protein structures predicted by AlphaFold2.
This notetook guide you through examples of sampling sequences, calculating conditional log-likelihoods, and extracting encoder output as structure representation.

### Supervised variant prediction - training a classifier on the embeddings

[<img src="https://colab.research.google.com/assets/colab-badge.svg">](https://colab.research.google.com/github/facebookresearch/esm/blob/main/examples/sup_variant_prediction.ipynb)


To help you get started with using the embeddings, this [jupyter notebook tutorial](examples/sup_variant_prediction.ipynb) shows how to train a supervised variant predictor using embeddings from ESM-1.
You can adopt a similar protocol to train a model for any downstream task, even with limited data.
First you can obtain the embeddings for ``examples/data/P62593.fasta`` either by [downloading the precomputed](https://dl.fbaipublicfiles.com/fair-esm/examples/P62593_reprs.tar.gz) embeddings
as instructed in the notebook or by running the following:

```bash
# Obtain the embeddings
python scripts/extract.py esm1v_t33_650M_UR90S_1 examples/data/P62593.fasta \
  examples/data/P62593_emb_esm1v --repr_layers 33 --include mean
```

Then, follow the remaining instructions in the tutorial. You can also run the tutorial in a [colab notebook](https://colab.research.google.com/github/facebookresearch/esm/blob/main/examples/sup_variant_prediction.ipynb).

**Note, alternatively use [the newer instructions for zero-shot variant prediction](examples/variant-prediction/),
which predicts mutational effects without any supervised training.**


### Unsupervised contact prediction
[<img src="https://colab.research.google.com/assets/colab-badge.svg">](https://colab.research.google.com/github/facebookresearch/esm/blob/main/examples/contact_prediction.ipynb)

This [jupyter notebook tutorial](examples/contact_prediction.ipynb) demonstrates contact prediction with both the ESM-2 and MSA Transformer (ESM-MSA-1) models.
Contact prediction is based on a logistic regression over the model's attention maps.
This methodology is based on our ICLR 2021 paper,
[Transformer protein language models are unsupervised structure learners. (Rao et al. 2020)](https://doi.org/10.1101/2020.12.15.422761)
The MSA Transformer (ESM-MSA-1) takes a multiple sequence alignment (MSA) as input, and uses the tied row self-attention maps in the same way.
See [MSA Transformer. (Rao et al. 2021)](https://www.biorxiv.org/content/10.1101/2021.02.12.430858v1).

To get unsupervised attention-based contacts, call `model.predict_contacts(tokens)` or `model(tokens, return_contacts=True)`.


### ESMStructuralSplitDataset and self-attention contact prediction
[<img src="https://colab.research.google.com/assets/colab-badge.svg">](https://colab.research.google.com/github/facebookresearch/esm/blob/main/examples/esm_structural_dataset.ipynb)

And this [jupyter notebook tutorial](examples/esm_structural_dataset.ipynb) shows how to load and index the `ESMStructuralSplitDataset`,
and computes the self-attention map unsupervised contact predictions using ESM-2.


## Available Models and Datasets <a name="available"></a>

### Pre-trained Models <a name="available-models"></a>

| Shorthand | `esm.pretrained.`           | #layers | #params | Dataset | Embedding Dim |  Model URL (automatically downloaded to `~/.cache/torch/hub/checkpoints`) |
|-----------|---------------------|---------|---------|---------|---------------|-----------------------------------------------------------------------|
| ESM-2     | `esm2_t48_15B_UR50D`         | 48 | 15B     | UR50/D 2021_04 | 5120 |  https://dl.fbaipublicfiles.com/fair-esm/models/esm2_t48_15B_UR50D.pt |
|           | `esm2_t36_3B_UR50D`          | 36 | 3B      | UR50/D 2021_04 | 2560 |  https://dl.fbaipublicfiles.com/fair-esm/models/esm2_t36_3B_UR50D.pt |
|           | `esm2_t33_650M_UR50D`        | 33 | 650M    | UR50/D 2021_04 | 1280 |  https://dl.fbaipublicfiles.com/fair-esm/models/esm2_t33_650M_UR50D.pt |
|           | `esm2_t30_150M_UR50D`        | 30 | 150M    | UR50/D 2021_04 | 640  |  https://dl.fbaipublicfiles.com/fair-esm/models/esm2_t30_150M_UR50D.pt |
|           | `esm2_t12_35M_UR50D`         | 12 | 35M     | UR50/D 2021_04 | 480  |  https://dl.fbaipublicfiles.com/fair-esm/models/esm2_t12_35M_UR50D.pt |
|           | `esm2_t6_8M_UR50D`           | 6  | 8M      | UR50/D 2021_04 | 320  |  https://dl.fbaipublicfiles.com/fair-esm/models/esm2_t6_8M_UR50D.pt |
| ESM-IF1   | `esm_if1_gvp4_t16_142M_UR50` | 20 | 124M    | CATH 4.3 + predicted structures for UR50 | 512  | https://dl.fbaipublicfiles.com/fair-esm/models/esm_if1_gvp4_t16_142M_UR50.pt   |
| ESM-1v    | `esm1v_t33_650M_UR90S_[1-5]` | 33 | 650M    | UR90/S 2020_03  | 1280 | https://dl.fbaipublicfiles.com/fair-esm/models/esm1v_t33_650M_UR90S_1.pt   |
| ESM-MSA-1b| `esm_msa1b_t12_100M_UR50S`   | 12 | 100M    | UR50/S + MSA 2018_03 | 768  | https://dl.fbaipublicfiles.com/fair-esm/models/esm_msa1b_t12_100M_UR50S.pt   |
| ESM-MSA-1 | `esm_msa1_t12_100M_UR50S`    | 12 | 100M    | UR50/S + MSA 2018_03 | 768  | https://dl.fbaipublicfiles.com/fair-esm/models/esm_msa1_t12_100M_UR50S.pt   |
| ESM-1b    | `esm1b_t33_650M_UR50S`       | 33 | 650M    | UR50/S 2018_03 | 1280 | https://dl.fbaipublicfiles.com/fair-esm/models/esm1b_t33_650M_UR50S.pt   |
| ESM-1     | `esm1_t34_670M_UR50S`        | 34 | 670M    | UR50/S 2018_03 | 1280 |  https://dl.fbaipublicfiles.com/fair-esm/models/esm1_t34_670M_UR50S.pt |
|           | `esm1_t34_670M_UR50D`        | 34 | 670M    | UR50/D 2018_03 | 1280 |  https://dl.fbaipublicfiles.com/fair-esm/models/esm1_t34_670M_UR50D.pt |
|           | `esm1_t34_670M_UR100`        | 34 | 670M    | UR100 2018_03  | 1280 |  https://dl.fbaipublicfiles.com/fair-esm/models/esm1_t34_670M_UR100.pt |
|           | `esm1_t12_85M_UR50S`         | 12 | 85M     | UR50/S 2018_03 | 768  |  https://dl.fbaipublicfiles.com/fair-esm/models/esm1_t12_85M_UR50S.pt  |
|           | `esm1_t6_43M_UR50S`          | 6  | 43M     | UR50/S 2018_03 | 768  |  https://dl.fbaipublicfiles.com/fair-esm/models/esm1_t6_43M_UR50S.pt   |


Here is a chronological list of the released models and the paper they were introduced in:

| Shorthand  | Release Notes |
|------------|---------------|
| ESM-1      | Released with Rives et al. 2019 (Aug 2020 update). |
| ESM-1b     | Released with Rives et al. 2019 (Dec 2020 update). See Appendix B. |
| ESM-MSA-1  | Released with Rao et al. 2021 (Preprint v1). |
| ESM-MSA-1b | Released with Rao et al. 2021 (ICML'21 version, June 2021). |
| ESM-1v     | Released with Meier et al. 2021. |
| ESM-IF1    | Released with Hsu et al. 2022. |
| ESM-2      | Released with Lin et al. 2022. |

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

### Pre-training Dataset Split  <a name="available-pretraining-split"></a>
The split files establishing which UniRef50 clusters were used as held-out evaluation set for pre-training
in [Rives et al. 2019](https://doi.org/10.1101/622803) and [Rao et al. 2021](https://doi.org/10.1101/2021.02.12.430858) can be found here:
* [UniRef50 IDs of evaluation set](https://dl.fbaipublicfiles.com/fair-esm/pretraining-data/uniref201803_ur50_valid_headers.txt.gz): 3.016 M clusters
* [UniRef100 IDs of evaluation set](https://dl.fbaipublicfiles.com/fair-esm/pretraining-data/uniref201803_ur100_valid_headers.txt.gz): 13.745 M proteins, expanding the same UniRef50 clusters.

These files only contain only the UniRef50 IDs and UniRef100 IDs corresponding to the [UniRef database, 2018-03 release](https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2018_03/uniref/)
which is released by the UniProt Consortium under a [Creative Commons Attribution (CC BY 4.0) License](https://www.uniprot.org/help/license).

## Citations <a name="citations"></a>

If you find the models useful in your research, we ask that you cite the relevant paper:

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

For the self-attention contact prediction:

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

For the MSA Transformer:

```bibtex
@article{rao2021msa,
  author = {Rao, Roshan and Liu, Jason and Verkuil, Robert and Meier, Joshua and Canny, John F. and Abbeel, Pieter and Sercu, Tom and Rives, Alexander},
  title={MSA Transformer},
  year={2021},
  doi={10.1101/2021.02.12.430858},
  url={https://www.biorxiv.org/content/10.1101/2021.02.12.430858v1},
  journal={bioRxiv}
}
```

For variant prediction using ESM-1v:

```bibtex
@article{meier2021language,
  author = {Meier, Joshua and Rao, Roshan and Verkuil, Robert and Liu, Jason and Sercu, Tom and Rives, Alexander},
  title = {Language models enable zero-shot prediction of the effects of mutations on protein function},
  year={2021},
  doi={10.1101/2021.07.09.450648},
  url={https://www.biorxiv.org/content/10.1101/2021.07.09.450648v1},
  journal={bioRxiv}
}
```

For inverse folding using ESM-IF1:

```bibtex
@article{hsu2022learning,
	author = {Hsu, Chloe and Verkuil, Robert and Liu, Jason and Lin, Zeming and Hie, Brian and Sercu, Tom and Lerer, Adam and Rives, Alexander},
	title = {Learning inverse folding from millions of predicted structures},
	year = {2022},
	doi = {10.1101/2022.04.10.487779},
	url = {https://www.biorxiv.org/content/early/2022/04/10/2022.04.10.487779},
	journal = {bioRxiv}
}

```

Much of this code builds on the [fairseq](https://github.com/pytorch/fairseq) sequence modeling framework. We use fairseq internally for our protein language modeling research. We highly recommend trying it out if you'd like to pre-train protein language models from scratch.

Additionally, if you would like to use the variant prediction benchmark from Meier et al. (2021), we provide a bibtex file with citations for all data in [./examples/variant-prediction/mutation_data.bib](./examples/variant-prediction/mutation_data.bib). You can cite each paper individually, or add all citations in bulk using the LaTeX command:

```tex
\nocite{wrenbeck2017deep,klesmith2015comprehensive,haddox2018mapping,romero2015dissecting,firnberg2014comprehensive,deng2012deep,stiffler2015evolvability,jacquier2013capturing,findlay2018comprehensive,mclaughlin2012spatial,kitzman2015massively,doud2016accurate,pokusaeva2019experimental,mishra2016systematic,kelsic2016rna,melnikov2014comprehensive,brenan2016phenotypic,rockah2015systematic,wu2015functional,aakre2015evolving,qi2014quantitative,matreyek2018multiplex,bandaru2017deconstruction,roscoe2013analyses,roscoe2014systematic,mavor2016determination,chan2017correlation,melamed2013deep,starita2013activity,araya2012fundamental}
```

## License <a name="license"></a>

This source code is licensed under the MIT license found in the `LICENSE` file
in the root directory of this source tree.
