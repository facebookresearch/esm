======================================================
Evolutionary Scale Modeling (esm): Pretrained language models for proteins
======================================================

This repository contains a PyTorch implementation of the transformer protein language models in
`"Biological structure and function emerge from scaling unsupervised learning
to 250 million protein sequences" (Rives et al., 2019)`__
from Facebook AI Research, along with pre-trained models.

__ https://doi.org/10.1101/622803

Quickstart
==========

As a prerequisite, you must have PyTorch 1.5 or later installed to use this repository.
A cuda device is optional and will be auto-detected.

Use this one-liner for installation:

.. code-block:: bash

    $ pip install git+https://github.com/facebookresearch/esm.git

Then, you can load and use a pretrained model as follows:

.. code-block:: python

    import torch
    import esm

    # Load 34 layer model
    model, alphabet = esm.pretrained.esm1_t34_670M_UR50S()
    batch_converter = alphabet.get_batch_converter()

    # Prepare data (two protein sequences)
    data = [("protein1", "MYLYQKIKN"), ("protein2", "MNAKYD")]
    batch_labels, batch_strs, batch_tokens = batch_converter(data)

    # Extract per-residue representations (on CPU)
    with torch.no_grad():
        results = model(batch_tokens, repr_layers=[34])
    token_representations = results["representations"][34]

    # Generate per-sequence representations via averaging
    # NOTE: token 0 is always a beginning-of-sequence token, so the first residue is token 1.
    sequence_representations = []
    for i, (_, seq) in enumerate(data):
        sequence_representations.append(token_representations[i, 1:len(seq) + 1].mean(0))


We also support PyTorch Hub, which removes the need to clone and/or install this repository yourself:

.. code-block:: python

    import torch

    model, alphabet = torch.hub.load("facebookresearch/esm", "esm1_t34_670M_UR50S")


FASTA representation extractor
------------------------------

For your convenience, we have provided a script that efficiently extracts representations in bulk from a FASTA file:

.. code-block:: bash

    # Extract final-layer representations for a FASTA file from a 34-layer model
    $ python extract.py esm1_t34_670M_UR50S examples/some_proteins.fasta my_reprs/ \
        --repr_layers 0 32 34 --include-per-tok --include-mean

    

    # my_reprs/ now contains one ".pt" file per FASTA sequence; use torch.load() to load them
    # extract.py has flags that determine what's included in the ".pt" file:
    # --repr-layers (default: final only) selects which layers to include representations from.
    # --include-per-tok includes the full sequence, with an embedding per amino acid (seq_len x hidden_dim).
    # --include-mean includes the embeddings per layer, averaged over the full sequence.
    # --include-bos includes the embeddings from the beginning-of-sequence token.

Available models
================

The following table lists the pretrained models available for use.
Names are self-explanatory corresponding to Table 1 in the updated paper 
(number of layers, number of params, training dataset).

* esm1_t34_670M_UR50S -- this is the best model and should be go-to.
* esm1_t34_670M_UR50D
* esm1_t34_670M_UR100
* esm1_t12_85M_UR50S
* esm1_t6_43M_UR50S

Reference
=========

If you find the model useful in your research, we ask that you cite the
following paper:

.. code-block:: bibtex

    @article{rives2019biological,
      author={Rives, Alexander and Meier, Joshua and Sercu, Tom and Goyal, Siddharth and Lin, Zeming and Guo, Demi and Ott, Myle and Zitnick, C. Lawrence and Ma, Jerry and Fergus, Rob},
      title={Biological Structure and Function Emerge from Scaling Unsupervised Learning to 250 Million Protein Sequences},
      year={2019},
      doi={10.1101/622803},
      url={https://www.biorxiv.org/content/10.1101/622803v3},
      journal={bioRxiv}
    }

Additionally, much of this code hails from the excellent `fairseq`__ sequence modeling framework; we have released this standalone model to facilitate more lightweight and flexible usage. We encourage those who wish to pretrain protein language models from scratch to use fairseq.

__ https://github.com/pytorch/fairseq

License
=======

This source code is licensed under the MIT license found in the ``LICENSE`` file
in the root directory of this source tree.
