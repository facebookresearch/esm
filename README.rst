======================================================
Evolutionary Scale Modeling (ESM)
======================================================
Pretrained language models for proteins
--------------

This repository contains a PyTorch implementation of and pre-trained weights for the transformer protein language models in
`"Biological structure and function emerge from scaling unsupervised learning to 250 million protein sequences" (Rives et al., 2019)`_
from Facebook AI Research:

.. code-block:: bibtex

    @article{rives2019biological,
      author={Rives, Alexander and Meier, Joshua and Sercu, Tom and Goyal, Siddharth and Lin, Zeming and Guo, Demi and Ott, Myle and Zitnick, C. Lawrence and Ma, Jerry and Fergus, Rob},
      title={Biological Structure and Function Emerge from Scaling Unsupervised Learning to 250 Million Protein Sequences},
      year={2019},
      doi={10.1101/622803},
      url={https://www.biorxiv.org/content/10.1101/622803v3},
      journal={bioRxiv}
    }


.. _"Biological structure and function emerge from scaling unsupervised learning to 250 million protein sequences" (Rives et al., 2019): https://doi.org/10.1101/622803

Quickstart
==========

As a prerequisite, you must have PyTorch 1.5 or later installed to use this repository.
A cuda device is optional and will be auto-detected.

You can either work in the root of this repository, or use this one-liner for installation:

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

    # Extract per-residue embeddings (on CPU)
    with torch.no_grad():
        results = model(batch_tokens, repr_layers=[34])
    token_embeddings = results["representations"][34]

    # Generate per-sequence embeddings via averaging
    # NOTE: token 0 is always a beginning-of-sequence token, so the first residue is token 1.
    sequence_embeddings = []
    for i, (_, seq) in enumerate(data):
        sequence_embeddings.append(token_embeddings[i, 1:len(seq) + 1].mean(0))


We also support PyTorch Hub, which removes the need to clone and/or install this repository yourself:

.. code-block:: python

    import torch

    model, alphabet = torch.hub.load("facebookresearch/esm", "esm1_t34_670M_UR50S")


FASTA embedding extractor
================

For your convenience, we have provided a script that efficiently extracts embeddings in bulk from a FASTA file:

.. code-block:: bash

    # Extract final-layer embedding for a FASTA file from a 34-layer model
    $ python extract.py esm1_t34_670M_UR50S examples/some_proteins.fasta my_reprs/ \
        --repr_layers 0 32 34 --include mean per_tok

    

    # my_reprs/ now contains one ".pt" file per FASTA sequence; use torch.load() to load them
    # extract.py has flags that determine what's included in the ".pt" file:
    # --repr-layers (default: final only) selects which layers to include embeddings from.
    # --include specifies what embeddings to save. You can use the following:
    # * per_tok includes the full sequence, with an embedding per amino acid (seq_len x hidden_dim).
    # * mean includes the embeddings averaged over the full sequence, per layer.
    # * bos includes the embeddings from the beginning-of-sequence token. 
    #    (NOTE: Don't use with the pre-trained models - we trained without bos-token supervision)

Tutorial
================

|ImageLink|_

.. |ImageLink| image:: https://colab.research.google.com/assets/colab-badge.svg
.. _ImageLink: https://colab.research.google.com/github/facebookresearch/esm/blob/master/examples/variant_prediction.ipynb


To help you get started, we `provide a jupyter notebook tutorial`__ demonstrating how to train a variant predictor using embeddings from ESM. You can adopt a similar protocol to train a model for any downstream task, even with limited data.
First you can obtain the embeddings for ``examples/P62593.fasta`` either by `downloading the precomputed`__ embeddings
as instructed in the notebook or by running the following:

.. code-block:: bash

    # Obtain the embeddings
    $ python extract.py esm1_t34_670M_UR50S examples/P62593.fasta examples/P62593_reprs/ \
        --repr_layers 34 --include mean
__ examples/variant_prediction.ipynb
__ https://dl.fbaipublicfiles.com/fair-esm/examples/P62593_reprs.tar.gz

Then, follow the remaining instructions in the tutorial. You can also run the tutorial in a `colab notebook`__.

__ https://colab.research.google.com/github/facebookresearch/esm/blob/master/examples/variant_prediction.ipynb


Available models
================

The following table lists the pretrained models available for use.
See also Table 1 in `the paper`_.

+-----------+---------------------+---------+---------+---------+---------------+----------------+-----------------------------------------------------------------------+
| Shorthand | Full Name           | #layers | #params | Dataset | Embedding Dim | Perplexity/ECE | Model URL                                                             |
+-----------+---------------------+---------+---------+---------+---------------+----------------+-----------------------------------------------------------------------+
| ESM1-main | esm1_t34_670M_UR50S | 34      | 670M    | UR50/S  | 1280          | 8.54           | https://dl.fbaipublicfiles.com/fair-esm/models/esm1_t34_670M_UR50S.pt |
+-----------+---------------------+---------+---------+---------+---------------+----------------+-----------------------------------------------------------------------+
|           | esm1_t34_670M_UR50D | 34      | 670M    | UR50/D  | 1280          | 8.46           | https://dl.fbaipublicfiles.com/fair-esm/models/esm1_t34_670M_UR50D.pt |
+-----------+---------------------+---------+---------+---------+---------------+----------------+-----------------------------------------------------------------------+
|           | esm1_t34_670M_UR100 | 34      | 670M    | UR100   | 1280          | 10.32          | https://dl.fbaipublicfiles.com/fair-esm/models/esm1_t34_670M_UR100.pt |
+-----------+---------------------+---------+---------+---------+---------------+----------------+-----------------------------------------------------------------------+
|           | esm1_t12_85M_UR50S  | 12      | 85M     | UR50/S  | 768           | 10.45          | https://dl.fbaipublicfiles.com/fair-esm/models/esm1_t12_85M_UR50S.pt  |
+-----------+---------------------+---------+---------+---------+---------------+----------------+-----------------------------------------------------------------------+
|           | esm1_t6_43M_UR50S   | 6       | 43M     | UR50/S  | 768           | 11.79          | https://dl.fbaipublicfiles.com/fair-esm/models/esm1_t6_43M_UR50S.pt   |
+-----------+---------------------+---------+---------+---------+---------------+----------------+-----------------------------------------------------------------------+


Comparison to related work
================
This table compares to related pre-training methods, and corresponds to Table 8 in `the paper`_.
The last 3 columns are the major benchmark results:

* RH: Remote Homology at the fold level, using Hit-10 metric on SCOP.
* SSP: Secondary structure Q8 accuracy on CB513. 
* Contact: Top-L long range contact precision on RaptorX test set from `Wang et al. (2017)`_.

.. _the paper: https://doi.org/10.1101/622803

+----------------+--------------+--------+------+------+---------+
| Model          | Pre-training | Params | RH   | SSP  | Contact |
+----------------+--------------+--------+------+------+---------+
| `UniRep`_      |              | 18M    | .527 | 58.4 | 21.9    |
+----------------+--------------+--------+------+------+---------+
| `SeqVec`_      |              | 93M    | .545 | 62.1 | 29.0    |
+----------------+--------------+--------+------+------+---------+
| `TAPE`_        |              | 38M    | .581 | 58.0 | 23.2    |
+----------------+--------------+--------+------+------+---------+
| LSTM biLM (S)  | UR50/S       | 28M    | .558 | 60.4 | 24.1    |
+----------------+--------------+--------+------+------+---------+
| LSTM biLM (L)  | UR50/S       | 113M   | .574 | 62.4 | 27.8    |
+----------------+--------------+--------+------+------+---------+
| Transformer-6  | UR50/S       | 43M    | .653 | 62.0 | 30.2    |
+----------------+--------------+--------+------+------+---------+
| Transformer-12 | UR50/S       | 85M    | .639 | 65.4 | 37.7    |
+----------------+--------------+--------+------+------+---------+
| Transformer-34 | UR100        | 670M   | .599 | 64.3 | 32.7    |
+----------------+--------------+--------+------+------+---------+
| Transformer-34 | UR50/S       | 670M   | .639 | 69.2 | 50.2    |
+----------------+--------------+--------+------+------+---------+
.. _Wang et al. (2017): https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005324
.. _UniRep: https://www.nature.com/articles/s41592-019-0598-1
.. _SeqVec: https://github.com/rostlab/SeqVec

Performance on TAPE benchmark
================

We evaluated our best performing model on the `TAPE`_ benchmark (Rao, et al. 2019), finding that our neural embeddings perform similarly to or better than alignment-based methods.

.. _TAPE: https://github.com/songlab-cal/tape

+--------------------+------+------+-----------------+--------------+-----------+-------------+
| Model              | SS3  | SS8  | Remote homology | Fluorescence | Stability | Contact     |
+--------------------+------+------+-----------------+--------------+-----------+-------------+
| ESM (best neural)  | 0.82 | 0.67 | 0.33            | 0.68         | 0.71      | (0.61)\*    |
+--------------------+------+------+-----------------+--------------+-----------+-------------+
| TAPE (best neural) | 0.75 | 0.59 | 0.26            | 0.68         | 0.73      | 0.4         |
+--------------------+------+------+-----------------+--------------+-----------+-------------+
| TAPE (alignment)   | 0.8  | 0.63 | 0.09            | N/A          | N/A       | 0.64        |
+--------------------+------+------+-----------------+--------------+-----------+-------------+
\* Not comparable: ESM (bests neural) uses a linear projection on the features (the contact head available in the PyTorch version of TAPE),
but the results from the TAPE paper use a ResNet head.
See the previous table for a rigorous comparison of ESM and TAPE in a fair benchmarking setup.

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

Additionally, much of this code hails from the excellent `fairseq`_ sequence modeling framework; we have released this standalone model to facilitate more lightweight and flexible usage. We encourage those who wish to pretrain protein language models from scratch to use fairseq.

.. _fairseq: https://github.com/pytorch/fairseq

License
=======

This source code is licensed under the MIT license found in the ``LICENSE`` file
in the root directory of this source tree.
