# Copyright (c) Meta Platforms, Inc. and affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.
#
from collections import Counter
from pathlib import Path
import pickle
import random
import time

from nltk import ngrams
import numpy as np


# Code for loading constants for preprocessing
seq_encode = ['L', 'A', 'G', 'V', 'S', 'E', 'R', 'T', 'I', 'D', 'P', 'K', 'Q', 'N', 'F', 'Y', 'M', 'H', 'W', 'C']
BASE = Path('./utils/ngram_stats/')
ngram_list = []
for fn in ["monogram_seg.p", "bigram_seg.p", "trigram_seg.p", "quadgram_seg.p"]:
    with open(BASE / fn, "rb") as f:
        ngram_list.append(pickle.load(f))

# Recompute ngram frequency based off the valid sequences used for design
for i, ngram_dict in enumerate(ngram_list):
    idx_dict = {}
    for k, v in ngram_dict.items():
        ids = []
        error = False

        for ki in k:
            if ki not in seq_encode:
                error = True
                break
            id = seq_encode.index(ki)
            ids.append(id)

        if error:
            continue

        ids = tuple(ids)
        idx_dict[ids] = v

    total = sum(idx_dict.values())
    # Min value for ngram is 1e-5
    idx_dict = {k: max(v / total, 1e-5) for k, v in idx_dict.items()}

    ngram_list[i] = idx_dict

def encode(seq):
    if isinstance(seq, np.ndarray):
        return seq  # already encoded
    elif isinstance(seq, str):
        return np.array([seq_encode.index(AA) for AA in seq])
    else:
        raise ValueError(f'Unknown seq type {seq}')


def compute_kl_div(seq, order):
   # Inputs Args:
   # Seq: N dimensional numpy array consisting of numbers between 0 and 19 (inclusive)
   # Order: integer for order of ngram used (should be between 0 and 3 for now)
   order_dict = ngram_list[order-1]
   seq = encode(seq) # this is not the problem.
 
   # Compute ngram frequency rate for the input sequence
   tup_dict = Counter(ngrams(seq,n=order))
   total = sum(tup_dict.values())
   tup_dict = {k: v / total for k, v in tup_dict.items()}
  
   p = np.array(list(tup_dict.values())) # observed probabilities of ngrams
   q = np.array([order_dict.get(k, 1e-5) for k in tup_dict.keys()]) # learned ngram probabilities
   return np.sum(p * np.log(p/q))
