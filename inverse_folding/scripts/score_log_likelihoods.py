# Copyright (c) Facebook, Inc. and its affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.
#
# usage: redesign_from_pdbfile.py [-h] [--chain CHAIN] [--temperature TEMPERATURE] PDBFILE

import argparse
from biotite.sequence.io.fasta import FastaFile, get_sequences
import numpy as np
from pathlib import Path
import torch
import torch.nn.functional as F
from tqdm import tqdm

from inverse_folding import (
    load_model_and_alphabet,
    load_coords,
    score_sequence,
)


def main():
    parser = argparse.ArgumentParser(
            description='Score sequences based on a given structure.'
    )
    parser.add_argument(
            'pdbfile', type=str,
            help='input filepath, either .pdb or .cif',
    )
    parser.add_argument(
            'seqfile', type=str,
            help='input filepath for variant sequences in a .fasta file',
    )
    parser.add_argument(
            '--outpath', type=str,
            help='output filepath for scores of variant sequences',
            default='output/sequence_scores.csv',
    )
    parser.add_argument(
            '--chain', type=str,
            help='chain id for the chain of interest', default='A',
    )
    args = parser.parse_args()

    model, alphabet = load_model_and_alphabet()
    coords, seq = load_coords(args.pdbfile, args.chain)
    print('Native sequence loaded from structure file:')
    print(seq)
    print('\n')

    ll, _ = score_sequence(model, alphabet, coords, seq) 
    print('Native sequence')
    print(f'Log likelihood: {ll:.2f}')
    print(f'Perplexity: {np.exp(-ll):.2f}')

    print('\nScoring variant sequences from sequence file..\n')
    infile = FastaFile()
    infile.read(args.seqfile)
    seqs = get_sequences(infile)
    Path(args.outpath).parent.mkdir(parents=True, exist_ok=True)
    with open(args.outpath, 'w') as fout:
        fout.write('seqid,log_likelihood\n')
        for header, seq in tqdm(seqs.items()):
            ll, _ = score_sequence(model, alphabet, coords, str(seq))
            fout.write(header + ',' + str(ll) + '\n')
    print(f'Results saved to {args.outpath}') 


if __name__ == '__main__':
    main()
