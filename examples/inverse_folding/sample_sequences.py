# Copyright (c) Facebook, Inc. and its affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.
#
# Sample sequences based on a given structure (greedy sampling, no beam search).
#
# usage: sample_sequences.py [-h] [--chain CHAIN] [--temperature TEMPERATURE]
# [--outpath OUTPATH] [--num-samples NUM_SAMPLES] pdbfile

import argparse
import numpy as np
from pathlib import Path

import esm
import esm.inverse_folding


def main():
    parser = argparse.ArgumentParser(
            description='Sample sequences based on a given structure.'
    )
    parser.add_argument(
            'pdbfile', type=str,
            help='input filepath, either .pdb or .cif',
    )
    parser.add_argument(
            '--chain', type=str,
            help='chain id for the chain of interest', default=None,
    )
    parser.add_argument(
            '--temperature', type=float,
            help='temperature for sampling, higher for more diversity',
            default=1.,
    )
    parser.add_argument(
            '--outpath', type=str,
            help='output filepath for saving sampled sequences',
            default='output/sampled_seqs.fasta',
    )
    parser.add_argument(
            '--num-samples', type=int,
            help='number of sequences to sample',
            default=1,
    )
    args = parser.parse_args()

    model, alphabet = esm.pretrained.esm_if1_gvp4_t16_142M_UR50()
    model = model.eval()
    coords, seq = esm.inverse_folding.util.load_coords(args.pdbfile, args.chain)
    print('Sequence loaded from file:')
    print(seq)

    print(f'Saving sampled sequences to {args.outpath}.')

    Path(args.outpath).parent.mkdir(parents=True, exist_ok=True)
    with open(args.outpath, 'w') as f:
        for i in range(args.num_samples):
            print(f'\nSampling.. ({i+1} of {args.num_samples})')
            sampled_seq = model.sample(coords, temperature=args.temperature)
            print('Sampled sequence:')
            print(sampled_seq)
            f.write(f'>sampled_seq_{i+1}\n')
            f.write(sampled_seq + '\n')

            recovery = np.mean([(a==b) for a, b in zip(seq, sampled_seq)])
            print('Sequence recovery:', recovery)


if __name__ == '__main__':
    main()
