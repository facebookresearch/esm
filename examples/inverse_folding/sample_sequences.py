# Copyright (c) Meta Platforms, Inc. and affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.
#
# Sample sequences based on a given structure (multinomial sampling, no beam search).
#
# usage: sample_sequences.py [-h] [--chain CHAIN] [--temperature TEMPERATURE]
# [--outpath OUTPATH] [--num-samples NUM_SAMPLES] pdbfile

import argparse
import numpy as np
from pathlib import Path
import torch

import esm
import esm.inverse_folding


def sample_seq_singlechain(model, alphabet, args):
    if torch.cuda.is_available() and not args.nogpu:
        model = model.cuda()
        print("Transferred model to GPU")
    coords, native_seq = esm.inverse_folding.util.load_coords(args.pdbfile, args.chain)
    print('Native sequence loaded from structure file:')
    print(native_seq)

    print(f'Saving sampled sequences to {args.outpath}.')

    Path(args.outpath).parent.mkdir(parents=True, exist_ok=True)
    with open(args.outpath, 'w') as f:
        for i in range(args.num_samples):
            print(f'\nSampling.. ({i+1} of {args.num_samples})')
            sampled_seq = model.sample(coords, temperature=args.temperature, device=torch.device('cuda'))
            print('Sampled sequence:')
            print(sampled_seq)
            f.write(f'>sampled_seq_{i+1}\n')
            f.write(sampled_seq + '\n')

            recovery = np.mean([(a==b) for a, b in zip(native_seq, sampled_seq)])
            print('Sequence recovery:', recovery)


def sample_seq_multichain(model, alphabet, args):
    if torch.cuda.is_available() and not args.nogpu:
        model = model.cuda()
        print("Transferred model to GPU")
    structure = esm.inverse_folding.util.load_structure(args.pdbfile)
    coords, native_seqs = esm.inverse_folding.multichain_util.extract_coords_from_complex(structure)
    target_chain_id = args.chain
    native_seq = native_seqs[target_chain_id]
    print('Native sequence loaded from structure file:')
    print(native_seq)
    print('\n')

    print(f'Saving sampled sequences to {args.outpath}.')

    Path(args.outpath).parent.mkdir(parents=True, exist_ok=True)
    with open(args.outpath, 'w') as f:
        for i in range(args.num_samples):
            print(f'\nSampling.. ({i+1} of {args.num_samples})')
            sampled_seq = esm.inverse_folding.multichain_util.sample_sequence_in_complex(
                    model, coords, target_chain_id, temperature=args.temperature)
            print('Sampled sequence:')
            print(sampled_seq)
            f.write(f'>sampled_seq_{i+1}\n')
            f.write(sampled_seq + '\n')

            recovery = np.mean([(a==b) for a, b in zip(native_seq, sampled_seq)])
            print('Sequence recovery:', recovery)


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
    parser.set_defaults(multichain_backbone=False)
    parser.add_argument(
            '--multichain-backbone', action='store_true',
            help='use the backbones of all chains in the input for conditioning'
    )
    parser.add_argument(
            '--singlechain-backbone', dest='multichain_backbone',
            action='store_false',
            help='use the backbone of only target chain in the input for conditioning'
    )
    parser.add_argument("--nogpu", action="store_true", help="Do not use GPU even if available")
  
    args = parser.parse_args()

    model, alphabet = esm.pretrained.esm_if1_gvp4_t16_142M_UR50()
    model = model.eval()

    if args.multichain_backbone:
        sample_seq_multichain(model, alphabet, args)
    else:
        sample_seq_singlechain(model, alphabet, args)


if __name__ == '__main__':
    main()
