# Copyright (c) Meta Platforms, Inc. and affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

import biotite.structure
import numpy as np
import torch
from typing import Sequence, Tuple, List

from esm.inverse_folding.util import (
    load_structure,
    extract_coords_from_structure,
    load_coords,
    get_sequence_loss,
    get_encoder_output,
)


def extract_coords_from_complex(structure: biotite.structure.AtomArray):
    """
    Args:
        structure: biotite AtomArray
    Returns:
        Tuple (coords_list, seq_list)
        - coords: Dictionary mapping chain ids to L x 3 x 3 array for N, CA, C
          coordinates representing the backbone of each chain
        - seqs: Dictionary mapping chain ids to native sequences of each chain
    """
    coords = {}
    seqs = {}
    all_chains = biotite.structure.get_chains(structure)
    for chain_id in all_chains:
        chain = structure[structure.chain_id == chain_id]
        coords[chain_id], seqs[chain_id] = extract_coords_from_structure(chain)
    return coords, seqs


def load_complex_coords(fpath, chains):
    """
    Args:
        fpath: filepath to either pdb or cif file
        chains: the chain ids (the order matters for autoregressive model)
    Returns:
        Tuple (coords_list, seq_list)
        - coords: Dictionary mapping chain ids to L x 3 x 3 array for N, CA, C
          coordinates representing the backbone of each chain
        - seqs: Dictionary mapping chain ids to native sequences of each chain
    """
    structure = load_structure(fpath, chains)
    return extract_coords_from_complex(structure)


def _concatenate_coords(coords, target_chain_id, padding_length=10):
    """
    Args:
        coords: Dictionary mapping chain ids to L x 3 x 3 array for N, CA, C
            coordinates representing the backbone of each chain
        target_chain_id: The chain id to sample sequences for
        padding_length: Length of padding between concatenated chains
    Returns:
        Tuple (coords, seq)
            - coords is an L x 3 x 3 array for N, CA, C coordinates, a
              concatenation of the chains with padding in between
            - seq is the extracted sequence, with padding tokens inserted
              between the concatenated chains
    """
    pad_coords = np.full((padding_length, 3, 3), np.nan, dtype=np.float32)
    # For best performance, put the target chain first in concatenation.
    coords_list = [coords[target_chain_id]]
    for chain_id in coords:
        if chain_id == target_chain_id:
            continue
        coords_list.append(pad_coords)
        coords_list.append(coords[chain_id])
    coords_concatenated = np.concatenate(coords_list, axis=0)
    return coords_concatenated


def sample_sequence_in_complex(model, coords, target_chain_id, temperature=1.,
        padding_length=10):
    """
    Samples sequence for one chain in a complex.
    Args:
        model: An instance of the GVPTransformer model
        coords: Dictionary mapping chain ids to L x 3 x 3 array for N, CA, C
            coordinates representing the backbone of each chain
        target_chain_id: The chain id to sample sequences for
        padding_length: padding length in between chains
    Returns:
        Sampled sequence for the target chain
    """
    target_chain_len = coords[target_chain_id].shape[0]
    all_coords = _concatenate_coords(coords, target_chain_id)
    device = next(model.parameters()).device

    # Supply padding tokens for other chains to avoid unused sampling for speed
    padding_pattern = ['<pad>'] * all_coords.shape[0]
    for i in range(target_chain_len):
        padding_pattern[i] = '<mask>'
    sampled = model.sample(all_coords, partial_seq=padding_pattern,
            temperature=temperature, device=device)
    sampled = sampled[:target_chain_len]
    return sampled


def score_sequence_in_complex(model, alphabet, coords, target_chain_id,
        target_seq, padding_length=10):
    """
    Scores sequence for one chain in a complex.
    Args:
        model: An instance of the GVPTransformer model
        alphabet: Alphabet for the model
        coords: Dictionary mapping chain ids to L x 3 x 3 array for N, CA, C
            coordinates representing the backbone of each chain
        target_chain_id: The chain id to sample sequences for
        target_seq: Target sequence for the target chain for scoring.
        padding_length: padding length in between chains
    Returns:
        Tuple (ll_fullseq, ll_withcoord)
        - ll_fullseq: Average log-likelihood over the full target chain
        - ll_withcoord: Average log-likelihood in target chain excluding those
            residues without coordinates
    """
    all_coords = _concatenate_coords(coords, target_chain_id)

    loss, target_padding_mask = get_sequence_loss(model, alphabet, all_coords,
            target_seq)
    ll_fullseq = -np.sum(loss * ~target_padding_mask) / np.sum(
            ~target_padding_mask)

    # Also calculate average when excluding masked portions
    coord_mask = np.all(np.isfinite(coords[target_chain_id]), axis=(-1, -2))
    ll_withcoord = -np.sum(loss * coord_mask) / np.sum(coord_mask)
    return ll_fullseq, ll_withcoord


def get_encoder_output_for_complex(model, alphabet, coords, target_chain_id):
    """
    Args:
        model: An instance of the GVPTransformer model
        alphabet: Alphabet for the model
        coords: Dictionary mapping chain ids to L x 3 x 3 array for N, CA, C
            coordinates representing the backbone of each chain
        target_chain_id: The chain id to sample sequences for
    Returns:
        Dictionary mapping chain id to encoder output for each chain
    """
    all_coords = _concatenate_coords(coords, target_chain_id)
    all_rep = get_encoder_output(model, alphabet, all_coords)
    target_chain_len = coords[target_chain_id].shape[0]
    return all_rep[:target_chain_len]
