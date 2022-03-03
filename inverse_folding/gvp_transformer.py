# Copyright (c) Facebook, Inc. and its affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

import argparse
from typing import Any, Dict, List, Optional, Tuple, NamedTuple
import torch
from torch import nn
from torch import Tensor
import torch.nn.functional as F
from scipy.spatial import transform

from esm.data import Alphabet

from .features import DihedralFeatures
from .gvp_encoder import GVPEncoder
from .gvp_utils import unflatten_graph
from .transformer_encoder import TransformerEncoder
from .transformer_decoder import TransformerDecoder
from .util import rotate, CoordBatchConverter 


class GVPTransformerModel(nn.Module):
    def __init__(self, args, encoder, decoder):
        super().__init__()
        self.args = args
        self.encoder = encoder
        self.decoder = decoder

    @classmethod
    def build_encoder(cls, args, src_dict, embed_tokens):
        encoder = TransformerEncoder(args, src_dict, embed_tokens)
        return encoder

    @classmethod
    def build_decoder(cls, args, tgt_dict, embed_tokens):
        decoder = TransformerDecoder(
            args,
            tgt_dict,
            embed_tokens,
        )
        return decoder

    @classmethod
    def build_model(cls, args, dictionary):
        encoder_embed_tokens = cls.build_embedding(
            args, dictionary, args.encoder_embed_dim,
        )
        decoder_embed_tokens = cls.build_embedding(
            args, dictionary, args.decoder_embed_dim, 
        )
        encoder = cls.build_encoder(args, dictionary, encoder_embed_tokens)
        decoder = cls.build_decoder(args, dictionary, decoder_embed_tokens)
        return cls(args, encoder, decoder)

    @classmethod
    def build_embedding(cls, args, dictionary, embed_dim):
        num_embeddings = len(dictionary)
        padding_idx = dictionary.padding_idx
        emb = nn.Embedding(num_embeddings, embed_dim, padding_idx)
        nn.init.normal_(emb.weight, mean=0, std=embed_dim ** -0.5)
        nn.init.constant_(emb.weight[padding_idx], 0)
        return emb

    def forward(
        self,
        coords,
        padding_mask,
        confidence,
        prev_output_tokens,
        return_all_hiddens: bool = False,
        features_only: bool = False,
    ):
        encoder_out = self.encoder(coords, padding_mask, confidence,
            return_all_hiddens=return_all_hiddens)
        logits, extra = self.decoder(
            prev_output_tokens,
            encoder_out=encoder_out,
            features_only=features_only,
            return_all_hiddens=return_all_hiddens,
        )
        return logits, extra
    
    def sample(self, coords, temperature=1.0, confidence=None):
        """
        Args:
            coords: L x 3 x 3 list representing one backbone
            temperature: sampling temperature, use low temperature for higher
                sequence recovery and high temperature for higher diversity
            confidence: optional length L list of confidence scores for coordinates
        """
        L = len(coords)
        # Convert to batch format
        batch_converter = CoordBatchConverter(self.decoder.dictionary)
        batch_coords, confidence, _, _, padding_mask = (
            batch_converter([(coords, confidence, None)])
        )
        
        # Start with prepend token
        sampled_tokens = torch.zeros(1, 1+L, dtype=int)
        sampled_tokens[0, 0] = self.decoder.dictionary.get_idx('<cath>')
            
        # Save incremental states for faster sampling
        incremental_state = dict()
        
        # Run encoder only once
        encoder_out = self.encoder(batch_coords, padding_mask, confidence)
        
        # Decode one token at a time
        for i in range(1, L+1):
            logits, _ = self.decoder(
                sampled_tokens[:, :i], 
                encoder_out,
                incremental_state=incremental_state,
            )
            logits = logits[0].transpose(0, 1)
            logits /= temperature
            probs = F.softmax(logits, dim=-1)
            sampled_tokens[:, i] = torch.multinomial(probs, 1).squeeze(-1)
        sampled_seq = sampled_tokens[0, 1:]
        
        # Convert back to string via lookup
        return ''.join([self.decoder.dictionary.get_tok(a) for a in sampled_seq])


def gvp_transformer_architecture(args):
    # Transformer args
    args.encoder_embed_dim = getattr(args, "encoder_embed_dim", 512)
    args.encoder_ffn_embed_dim = getattr(args, "encoder_ffn_embed_dim", 2048)
    args.encoder_attention_heads = getattr(args, "encoder_attention_heads", 8)
    args.encoder_layers = getattr(args, "encoder_layers", 8)
    args.decoder_embed_dim = getattr(args, "decoder_embed_dim", 512)
    args.decoder_ffn_embed_dim = getattr(args, "decoder_ffn_embed_dim", 2048)
    args.decoder_attention_heads = getattr(args, "decoder_attention_heads", 8)
    args.decoder_layers = getattr(args, "decoder_layers", 8)
    args.attention_dropout = getattr(args, "attention_dropout", 0.1)
    args.dropout = getattr(args, "dropout", 0.1)

    # GVP encoder args
    args.gvp_num_encoder_layers = getattr(args, "gvp_num_encoder_layers", 4)
    args.gvp_dropout = getattr(args, "gvp_dropout", 0.1)
    args.gvp_top_k_neighbors = getattr(args, "gvp_top_k_neighbors", 30)
    args.gvp_num_positional_embeddings = getattr(args,
            "gvp_num_positional_embeddings", 16)
    args.gvp_remove_edges_without_coords = getattr(args,
            "gvp_remove_edges_without_coords", True)
    args.gvp_node_hidden_dim_scalar = getattr(args,
            "gvp_node_hidden_dim_scalar", 1024)
    args.gvp_node_hidden_dim_vector = getattr(args,
            "gvp_node_hidden_dim_vector", 256)
    args.gvp_edge_hidden_dim_scalar = getattr(args,
            "gvp_edge_hidden_dim_scalar", 32)
    args.gvp_edge_hidden_dim_vector = getattr(args,
            "gvp_edge_hidden_dim_vector", 1)


def load_model_and_alphabet(eval_mode=True,
        path='/home/chloehsu/inverse_opensource/model.pt'):
    # TODO: update path to torchhub
    args = argparse.Namespace()
    gvp_transformer_architecture(args)
    alphabet = Alphabet.from_architecture("inverse_folding")
    model = GVPTransformerModel.build_model(args, alphabet)
    # TODO: update this to torchhub path
    path = '/home/chloehsu/inverse_opensource/model.pt'
    with open(path, "rb") as f:
        state = torch.load(f, map_location=torch.device("cpu"))
    model.load_state_dict(state, strict=True)
    model = model.cpu()
    if eval_mode:
        model = model.eval()
    return model, alphabet


def test_model():
    import json
    import numpy as np
    from tqdm import tqdm
    from scipy.stats import special_ortho_group
    from pathlib import Path
    example_file = Path(__file__).absolute().parent / "example.json"
    with open(example_file) as f:
        examples = json.load(f)

    model, alphabet = load_model_and_alphabet()
    batch_converter = CoordBatchConverter(alphabet)

    with torch.no_grad():
        print('Testing batch inference on 3 examples...')
        # Test batch with multiple examples
        batch = [(e["coords"], None, e["seq"]) for e in examples[:3]]
        coords, confidence, strs, tokens, padding_mask = (
            batch_converter(batch)
        )
        prev_output_tokens = tokens[:, :-1]
        target = tokens[:, 1:]
        logits, _ = model.forward(coords, padding_mask, confidence,
                prev_output_tokens)
        loss = torch.nn.functional.cross_entropy(logits, target, reduction='none')
        coord_mask = torch.all(torch.all(torch.isfinite(coords), dim=-1), dim=-1)
        coord_mask = coord_mask[:, 1:-1]
        avgloss = torch.sum(loss * coord_mask) / torch.sum(coord_mask)
        print('ppl:', torch.exp(avgloss).item())

        print('Testing on 10 examples from validation set...')
        # Test batch with single example
        for example in tqdm(examples):
            batch = [(example["coords"], None, example["seq"])]
            coords, confidence, strs, tokens, padding_mask = (
                batch_converter(batch)
            )
            prev_output_tokens = tokens[:, :-1]
            target = tokens[:, 1:]
            logits, _ = model.forward(coords, padding_mask, confidence,
                    prev_output_tokens)
            assert torch.any(torch.isnan(logits)) == False

            # Test equivariance
            R = special_ortho_group.rvs(3)
            R = torch.tensor(R, dtype=torch.float32)
            coords = torch.matmul(coords, R)
            logits_rotated, _ = model.forward(coords, padding_mask,
                    confidence, prev_output_tokens)
            np.testing.assert_allclose(
                    logits.detach().numpy(), 
                    logits_rotated.detach().numpy(), 
                    atol=1e-01
            )


if  __name__ == "__main__":
    test_model()
