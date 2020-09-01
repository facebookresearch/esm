# Copyright (c) Facebook, Inc. and its affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

import math

import torch
import torch.nn as nn
import torch.nn.functional as F

from .modules import TransformerLayer, PositionalEmbedding  # noqa


class ProteinBertModel(nn.Module):
    @classmethod
    def add_args(cls, parser):
        parser.add_argument(
            "--num_layers", default=36, type=int, metavar="N", help="number of layers"
        )
        parser.add_argument(
            "--embed_dim", default=1280, type=int, metavar="N", help="embedding dimension"
        )
        parser.add_argument(
            "--logit_bias", action="store_true", help="whether to apply bias to logits"
        )
        parser.add_argument(
            "--ffn_embed_dim",
            default=5120,
            type=int,
            metavar="N",
            help="embedding dimension for FFN",
        )
        parser.add_argument(
            "--attention_heads",
            default=20,
            type=int,
            metavar="N",
            help="number of attention heads",
        )

    def __init__(self, args, alphabet_size, padding_idx):
        super().__init__()
        self.args = args
        self.alphabet_size = alphabet_size
        self.padding_idx = padding_idx
        self.embed_scale = math.sqrt(self.args.embed_dim)
        self._init_submodules()

    def _init_submodules(self):
        self.embed_tokens = nn.Embedding(
            self.alphabet_size, self.args.embed_dim, padding_idx=self.padding_idx
        )
        self.embed_positions = PositionalEmbedding(self.args.embed_dim, self.padding_idx)
        self.layers = nn.ModuleList(
            [
                TransformerLayer(
                    self.args.embed_dim, self.args.ffn_embed_dim, self.args.attention_heads
                )
                for _ in range(self.args.layers)
            ]
        )
        self.embed_out = nn.Parameter(
            torch.zeros((self.alphabet_size, self.args.embed_dim))
        )
        self.embed_out_bias = None
        if self.args.final_bias:
            self.embed_out_bias = nn.Parameter(torch.zeros(self.alphabet_size))

    def forward(self, tokens, repr_layers=[]):
        assert tokens.ndim == 2
        padding_mask = tokens.eq(self.padding_idx)
        if not padding_mask.any():
            padding_mask = None

        x = self.embed_scale * self.embed_tokens(tokens)
        x = x + self.embed_positions(tokens)

        repr_layers = set(repr_layers)
        hidden_representations = {}
        if 0 in repr_layers:
            hidden_representations[0] = x

        # (B, T, E) => (T, B, E)
        x = x.transpose(0, 1)

        for layer_idx, layer in enumerate(self.layers):
            x, _ = layer(x, self_attn_padding_mask=padding_mask)
            if (layer_idx + 1) in repr_layers:
                hidden_representations[layer_idx + 1] = x.transpose(0, 1)

        x = F.linear(x, self.embed_out, bias=self.embed_out_bias)

        # (T, B, E) => (B, T, E)
        x = x.transpose(0, 1)

        result = {"logits": x, "representations": hidden_representations}

        return result

    @property
    def num_layers(self):
        return self.args.layers
