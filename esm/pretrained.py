# Copyright (c) Facebook, Inc. and its affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

import esm
import torch
from argparse import Namespace
from .constants import proteinseq_toks

def load_model_and_alphabet(model_name):
    if model_name.endswith(".pt"):  # treat as filepath
        return load_model_and_alphabet_local(model_name)
    else:
        return load_model_and_alphabet_hub(model_name)

def load_model_and_alphabet_hub(model_name):
    alphabet = esm.Alphabet.from_dict(proteinseq_toks)

    url = f"https://dl.fbaipublicfiles.com/fair-esm/models/{model_name}.pt"
    if torch.cuda.is_available():
        model_data = torch.hub.load_state_dict_from_url(url, progress=False)
    else:
        model_data = torch.hub.load_state_dict_from_url(url, progress=False, map_location=torch.device('cpu'))

    # upgrade state dict
    pra = lambda s: ''.join(s.split('decoder_')[1:] if 'decoder' in s else s)
    prs = lambda s: ''.join(s.split('decoder.')[1:] if 'decoder' in s else s)
    model_args = {pra(arg[0]): arg[1] for arg in vars(model_data["args"]).items()}
    model_state = {prs(arg[0]): arg[1] for arg in model_data["model"].items()}

    model = esm.ProteinBertModel(
        Namespace(**model_args), len(alphabet), padding_idx=alphabet.padding_idx
    )
    model.load_state_dict(model_state)

    return model, alphabet

def load_model_and_alphabet_local(model_location):
    alphabet = esm.Alphabet.from_dict(proteinseq_toks)

    model_data = torch.load(model_location)

    # upgrade state dict
    pra = lambda s: ''.join(s.split('decoder_')[1:] if 'decoder' in s else s)
    prs = lambda s: ''.join(s.split('decoder.')[1:] if 'decoder' in s else s)
    model_args = {pra(arg[0]): arg[1] for arg in vars(model_data["args"]).items()}
    model_state = {prs(arg[0]): arg[1] for arg in model_data["model"].items()}

    model = esm.ProteinBertModel(
        Namespace(**model_args), len(alphabet), padding_idx=alphabet.padding_idx
    )
    model.load_state_dict(model_state)
    return model, alphabet

def esm1_t34_670M_UR50S_local():
    model_location = '/checkpoint/bioseq_nonsecure/br2020/br4/checkpoint94.pt'
    model, alphabet = load_model_and_alphabet_local(model_location)

    return model, alphabet

def esm1_t34_670M_UR50S_hub():
    return load_model_and_alphabet_hub("esm1_t34_670M_UR50S")

def esm1_t34_670M_UR50S():
    """ 34 layer transformer model with 670M params, trained on Uniref50 Sparse.

    Returns a tuple of (ProteinBertModel, Alphabet).
    """
    #return esm1_t34_670M_UR50S_hub()
    #return esm1_t34_670M_UR50S_local()
    return load_model_and_alphabet_hub("esm1_t34_670M_UR50S")

def esm1_t34_670M_UR50D():
    """ 34 layer transformer model with 670M params, trained on Uniref50 Dense.

    Returns a tuple of (ProteinBertModel, Alphabet).
    """
    return load_model_and_alphabet_hub("esm1_t34_670M_UR50D")

def esm1_t34_670M_UR100():
    """ 34 layer transformer model with 670M params, trained on Uniref100.

    Returns a tuple of (ProteinBertModel, Alphabet).
    """
    return load_model_and_alphabet_hub("esm1_t34_670M_UR100")

def esm1_t12_85M_UR50S():
    """ 12 layer transformer model with 85M params, trained on Uniref50 Sparse.

    Returns a tuple of (ProteinBertModel, Alphabet).
    """
    return load_model_and_alphabet_hub("esm1_t12_85M_UR50S")

def esm1_t6_43M_UR50S():
    """ 6 layer transformer model with 43M params, trained on Uniref50 Sparse.

    Returns a tuple of (ProteinBertModel, Alphabet).
    """
    return load_model_and_alphabet_hub("esm1_t6_43M_UR50S")
