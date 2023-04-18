# Copyright (c) Meta Platforms, Inc. and affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

import re
import urllib
import warnings
from argparse import Namespace
from pathlib import Path

import torch

import esm
from esm.model.esm2 import ESM2


def _has_regression_weights(model_name):
    """Return whether we expect / require regression weights;
    Right now that is all models except ESM-1v, ESM-IF, and partially trained ESM2 models"""
    return not ("esm1v" in model_name or "esm_if" in model_name or "270K" in model_name or "500K" in model_name)


def load_model_and_alphabet(model_name):
    if model_name.endswith(".pt"):  # treat as filepath
        return load_model_and_alphabet_local(model_name)
    else:
        return load_model_and_alphabet_hub(model_name)


def load_hub_workaround(url):
    try:
        data = torch.hub.load_state_dict_from_url(url, progress=False, map_location="cpu")
    except RuntimeError:
        # Pytorch version issue - see https://github.com/pytorch/pytorch/issues/43106
        fn = Path(url).name
        data = torch.load(
            f"{torch.hub.get_dir()}/checkpoints/{fn}",
            map_location="cpu",
        )
    except urllib.error.HTTPError as e:
        raise Exception(f"Could not load {url}, check if you specified a correct model name?")
    return data


def load_regression_hub(model_name):
    url = f"https://dl.fbaipublicfiles.com/fair-esm/regression/{model_name}-contact-regression.pt"
    regression_data = load_hub_workaround(url)
    return regression_data


def _download_model_and_regression_data(model_name):
    url = f"https://dl.fbaipublicfiles.com/fair-esm/models/{model_name}.pt"
    model_data = load_hub_workaround(url)
    if _has_regression_weights(model_name):
        regression_data = load_regression_hub(model_name)
    else:
        regression_data = None
    return model_data, regression_data


def load_model_and_alphabet_hub(model_name):
    model_data, regression_data = _download_model_and_regression_data(model_name)
    return load_model_and_alphabet_core(model_name, model_data, regression_data)


def load_model_and_alphabet_local(model_location):
    """Load from local path. The regression weights need to be co-located"""
    model_location = Path(model_location)
    model_data = torch.load(str(model_location), map_location="cpu")
    model_name = model_location.stem
    if _has_regression_weights(model_name):
        regression_location = str(model_location.with_suffix("")) + "-contact-regression.pt"
        regression_data = torch.load(regression_location, map_location="cpu")
    else:
        regression_data = None
    return load_model_and_alphabet_core(model_name, model_data, regression_data)


def has_emb_layer_norm_before(model_state):
    """Determine whether layer norm needs to be applied before the encoder"""
    return any(k.startswith("emb_layer_norm_before") for k, param in model_state.items())


def _load_model_and_alphabet_core_v1(model_data):
    import esm  # since esm.inverse_folding is imported below, you actually have to re-import esm here

    alphabet = esm.Alphabet.from_architecture(model_data["args"].arch)

    if model_data["args"].arch == "roberta_large":
        # upgrade state dict
        pra = lambda s: "".join(s.split("encoder_")[1:] if "encoder" in s else s)
        prs1 = lambda s: "".join(s.split("encoder.")[1:] if "encoder" in s else s)
        prs2 = lambda s: "".join(
            s.split("sentence_encoder.")[1:] if "sentence_encoder" in s else s
        )
        model_args = {pra(arg[0]): arg[1] for arg in vars(model_data["args"]).items()}
        model_state = {prs1(prs2(arg[0])): arg[1] for arg in model_data["model"].items()}
        model_state["embed_tokens.weight"][alphabet.mask_idx].zero_()  # For token drop
        model_args["emb_layer_norm_before"] = has_emb_layer_norm_before(model_state)
        model_type = esm.ProteinBertModel

    elif model_data["args"].arch == "protein_bert_base":

        # upgrade state dict
        pra = lambda s: "".join(s.split("decoder_")[1:] if "decoder" in s else s)
        prs = lambda s: "".join(s.split("decoder.")[1:] if "decoder" in s else s)
        model_args = {pra(arg[0]): arg[1] for arg in vars(model_data["args"]).items()}
        model_state = {prs(arg[0]): arg[1] for arg in model_data["model"].items()}
        model_type = esm.ProteinBertModel
    elif model_data["args"].arch == "msa_transformer":

        # upgrade state dict
        pra = lambda s: "".join(s.split("encoder_")[1:] if "encoder" in s else s)
        prs1 = lambda s: "".join(s.split("encoder.")[1:] if "encoder" in s else s)
        prs2 = lambda s: "".join(
            s.split("sentence_encoder.")[1:] if "sentence_encoder" in s else s
        )
        prs3 = lambda s: s.replace("row", "column") if "row" in s else s.replace("column", "row")
        model_args = {pra(arg[0]): arg[1] for arg in vars(model_data["args"]).items()}
        model_state = {prs1(prs2(prs3(arg[0]))): arg[1] for arg in model_data["model"].items()}
        if model_args.get("embed_positions_msa", False):
            emb_dim = model_state["msa_position_embedding"].size(-1)
            model_args["embed_positions_msa_dim"] = emb_dim  # initial release, bug: emb_dim==1

        model_type = esm.MSATransformer

    elif "invariant_gvp" in model_data["args"].arch:
        import esm.inverse_folding

        model_type = esm.inverse_folding.gvp_transformer.GVPTransformerModel
        model_args = vars(model_data["args"])  # convert Namespace -> dict

        def update_name(s):
            # Map the module names in checkpoints trained with internal code to
            # the updated module names in open source code
            s = s.replace("W_v", "embed_graph.embed_node")
            s = s.replace("W_e", "embed_graph.embed_edge")
            s = s.replace("embed_scores.0", "embed_confidence")
            s = s.replace("embed_score.", "embed_graph.embed_confidence.")
            s = s.replace("seq_logits_projection.", "")
            s = s.replace("embed_ingraham_features", "embed_dihedrals")
            s = s.replace("embed_gvp_in_local_frame.0", "embed_gvp_output")
            s = s.replace("embed_features_in_local_frame.0", "embed_gvp_input_features")
            return s

        model_state = {
            update_name(sname): svalue
            for sname, svalue in model_data["model"].items()
            if "version" not in sname
        }

    else:
        raise ValueError("Unknown architecture selected")

    model = model_type(
        Namespace(**model_args),
        alphabet,
    )

    return model, alphabet, model_state


def _load_model_and_alphabet_core_v2(model_data):
    def upgrade_state_dict(state_dict):
        """Removes prefixes 'model.encoder.sentence_encoder.' and 'model.encoder.'."""
        prefixes = ["encoder.sentence_encoder.", "encoder."]
        pattern = re.compile("^" + "|".join(prefixes))
        state_dict = {pattern.sub("", name): param for name, param in state_dict.items()}
        return state_dict

    cfg = model_data["cfg"]["model"]
    state_dict = model_data["model"]
    state_dict = upgrade_state_dict(state_dict)
    alphabet = esm.data.Alphabet.from_architecture("ESM-1b")
    model = ESM2(
        num_layers=cfg.encoder_layers,
        embed_dim=cfg.encoder_embed_dim,
        attention_heads=cfg.encoder_attention_heads,
        alphabet=alphabet,
        token_dropout=cfg.token_dropout,
    )
    return model, alphabet, state_dict


def load_model_and_alphabet_core(model_name, model_data, regression_data=None):
    if regression_data is not None:
        model_data["model"].update(regression_data["model"])

    if model_name.startswith("esm2"):
        model, alphabet, model_state = _load_model_and_alphabet_core_v2(model_data)
    else:
        model, alphabet, model_state = _load_model_and_alphabet_core_v1(model_data)

    expected_keys = set(model.state_dict().keys())
    found_keys = set(model_state.keys())

    if regression_data is None:
        expected_missing = {"contact_head.regression.weight", "contact_head.regression.bias"}
        error_msgs = []
        missing = (expected_keys - found_keys) - expected_missing
        if missing:
            error_msgs.append(f"Missing key(s) in state_dict: {missing}.")
        unexpected = found_keys - expected_keys
        if unexpected:
            error_msgs.append(f"Unexpected key(s) in state_dict: {unexpected}.")

        if error_msgs:
            raise RuntimeError(
                "Error(s) in loading state_dict for {}:\n\t{}".format(
                    model.__class__.__name__, "\n\t".join(error_msgs)
                )
            )
        if expected_missing - found_keys:
            warnings.warn(
                "Regression weights not found, predicting contacts will not produce correct results."
            )

    model.load_state_dict(model_state, strict=regression_data is not None)

    return model, alphabet


def esm1_t34_670M_UR50S():
    """34 layer transformer model with 670M params, trained on Uniref50 Sparse.

    Returns a tuple of (Model, Alphabet).
    """
    return load_model_and_alphabet_hub("esm1_t34_670M_UR50S")


def esm1_t34_670M_UR50D():
    """34 layer transformer model with 670M params, trained on Uniref50 Dense.

    Returns a tuple of (Model, Alphabet).
    """
    return load_model_and_alphabet_hub("esm1_t34_670M_UR50D")


def esm1_t34_670M_UR100():
    """34 layer transformer model with 670M params, trained on Uniref100.

    Returns a tuple of (Model, Alphabet).
    """
    return load_model_and_alphabet_hub("esm1_t34_670M_UR100")


def esm1_t12_85M_UR50S():
    """12 layer transformer model with 85M params, trained on Uniref50 Sparse.

    Returns a tuple of (Model, Alphabet).
    """
    return load_model_and_alphabet_hub("esm1_t12_85M_UR50S")


def esm1_t6_43M_UR50S():
    """6 layer transformer model with 43M params, trained on Uniref50 Sparse.

    Returns a tuple of (Model, Alphabet).
    """
    return load_model_and_alphabet_hub("esm1_t6_43M_UR50S")


def esm1b_t33_650M_UR50S():
    """33 layer transformer model with 650M params, trained on Uniref50 Sparse.
    This is our best performing model, which will be described in a future publication.

    Returns a tuple of (Model, Alphabet).
    """
    return load_model_and_alphabet_hub("esm1b_t33_650M_UR50S")


def esm_msa1_t12_100M_UR50S():
    warnings.warn(
        "This model had a minor bug in the positional embeddings, "
        "please use ESM-MSA-1b: esm.pretrained.esm_msa1b_t12_100M_UR50S()",
    )
    return load_model_and_alphabet_hub("esm_msa1_t12_100M_UR50S")


def esm_msa1b_t12_100M_UR50S():
    return load_model_and_alphabet_hub("esm_msa1b_t12_100M_UR50S")


def esm1v_t33_650M_UR90S():
    """33 layer transformer model with 650M params, trained on Uniref90.
    This is model 1 of a 5 model ensemble.

    Returns a tuple of (Model, Alphabet).
    """
    return load_model_and_alphabet_hub("esm1v_t33_650M_UR90S_1")


def esm1v_t33_650M_UR90S_1():
    """33 layer transformer model with 650M params, trained on Uniref90.
    This is model 1 of a 5 model ensemble.

    Returns a tuple of (Model, Alphabet).
    """
    return load_model_and_alphabet_hub("esm1v_t33_650M_UR90S_1")


def esm1v_t33_650M_UR90S_2():
    """33 layer transformer model with 650M params, trained on Uniref90.
    This is model 2 of a 5 model ensemble.

    Returns a tuple of (Model, Alphabet).
    """
    return load_model_and_alphabet_hub("esm1v_t33_650M_UR90S_2")


def esm1v_t33_650M_UR90S_3():
    """33 layer transformer model with 650M params, trained on Uniref90.
    This is model 3 of a 5 model ensemble.

    Returns a tuple of (Model, Alphabet).
    """
    return load_model_and_alphabet_hub("esm1v_t33_650M_UR90S_3")


def esm1v_t33_650M_UR90S_4():
    """33 layer transformer model with 650M params, trained on Uniref90.
    This is model 4 of a 5 model ensemble.

    Returns a tuple of (Model, Alphabet).
    """
    return load_model_and_alphabet_hub("esm1v_t33_650M_UR90S_4")


def esm1v_t33_650M_UR90S_5():
    """33 layer transformer model with 650M params, trained on Uniref90.
    This is model 5 of a 5 model ensemble.

    Returns a tuple of (Model, Alphabet).
    """
    return load_model_and_alphabet_hub("esm1v_t33_650M_UR90S_5")


def esm_if1_gvp4_t16_142M_UR50():
    """Inverse folding model with 142M params, with 4 GVP-GNN layers, 8
    Transformer encoder layers, and 8 Transformer decoder layers, trained on
    CATH structures and 12 million alphafold2 predicted structures from UniRef50
    sequences.

    Returns a tuple of (Model, Alphabet).
    """
    return load_model_and_alphabet_hub("esm_if1_gvp4_t16_142M_UR50")


def esm2_t6_8M_UR50D():
    """6 layer ESM-2 model with 8M params, trained on UniRef50.

    Returns a tuple of (Model, Alphabet).
    """
    return load_model_and_alphabet_hub("esm2_t6_8M_UR50D")


def esm2_t12_35M_UR50D():
    """12 layer ESM-2 model with 35M params, trained on UniRef50.

    Returns a tuple of (Model, Alphabet).
    """
    return load_model_and_alphabet_hub("esm2_t12_35M_UR50D")


def esm2_t30_150M_UR50D():
    """30 layer ESM-2 model with 150M params, trained on UniRef50.

    Returns a tuple of (Model, Alphabet).
    """
    return load_model_and_alphabet_hub("esm2_t30_150M_UR50D")


def esm2_t33_650M_UR50D():
    """33 layer ESM-2 model with 650M params, trained on UniRef50.

    Returns a tuple of (Model, Alphabet).
    """
    return load_model_and_alphabet_hub("esm2_t33_650M_UR50D")


def esm2_t36_3B_UR50D():
    """36 layer ESM-2 model with 3B params, trained on UniRef50.

    Returns a tuple of (Model, Alphabet).
    """
    return load_model_and_alphabet_hub("esm2_t36_3B_UR50D")


def esm2_t48_15B_UR50D():
    """48 layer ESM-2 model with 15B params, trained on UniRef50.
    If you have OOM while loading this model, please refer to README
    on how to employ FSDP and ZeRO CPU offloading

    Returns a tuple of (Model, Alphabet).
    """
    return load_model_and_alphabet_hub("esm2_t48_15B_UR50D")


def esmfold_v0():
    """
    ESMFold v0 model with 3B ESM-2, 48 folding blocks.
    This version was used for the paper (Lin et al, 2022). It was trained 
    on all PDB chains until 2020-05, to ensure temporal holdout with CASP14
    and the CAMEO validation and test set reported there.
    """
    import esm.esmfold.v1.pretrained
    return esm.esmfold.v1.pretrained.esmfold_v0()


def esmfold_v1():
    """
    ESMFold v1 model using 3B ESM-2, 48 folding blocks.
    ESMFold provides fast high accuracy atomic level structure prediction
    directly from the individual sequence of a protein. ESMFold uses the ESM2
    protein language model to extract meaningful representations from the
    protein sequence.
    """
    import esm.esmfold.v1.pretrained
    return esm.esmfold.v1.pretrained.esmfold_v1()

def esmfold_structure_module_only_8M():
    """
    ESMFold baseline model using 8M ESM-2, 0 folding blocks.
    ESM-2 here is trained out to 500K updates.
    This is a model designed to test the capabilities of the language model
    when ablated for number of parameters in the language model.
    See table S1 in (Lin et al, 2022).
    """
    import esm.esmfold.v1.pretrained
    return esm.esmfold.v1.pretrained.esmfold_structure_module_only_8M()


def esmfold_structure_module_only_8M_270K():
    """
    ESMFold baseline model using 8M ESM-2, 0 folding blocks.
    ESM-2 here is trained out to 270K updates.
    This is a model designed to test the capabilities of the language model
    when ablated for number of parameters in the language model.
    See table S1 in (Lin et al, 2022).
    """
    import esm.esmfold.v1.pretrained
    return esm.esmfold.v1.pretrained.esmfold_structure_module_only_8M_270K()


def esmfold_structure_module_only_35M():
    """
    ESMFold baseline model using 35M ESM-2, 0 folding blocks.
    ESM-2 here is trained out to 500K updates.
    This is a model designed to test the capabilities of the language model
    when ablated for number of parameters in the language model.
    See table S1 in (Lin et al, 2022).
    """
    import esm.esmfold.v1.pretrained
    return esm.esmfold.v1.pretrained.esmfold_structure_module_only_35M()


def esmfold_structure_module_only_35M_270K():
    """
    ESMFold baseline model using 35M ESM-2, 0 folding blocks.
    ESM-2 here is trained out to 270K updates.
    This is a model designed to test the capabilities of the language model
    when ablated for number of parameters in the language model.
    See table S1 in (Lin et al, 2022).
    """
    import esm.esmfold.v1.pretrained
    return esm.esmfold.v1.pretrained.esmfold_structure_module_only_35M_270K()


def esmfold_structure_module_only_150M():
    """
    ESMFold baseline model using 150M ESM-2, 0 folding blocks.
    ESM-2 here is trained out to 500K updates.
    This is a model designed to test the capabilities of the language model
    when ablated for number of parameters in the language model.
    See table S1 in (Lin et al, 2022).
    """
    import esm.esmfold.v1.pretrained
    return esm.esmfold.v1.pretrained.esmfold_structure_module_only_150M()


def esmfold_structure_module_only_150M_270K():
    """
    ESMFold baseline model using 150M ESM-2, 0 folding blocks.
    ESM-2 here is trained out to 270K updates.
    This is a model designed to test the capabilities of the language model
    when ablated for number of parameters in the language model.
    See table S1 in (Lin et al, 2022).
    """
    import esm.esmfold.v1.pretrained
    return esm.esmfold.v1.pretrained.esmfold_structure_module_only_150M_270K()


def esmfold_structure_module_only_650M():
    """
    ESMFold baseline model using 650M ESM-2, 0 folding blocks.
    ESM-2 here is trained out to 500K updates.
    This is a model designed to test the capabilities of the language model
    when ablated for number of parameters in the language model.
    See table S1 in (Lin et al, 2022).
    """
    import esm.esmfold.v1.pretrained
    return esm.esmfold.v1.pretrained.esmfold_structure_module_only_650M()


def esmfold_structure_module_only_650M_270K():
    """
    ESMFold baseline model using 650M ESM-2, 0 folding blocks.
    ESM-2 here is trained out to 270K updates.
    This is a model designed to test the capabilities of the language model
    when ablated for number of parameters in the language model.
    See table S1 in (Lin et al, 2022).
    """
    import esm.esmfold.v1.pretrained
    return esm.esmfold.v1.pretrained.esmfold_structure_module_only_650M_270K()


def esmfold_structure_module_only_3B():
    """
    ESMFold baseline model using 3B ESM-2, 0 folding blocks.
    ESM-2 here is trained out to 500K updates.
    This is a model designed to test the capabilities of the language model
    when ablated for number of parameters in the language model.
    See table S1 in (Lin et al, 2022).
    """
    import esm.esmfold.v1.pretrained
    return esm.esmfold.v1.pretrained.esmfold_structure_module_only_3B()


def esmfold_structure_module_only_3B_270K():
    """
    ESMFold baseline model using 3B ESM-2, 0 folding blocks.
    ESM-2 here is trained out to 270K updates.
    This is a model designed to test the capabilities of the language model
    when ablated for number of parameters in the language model.
    See table S1 in (Lin et al, 2022).
    """
    import esm.esmfold.v1.pretrained
    return esm.esmfold.v1.pretrained.esmfold_structure_module_only_3B_270K()


def esmfold_structure_module_only_15B():
    """
    ESMFold baseline model using 15B ESM-2, 0 folding blocks.
    ESM-2 here is trained out to 270K updates.
    The 15B parameter ESM-2 was not trained out to 500K updates
    This is a model designed to test the capabilities of the language model
    when ablated for number of parameters in the language model.
    See table S1 in (Lin et al, 2022).
    """
    import esm.esmfold.v1.pretrained
    return esm.esmfold.v1.pretrained.esmfold_structure_module_only_15B()
