import torch
from fairscale.nn.data_parallel import FullyShardedDataParallel as FSDP
from fairscale.nn.wrap import enable_wrap, wrap

import esm

# init the distributed world with world_size 1
url = "tcp://localhost:23456"
torch.distributed.init_process_group(backend="nccl", init_method=url, world_size=1, rank=0)

# download model data from the hub
model_data, regression_data = esm.pretrained._download_model_and_regression_data(
    "esm2_t48_15B_UR50D"
)
if regression_data is not None:
    model_data["model"].update(regression_data["model"])

# initialize the model with FSDP wrapper
fsdp_params = dict(
    mixed_precision=True,
    flatten_parameters=True,
    state_dict_device=torch.device("cpu"),  # reduce GPU mem usage
    cpu_offload=True,  # enable cpu offloading
)
with enable_wrap(wrapper_cls=FSDP, **fsdp_params):
    model, vocab, _ = esm.pretrained._load_model_and_alphabet_core_v2(model_data)
    batch_converter = vocab.get_batch_converter()
    model.eval()

    # Wrap each layer in FSDP separately
    for name, child in model.named_children():
        if name == "layers":
            for layer_name, layer in child.named_children():
                wrapped_layer = wrap(layer)
                setattr(child, layer_name, wrapped_layer)
    model = wrap(model)

data = [
    ("protein1", "MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG"),
    ("protein2", "KALTARQQEVFDLIRDHISQTGMPPTRAEIAQRLGFRSPNAAEEHLKALARKGVIEIVSGASRGIRLLQEE"),
    (
        "protein2 with mask",
        "KALTARQQEVFDLIRD<mask>ISQTGMPPTRAEIAQRLGFRSPNAAEEHLKALARKGVIEIVSGASRGIRLLQEE",
    ),
    ("protein3", "K A <mask> I S Q"),
]

batch_labels, batch_strs, batch_tokens = batch_converter(data)
batch_tokens = batch_tokens.cuda()
with torch.no_grad():
    results = model(batch_tokens, repr_layers=[33], return_contacts=True)
print(results)
