# Copyright (c) Meta Platforms, Inc. and affiliates.

# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.


from pathlib import Path
import sys,os
import argparse
import logging
import sys
import typing as T
from pathlib import Path
from timeit import default_timer as timer

import torch

import esm
from esm.data import read_fasta

logger = logging.getLogger()
logger.setLevel(logging.INFO)

formatter = logging.Formatter(
    "%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    datefmt="%y/%m/%d %H:%M:%S",
)

console_handler = logging.StreamHandler(sys.stdout)
console_handler.setLevel(logging.INFO)
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)


PathLike = T.Union[str, Path]


def enable_cpu_offloading(model):
    from torch.distributed.fsdp import CPUOffload, FullyShardedDataParallel
    from torch.distributed.fsdp.wrap import enable_wrap, wrap

    torch.distributed.init_process_group(
        backend="nccl", init_method="tcp://localhost:9999", world_size=1, rank=0
    )

    wrapper_kwargs = dict(cpu_offload=CPUOffload(offload_params=True))

    with enable_wrap(wrapper_cls=FullyShardedDataParallel, **wrapper_kwargs):
        for layer_name, layer in model.layers.named_children():
            wrapped_layer = wrap(layer)
            setattr(model.layers, layer_name, wrapped_layer)
        model = wrap(model)

    return model


def init_model_on_gpu_with_cpu_offloading(model):
    model = model.eval()
    model_esm = enable_cpu_offloading(model.esm)
    del model.esm
    model.cuda()
    model.esm = model_esm
    return model


def create_batched_sequence_datasest(
    sequences: T.List[T.Tuple[str, str]], max_tokens_per_batch: int = 1024
) -> T.Generator[T.Tuple[T.List[str], T.List[str]], None, None]:

    batch_headers, batch_sequences, num_tokens = [], [], 0
    for header, seq in sequences:
        if (len(seq) + num_tokens > max_tokens_per_batch) and num_tokens > 0:
            yield batch_headers, batch_sequences
            batch_headers, batch_sequences, num_tokens = [], [], 0
        batch_headers.append(header)
        batch_sequences.append(seq)
        num_tokens += len(seq)

    yield batch_headers, batch_sequences


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--fasta",
        help="Path to input FASTA file",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "-o", "--pdb", help="Path to output PDB directory", type=Path, required=True
    )
    parser.add_argument(
        "-m", "--model-dir", help="Parent path to Pretrained ESM data directory. ", type=Path, default=None
    )
    parser.add_argument(
        "--num-recycles",
        type=int,
        default=None,
        help="Number of recycles to run. Defaults to number used in training (4).",
    )
    parser.add_argument(
        "--max-tokens-per-batch",
        type=int,
        default=1024,
        help="Maximum number of tokens per gpu forward-pass. This will group shorter sequences together "
        "for batched prediction. Lowering this can help with out of memory issues, if these occur on "
        "short sequences.",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=None,
        help="Chunks axial attention computation to reduce memory usage from O(L^2) to O(L). "
        "Equivalent to running a for loop over chunks of of each dimension. Lower values will "
        "result in lower memory usage at the cost of speed. Recommended values: 128, 64, 32. "
        "Default: None.",
    )
    parser.add_argument("--cpu-only", help="CPU only", action="store_true")
    parser.add_argument("--cpu-offload", help="Enable CPU offloading", action="store_true")
    return parser


def run(args):
    if not args.fasta.exists():
        raise FileNotFoundError(args.fasta)

    args.pdb.mkdir(exist_ok=True)

    # Read fasta and sort sequences by length
    logger.info(f"Reading sequences from {args.fasta}")
    all_sequences = sorted(read_fasta(args.fasta), key=lambda header_seq: len(header_seq[1]))
    logger.info(f"Loaded {len(all_sequences)} sequences from {args.fasta}")

    logger.info("Loading model")

    # Use pre-downloaded ESM weights from model_pth.
    if args.model_dir is not None:
        # if pretrained model path is available
        torch.hub.set_dir(args.model_dir)

    model = esm.pretrained.esmfold_v1()


    model = model.eval()
    model.set_chunk_size(args.chunk_size)

    if args.cpu_only:
        model.esm.float()  # convert to fp32 as ESM-2 in fp16 is not supported on CPU
        model.cpu()
    elif args.cpu_offload:
        model = init_model_on_gpu_with_cpu_offloading(model)
    else:
        model.cuda()
    logger.info("Starting Predictions")
    batched_sequences = create_batched_sequence_datasest(all_sequences, args.max_tokens_per_batch)

    num_completed = 0
    num_sequences = len(all_sequences)
    for headers, sequences in batched_sequences:
        start = timer()
        try:
            output = model.infer(sequences, num_recycles=args.num_recycles)
        except RuntimeError as e:
            if e.args[0].startswith("CUDA out of memory"):
                if len(sequences) > 1:
                    logger.info(
                        f"Failed (CUDA out of memory) to predict batch of size {len(sequences)}. "
                        "Try lowering `--max-tokens-per-batch`."
                    )
                else:
                    logger.info(
                        f"Failed (CUDA out of memory) on sequence {headers[0]} of length {len(sequences[0])}."
                    )

                continue
            raise

        output = {key: value.cpu() for key, value in output.items()}
        pdbs = model.output_to_pdb(output)
        tottime = timer() - start
        time_string = f"{tottime / len(headers):0.1f}s"
        if len(sequences) > 1:
            time_string = time_string + f" (amortized, batch size {len(sequences)})"
        for header, seq, pdb_string, mean_plddt, ptm in zip(
            headers, sequences, pdbs, output["mean_plddt"], output["ptm"]
        ):
            output_file = args.pdb / f"{header}.pdb"
            output_file.write_text(pdb_string)
            num_completed += 1
            logger.info(
                f"Predicted structure for {header} with length {len(seq)}, pLDDT {mean_plddt:0.1f}, "
                f"pTM {ptm:0.3f} in {time_string}. "
                f"{num_completed} / {num_sequences} completed."
            )


def main():
    parser = create_parser()
    args = parser.parse_args()
    run(args)

if __name__ == "__main__":
    main()
