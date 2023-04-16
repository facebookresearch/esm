# ESM Metagenomic Atlas

The [ESM Metagenomic Atlas](https://esmatlas.com) is a resource of more than 700 million metagenomic protein structures predicted by ESMFold.
See our [blog post](https://ai.facebook.com/blog/protein-folding-esmfold-metagenomics/) to learn more.

The first `v0` version of the Atlas was released on November 1st 2022,
corresponding to the sequences in the `2022_05` release of the MGnify protein database described [here](https://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/2022_05/README.txt).
An update `v2023_02` was released on March 17th 2023, corresponding to the sequences in the `2023_02` release of the MGnify protein database described [here](https://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/2023_02/).

Here we provide bulk download instructions for the tarballs of predicted structures, as well as foldseek databases, and LM embeddings for every sequence in the Atlas.

The structures in the ESM Metagenomic Atlas were predicted with `esm.pretrained.esmfold_v0()` for Atlas `v0`, while Atlas `v2023_02` used `esm.pretrained.esmfold_v1()`.
We find that protein structures with predicted LDDT > 0.7 and predict TM > 0.7 to be both reasonably well structured and interesting.
Therefore, we provide both the small set of "high confidence" metagenomic structures from Atlas `v0`, as well as the full set.
The small set of structures is built from taking a 30% sequence identity clustering of MGnify90 `2022_05`, and using the best structure from each cluster.
The best structure is selected using the score pTM * pLDDT.

The high confidence structures are around 1TB in size.

The full database is available as PDB structures and is 15TB in size for `v0`.

As an entrypoint into the database, the metadata is available as a small (~25 GB) file [metadata.sqlite](https://dl.fbaipublicfiles.com/esmatlas/v2023_02/metadata-rc2.sqlite.gz) / [metadata.parquet](https://dl.fbaipublicfiles.com/esmatlas/v2023_02/metadata-rc2.parquet).
You can load the file with pandas: `df = pd.read_parquet('metadata.parquet')`.
The dataframe has 773,846,840 records, and the file size is around 16GB.
This dataframe has 10 columns:
- `id` is the MGnify ID
- `ptm` is the predicted TM score
- `plddt` is the predicted average lddt
- `num_conf` is the number of residues with plddt > 0.7
- `len` is the total residues in the protein
- `is_fragment` indicates whether the protein sequence is identified as a fragment in the MGnify90 sequence database.
- `sequenceChecksum` is the CRC64 hash of the sequence. Can be used for cheaper lookups.
- `esmfold_version` is the version of ESMFold, matching the model accessible as `esm.pretrained.esmfold_v{0,1}`
- `atlas_version` is the Atlas version where this structure first appeared. Note: some of the predictions appearing for the first time in `v0` are also part of Atlas `v2023_02`.
- `sequence_dbs` is the metagenomic source databases where this structure is part of, as `MGnify90_2022_05`, comma-separated if it exists in more than one release, as `MGnify90_2022_05,MGnify90_2023_02`.

Proteins longer than 1280 residues were not folded. To select only the records without missing structures use  `df[~df.plddt.isna()]`.

We recommend using `s5cmd` or `aria2c` to download files (can be installed via conda package manager).
To download any of the structures provided, you can use this `aria2c` command:
```
aria2c --dir=/path/to/download/to --input-file=url-file-provided.txt
```


# Full database

We separate the database by pTM and pLDDT.
This will allow you to download structures based on your needs - for example you can choose to download only the most high confidence structures.
Structures are provided by the bins given in this repo under [v0/full/bins.txt](v0/full/bins.txt) and [v2023_02/full/bins.txt](v2023_02/full/bins.txt).
For example, the foldseek database containing ptm from 0.60 to 0.70 and plddt from 0.80 to 0.90 is named `tm_.60_.70_plddt_.80_.90.DB`.
The data are given as bundles of 500k or 1M structures each.

* **Bulk Predicted structures**: see `{v0,v2023_02}/full/tarballs/`. The URLs for all shards of a single bundle will be in e.g. [v0/full/tarballs/tm_.60_.70_plddt_.80_.90.txt](v0/full/tarballs/tm_.60_.70_plddt_.80_.90.txt). The URLs for all tarballs across bins and their shards is available under `{v0,v2023_02}/full/tarballs.txt`.
* **Foldseek DBs** `{v0,v2023_02}/full/foldseekdb/` and `{v0,v2023_02}/full/foldseekdb.txt`
* **Bulk Embeddings (NEW)** under `{v0,v2023_02}/full/esm2_embeddings/` and `{v0,v2023_02}/full/esm2_embeddings.txt`.
  * Note: Individual LM embeddings can be fetched via the API endpoint `/fetchEmbedding/ESM2/:id` as described on <https://esmatlas.com/about#api>


# High confidence MGnify30 structures

The high confidence MGnify30 structures are built using this procedure:
1. MGnify90 is clustered down to 30% sequence similarity with `mmseqs easy-linclust --kmer-per-seq 100 -cluster-mode 2 --cov-mode 1 -c 0.8`.
1. Structures are filtered to >0.7 pTM and pLDDT
1. Structures are sorted by `pTM * pLDDT` and the best from each cluster is chosen as the representative.


We provide the urls to v0 in [v0/highquality_clust30/tarballs.txt](v0/highquality_clust30/tarballs.txt) and [v0/highquality_clust30/foldseekdb.txt](v0/highquality_clust30/foldseekdb.txt).

As a reminder, please use `aria2c --dir=/path/to/download/to --input-file=v0/highquality_clust30/foldseekdb.txt` to download the foldseek database

The sequences corresponding to the high confidence structures can be downloaded from <https://dl.fbaipublicfiles.com/esmatlas/v0/highquality_clust30/highquality_clust30.fasta>.

# Citation
If you use any of the ESM Metagenomic Atlas data in your work, please cite


```bibtex
@article{lin2022evolutionary,
  title={Evolutionary-scale prediction of atomic level protein structure with a language model},
  author={Lin, Zeming and Akin, Halil and Rao, Roshan and Hie, Brian and Zhu, Zhongkai and Lu, Wenting and Smetanin, Nikita and Verkuil, Robert and Kabeli, Ori and Shmueli, Yaniv and dos Santos Costa, Allan and Fazel-Zarandi, Maryam and Sercu, Tom and Candido, Salvatore and Rives, Alexander},
  year={2022},
  journal={bioRxiv},
  note={bioRxiv 2022.07.20.500902},
  url={https://doi.org/10.1101/2022.07.20.500902},
}
```

## License <a name="license"></a>

ESM Metagenomic Atlas (also referred to as “ESM Metagenomic Structure Atlas” or “ESM Atlas”) data is available under a CC BY 4.0 license for academic and commercial use. Copyright (c) Meta Platforms, Inc. All Rights Reserved. Use of the ESM Metagenomic Atlas data is subject to the Meta Open Source [Terms of Use](https://opensource.fb.com/legal/terms/) and [Privacy Policy](https://opensource.fb.com/legal/privacy/).

If you make use of a structure predicted by ESMFold, please cite the following paper:

```
Zeming Lin, Halil Akin, Roshan Rao, Brian Hie, Zhongkai Zhu, Wenting Lu, Nikita Smetanin, Robert Verkuil, Ori Kabeli, Yaniv Shmueli, Allan dos Santos Costa, Maryam Fazel-Zarandi, Tom Sercu, Sal Candido, Alexander Rives.
Evolutionary-scale prediction of atomic level protein structure with a language model.
Science 379, 6637 (2023).
https://www.science.org/doi/abs/10.1126/science.ade2574
```
