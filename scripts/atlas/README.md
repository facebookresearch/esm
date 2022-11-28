# ESM Metagenomic Atlas

The [ESM Metagenomic Atlas](https://esmatlas.com) is a repository of over database of more than 600 million metagenomic protein structures predicted by ESMFold.

Bulk download instructions are available here, as well as foldseek databases available for download.

All structures in the ESM Metagenomic Atlas were predicted with the ESMFold model released as `esm.pretrained.esmfold_v0()`.
We find that protein structures with predicted LDDT > 0.7 and predict TM > 0.7 to be both reasonably well structured and interesting.
Therefore, we provide both the small set of "high quality" metagenomic structures, as well as the full set.
The small set of structures is built from taking a 30% sequence identity clustering of MGnify90, and using the best structure from each cluster.
The best structure is selected using the score pTM * pLDDT.

The high quality structures are around 1TB in size.

The full database is available as PDB structures and is 15TB in size.

We also provide a metadata dataframe: <https://dl.fbaipublicfiles.com/esmatlas/v0/stats.parquet>.
You can load the file with pandas: `df = pd.read_parquet('stats.parquet')`.
The dataframe has length `617051007`, the file size is 6.0GB and has md5 hash `3948a44562b6bd4c184167465eec17de`.
This dataframe has 4 columns:
- `id` is the MGnify ID
- `ptm` is the predicted TM score
- `plddt` is the predicted average lddt
- `num_conf` is the number of residues with plddt > 0.7
- `len` is the total residues in the protein

In parallel with `stats.parquet`, the sequences can be downloaded as fasta file from: <https://dl.fbaipublicfiles.com/esmatlas/v0/full/atlas.fasta>.
The fasta file has `617051007` records matching the stats file, has file size 114GB, and has md5 hash `dc45f4383536c93f9d871facac7cca93`.

We recommend using `s5cmd` or `aria2c` to download files (installable via anaconda).

**To download any of the structures provided, please use this `aria2c` command**
```
aria2c --dir=/path/to/download/to --input-file=url-file-provided.txt
```

NOTE: Predicted Aligned errors will be uploaded at a later date, stay tuned!

# High quality MGnify30 structures

The high quality MGnify30 structures are built using this procedure:
1. MGnify90 is clustered down to 30% sequence similarity with `mmseqs easy-linclust --kmer-per-seq 100 -cluster-mode 2 --cov-mode 1 -c 0.8`.
1. Structures are filtered to >0.7 pTM and pLDDT
1. Structures are sorted by `pTM * pLDDT` and the best from each cluster is chosen as the representative.


The high quality structures can be downloaded from the s3 paths below:

```
s3://dl.fbaipublicfiles.com/esmatlas/v0/highquality_clust30/tarballs/
s3://dl.fbaipublicfiles.com/esmatlas/v0/highquality_clust30/foldseekdb/
s3://dl.fbaipublicfiles.com/esmatlas/v0/highquality_clust30/highquality_clust30.fasta
```

We provide the urls to v0 in [v0/highquality_clust30/tarballs.txt](v0/highquality_clust30/tarballs.txt) and [v0/highquality_clust30/foldseekdb.txt](v0/highquality_clust30/foldseekdb.txt).

As a reminder, please use `aria2c --dir=/path/to/download/to --input-file=v0/highquality_clust30/foldseekdb.txt` to download the foldseek database

# Full database

We separate the database by pTM and pLDDT.
This will allow you to download structures based on your needs - for example you can choose to download only the most high quality structures.
We provide each stratification in `v0/full/tarballs/` and `v0/full/foldseekdb/`
Structures are provided by the bins given in this repo under [v0/full/bins.txt](v0/full/bins.txt).

For example, the foldseek database containing ptm from 0.60 to 0.70 and plddt from 0.80 to 0.90 is named `tm_.60_.70_plddt_.80_.90.DB`.
The PDBs are given as bundles of 1 million structures each.
The URLs for that bundle will be in `v0/full/tarballs/tm_.60_.70_plddt_.80_.90.txt`

We also provide the fasta file for the full dataset:
`s3://dl.fbaipublicfiles.com/esmatlas/v0/full/mgnify90.fasta`
