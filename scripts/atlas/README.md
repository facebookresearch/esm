# ESM Metagenomic Atlas

The [ESM Metagenomic Atlas](https://esmatlas.com) is a repository of over database of more than 600 million metagenomic protein structures predicted by ESMFold.
See our [blog post](https://ai.facebook.com/blog/protein-folding-esmfold-metagenomics/) to learn more.

The first `v0` version of the Atlas was released on November 1st 2022,
corresponding to the sequences in the `2022_05` release of the MGnify protein database described [here](https://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/2022_05/README.txt).
An update `v2023_02` was released on March 17th 2023, corresponding to the sequences in the `2023_02` release of the MGnify protein database described [here](https://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/2022_11/).

Bulk download instructions are available here, as well as foldseek databases available for download.

The structures in the ESM Metagenomic Atlas were predicted with `esm.pretrained.esmfold_v0()` for Atlas `v0`, while Atlas `v2023_02` used `esm.pretrained.esmfold_v1()`.
We find that protein structures with predicted LDDT > 0.7 and predict TM > 0.7 to be both reasonably well structured and interesting.
Therefore, we provide both the small set of "high confidence" metagenomic structures from Atlas `v0`, as well as the full set.
The small set of structures is built from taking a 30% sequence identity clustering of MGnify90 `2022_05`, and using the best structure from each cluster.
The best structure is selected using the score pTM * pLDDT.

The high confidence structures are around 1TB in size.

The full database is available as PDB structures and is 15TB in size for `v0`.

We also provide a metadata dataframe: <https://dl.fbaipublicfiles.com/esmatlas/v2023_02/metadata.parquet>.
You can load the file with pandas: `df = pd.read_parquet('metadata.parquet')`.
The dataframe has length `TODO`, the file size is around 25GB and has md5 hash `TODO`.
This dataframe has TODO columns:
- `id` is the MGnify ID
- `ptm` is the predicted TM score
- `plddt` is the predicted average lddt
- `num_conf` is the number of residues with plddt > 0.7
- `len` is the total residues in the protein

In parallel with `v0/stats.parquet`, the sequences can be downloaded as fasta file from: <https://dl.fbaipublicfiles.com/esmatlas/v0/full/atlas.fasta>.
The fasta file has `617051007` records matching the stats file, has file size 114GB, and has md5 hash `dc45f4383536c93f9d871facac7cca93`.

We recommend using `s5cmd` or `aria2c` to download files (installable via anaconda).

**To download any of the structures provided, please use this `aria2c` command**
```
aria2c --dir=/path/to/download/to --input-file=url-file-provided.txt
```

NOTE: Predicted Aligned errors will be uploaded at a later date, stay tuned!

# High confidence MGnify30 structures

The high confidence MGnify30 structures are built using this procedure:
1. MGnify90 is clustered down to 30% sequence similarity with `mmseqs easy-linclust --kmer-per-seq 100 -cluster-mode 2 --cov-mode 1 -c 0.8`.
1. Structures are filtered to >0.7 pTM and pLDDT
1. Structures are sorted by `pTM * pLDDT` and the best from each cluster is chosen as the representative.


The high confidence structures can be downloaded from the s3 paths below:

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
We provide each stratification in `v0/full/tarballs/`,  `v0/full/foldseekdb/`, and `v2023_02/full/tarballs/`.
Structures are provided by the bins given in this repo under [v0/full/bins.txt](v0/full/bins.txt) and TODO [v2023_02/full/bins.txt](v2023_02/full/bins.txt).

For example, the foldseek database containing ptm from 0.60 to 0.70 and plddt from 0.80 to 0.90 is named `tm_.60_.70_plddt_.80_.90.DB`.
The PDBs are given as bundles of 500k or 1M structures each.
The URLs for that bundle will be in `v0/full/tarballs/tm_.60_.70_plddt_.80_.90.txt`


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
