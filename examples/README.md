# What's in this directory

* The notebooks are introduced and summarized in `../README.md`
* `data/some_proteins.fasta` and its smaller version, `data/few_proteins.fasta` are a random selection of UniRef50 sequences used in the second example of `../README.md`
* `data/1a3a_1_A.a3m`, `data/1xcr_1_A.a3m`, `data/5ahw_1_A.a3m` are MSAs distributed with trRosetta, used in `contact_prediction.ipynb`
* `data/P62593.fasta` is introduced and used in `sup_variant_prediction.ipynb`
* Example MSAs genereated in the same way as the MSAs used for MSA Transformer pre-training:
  - `data/UniRef50_E9K9Y4.a3m`, `data/UniRef50_UPI0003108055.a3m`, `data/UniRef50_UPI0003674933.a3m`, from the same sequences as trRosetta:
  `data/hhblits_uniclust_2017_10_1a3a_1_A.a3m`, `data/hhblits_uniclust_2017_10_1xcr_1_A.a3m`, `data/hhblits_uniclust_2017_10_5ahw_1_A.a3m`.
  - Generated with: `hhblits -i UniRef50_$id.fas -oa3m UniRef50_$id.a3m -n 3 -d /uniclust30_2017_10/uniclust30_2017_10`.
* `esm2_infer_fairscale_fsdp_cpu_offloading.py` shows how to load the ESM-2 15B model with Fairscale's FSDP's CPU offloading capability on a single GPU