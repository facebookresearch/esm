# Zero-shot variant prediction with protein language models

This folder contains code and pre-trained weights for the ESM-1v models described in 
[Language models enable zero-shot prediction of the effects of mutations on protein function. (Meier et al. 2021)](https://doi.org/10.1101/2021.07.09.450648).

### Labeling a deep mutational scan with model predictions

Given a deep mutational scan and its associated sequence, the effects of mutations can be predicted using an ensemble of five ESM-1v models:
```
python predict.py \
    --model-location esm1v_t33_650M_UR90S_1 esm1v_t33_650M_UR90S_2 esm1v_t33_650M_UR90S_3 esm1v_t33_650M_UR90S_4 esm1v_t33_650M_UR90S_5 \
    --sequence HPETLVKVKDAEDQLGARVGYIELDLNSGKILESFRPEERFPMMSTFKVLLCGAVLSRVDAGQEQLGRRIHYSQNDLVEYSPVTEKHLTDGMTVRELCSAAITMSDNTAANLLLTTIGGPKELTAFLHNMGDHVTRLDRWEPELNEAIPNDERDTTMPAAMATTLRKLLTGELLTLASRQQLIDWMEADKVAGPLLRSALPAGWFIADKSGAGERGSRGIIAALGPDGKPSRIVVIYTTGSQATMDERNRQIAEIGASLIKHW \
    --dms-input ./data/BLAT_ECOLX_Ranganathan2015.csv \
    --mutation-col mutant \
    --dms-output ./data/BLAT_ECOLX_Ranganathan2015_labeled.csv \
    --offset-idx 24 \
    --scoring-strategy wt-marginals
```

Similarly, one could use the MSA Transformer:
```
python predict.py \
    --model-location esm_msa1b_t12_100M_UR50S \
    --sequence HPETLVKVKDAEDQLGARVGYIELDLNSGKILESFRPEERFPMMSTFKVLLCGAVLSRVDAGQEQLGRRIHYSQNDLVEYSPVTEKHLTDGMTVRELCSAAITMSDNTAANLLLTTIGGPKELTAFLHNMGDHVTRLDRWEPELNEAIPNDERDTTMPAAMATTLRKLLTGELLTLASRQQLIDWMEADKVAGPLLRSALPAGWFIADKSGAGERGSRGIIAALGPDGKPSRIVVIYTTGSQATMDERNRQIAEIGASLIKHW \
    --dms-input ./data/BLAT_ECOLX_Ranganathan2015.csv \
    --mutation-col mutant \
    --dms-output ./data/BLAT_ECOLX_Ranganathan2015_labeled.csv \
    --offset-idx 24 \
    --scoring-strategy masked-marginals \
    --msa-path ./data/BLAT_ECOLX_1_b0.5.a3m
```

### Per-task performance
In `data/` we also release result data files of model predictions of the 41 Deep Mutational Scanning datasets reported in the paper
[Language models enable zero-shot prediction of the effects of mutations on protein function. (Meier et al. 2021)](https://doi.org/10.1101/2021.07.09.450648).

The order of the datasets matches the paper Figure 3 and Figure 8;
the first 10 proteins are validation proteins used during method development, and the next 31 are test proteins.

Note the three different levels of aggregation:
* `raw_df` \[[LINK](https://dl.fbaipublicfiles.com/fair-esm/examples/variant-prediction/data/raw_df.csv)\]: every row contains one mutation for a specific protein (multi-index: ["protein_name", "mutant"]) and columns are different prediction methods.
* `rho_pp` dataframe (rho per protein): contains the main metric absolute value of spearman rho per protein summarizing how good a prediction method performs on that protein. Other two fields `rho_boot_mean, rho_boot_std` are the mean and standard deviation of 20 bootstrapped samples.
* `aggregated_rho`: the performance metrics from `rho_pp` averaged over the proteins from the valid / full / test set. There is a multi-index header `(valid / full / test) x (rho / rho_boot_mean / rho_boot_std)`

