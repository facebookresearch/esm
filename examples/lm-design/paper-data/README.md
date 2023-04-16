Here we provide the data associated with the paper
["Language models generalize beyond natural proteins" (2022)](https://doi.org/10.1101/2022.12.21.521521) by 
Robert Verkuil\*, Ori Kabeli\*, Yilun Du, Basile I. M. Wicky, Lukas F. Milles, Justas Dauparas, David Baker, Sergey Ovchinnikov, Tom Sercu, and Alexander Rives.

## Free Generations (Section 3)
Designs from the Free Generations section of the paper (Section 3) along with their statistics and pdb files can be found at [free_generations_full.db](https://dl.fbaipublicfiles.com/fair-esm/examples/lm_design/free_generations_full.db) and can be loaded using:
```
import pandas as pd; pd.read_sql('free_generations_full', 'sqlite:///free_generations_full.db')
```


## Designs with wetlab validation (Section 2 and 3)
* [data.csv](./data.csv) - Load scalar data in `data.csv` with `pd.read_csv`.
* [data.hdf5](https://dl.fbaipublicfiles.com/fair-esm/examples/lm_design/design_lm_data_2022_v1.hdf5)
For long-form data, download `data.hdf5` from [this link](https://dl.fbaipublicfiles.com/fair-esm/examples/lm_design/design_lm_data_2022_v1.hdf5) and load with `pd.read_hdf`.
```
# Design information
Design ID - {F,G}{0-267} unique identifier for each (LM or AlphaFold) design evaluated. 8 Nan values correspond to 8 ground truth sequences tested.
Experiment Name - Label for the testing pool to which the design / ground-truth sequence belongs. See Supplement; Section 1.6 for a full description of submitted sequences. These pools (minus ground-truth sequences) have experimental results shown in fig. S11.
Design Model - 228x LM, 20x AlphaFold, 20x AF+ngram, 8x Ground Truth.
Target ID - PDB ID of de novo target for all fixed backbone designs, 'Generation' for all free generations.
Sequence - Designed sequence

# In Silico Evaluation
*AlphaFold predicted PDB file - Structure prediction from AlphaFold (5x pTM models, select best by pLDDT -> Amber Relax).
AlphaFold RMSD - (AlphaFold-predicted) RMSD to target backbone for fixed backbone designs, Nan for free generations
AlphaFold pLDDT - (AlphaFold-predicted) Avg pLDDT for the predicted structure

# Experimental Evaluation
# Results from experimental testing.  Final classifications are in the booleans: {Soluble, Success, Success+Monodisperse}.
Total Yield - Actual total soluble yield (in mg) from the 4x1mL prep. (Actual yield is closer to ~2x, we can only inject 1/2 of the total product onto the column.)
yield_per_Leq - Total Yield, adjusted to 1 L of culture equivalent
*Elution Volume (mL) - Array of x-values for plotting of the SEC trace.
*Chromatographic Absorbance at 280nm - Array of y-values for plotting of the SEC trace.
*Elution Volume (mL) (raw) - Raw version, data is not truncated, lengths may differ between rows.
*Chromatographic Absorbance at 280nm (raw) - Raw version, data is not truncated, lengths may differ between rows.
Soluble - Total Yield > 0.05 mg.
Success - Soluble and SEC peak at the expected elution volume.
Success+Monodisperse - SEC peak *only* at the expected elution volume.

# Jackhmmer results
# See Supplement, Section 1.5 for verbose details.
# In short: Summary statistics of Jackhmmer searches (-n 1 --seed 0) of the designed sequence against UniRef90.  Hits that were removed from ESM2's train set were removed from consideration here.  See `.txt` files for ID's of these omitted sequences.
min Jackhmmer E-value - Minimum (best-domain) E-value
max Jackhmmer Seq-id (significant hits only) - Maximum Sequence identity over all significant (best domain E-value < 1) hits.
max Jackhmmer TM-score (top-10 hits only) - Maximum TM-score of the ≈top-10 (by best-domain E-value) hits.  (Purging was applied after top-10, so the number considered may be slightly lower, counts were rarely reduced below 7).

(* denotes long-form data only available in data.hdf5)
```

## `artificial_sequence_purge_ids.txt`
ID's of sequences removed due to being annotateed "artificial sequence" by the UniProt website when `2021_04` was the latest release.

## `uniref90_jackhmmer_purge_ids.txt`
ID's of sequences removed by Jackhmmer search (`-n 1 --seed 0`) of UniRef90 when given the de novo target set as queries.

## Minimal structure projection
A small new model head was constructed on top of ESM2, which is [this linear projection layer](https://dl.fbaipublicfiles.com/fair-esm/examples/lm_design/linear_projection_model.pt).
For a given sequence the projection measures the compatibility of the internal representations of the language model with a structure.
The linear projection layer is automatically downloaded when running the `lm_design` code.

## Reference

If using this work, please cite:
```bibtex
@article{verkuil2022language,
  author={Robert Verkuil\*, Ori Kabeli\*, Yilun Du, Basile I. M. Wicky, Lukas F. Milles, Justas Dauparas, David Baker, Sergey Ovchinnikov, Tom Sercu, and Alexander Rives},
  title={Language models generalize beyond natural proteins},
  year={2022},
  journal={bioRxiv},
  note={bioRxiv 2022.12.21.521521},
  url={https://doi.org/10.1101/2022.12.21.521521},
}
```
