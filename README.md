This repository contains code for preprocessing epigenetic data to learn chromatin states across individuals.

# Usage

Covariate correction reduces the number of chromatin states corresponding to sample-level confounders.  This process should be performed within each chromosome and mark.

To perform covariate correction on the epigenetic data across individuals, the regress_covariates_signal.R script can be used.  This script has 5 positional arguments

1. A uncorrected signal input file formatted for ChromHMM BinarizeBed
2  A table of covariates with covariates in the columns and samples in the rows
3. Output file for corrected signal
4. Comma separated string of marks (mark ids should be in the sample ids)
5. chromosome id (used for reformatting for ChromHMM BinarizeBed)

