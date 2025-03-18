# co-transcriptome

## model_means.R - Run Generalized Linear Mixed Model with Negative Binomial on host reads

This script inputs counts, sample IDs, and batch information for one host species and outputs for each gene:

- anova results and variances for each model term
- emmeans and standard errors
- DEG data for the infected term, including log2FCs and p values

Run your sbatch command like this (you don't need to hard-code input files into the script)

`sbatch sbatch_modelmeans.sh /path/to/counts.csv /path/to/sampleIDs.csv /path/to/batch.csv /path/to/output_dir/`

## bcin_genemodel.R - Run Generalized Linear Mixed Model with Negatie Binomial on botrytis reads

This script inputs counts, sample IDs, and batch information for multiple host species and outputs for each gene:

- anova results and variances for each model term
- emmeans and standard errors

Currently written for comparing only two host species, can easily adapt the dataset setup to accomodate more.
