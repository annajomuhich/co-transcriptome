# co-transcriptome

## Step 1: remove_spls.R

This is to remove bad QC or outlier samples before further processing/modeling. Use results from MultiQC and PCAs to determine what to remove. Needs to be adjusted per dataset

input: /raw_reads/all_samples/readcounts.csv
output: /raw_reads/badQC_removed/readcounts.csv

## Step 2: readcount_transf_norm.R

This script cleans up the readcount data, performs TMM normalization, and calculates cpm for each organism within each sample. Needs to be adjusted per dataset.

input: /raw_reads/badQC_removed/readcounts.csv
output: /normalized_counts/
- host_norm_counts_all.csv (zero-read genes removed)
- host_norm_counts_expressed.csv (low read count genes removed, 0HAI Mock removed)
- host_removed_genes.csv (list of genes that were removed in this script)
- bcin_norm_counts_all.csv
- bcin_norm_counts_expressed.csv
- bcin_removed_genes.csv

## model_means.R - Run Generalized Linear Mixed Model with Negative Binomial on host reads

This script inputs counts, sample IDs, and batch information for one host species and outputs for each gene:

- anova results and variances for each model term
- emmeans and standard errors
- DEG data for the infected term, including log2FCs and p values

Run your sbatch command like this (you don't need to hard-code input files into the script)

`sbatch sbatch_modelmeans.sh /path/to/counts.csv /path/to/sampleIDs.csv /path/to/batch.csv /path/to/output_dir/`

## bcin_genemodel.R - Run Generalized Linear Mixed Model with Negatie Binomial on botrytis reads

This script inputs counts, sample IDs, and batch information for 2 host species and outputs for each gene:

- anova results and variances for each model term
- emmeans and standard errors
- Botrytis DEGs between the hosts

Currently written for comparing only two host species, can easily adapt the dataset setup to accomodate more.

Run your sbatch command like this (where 1 refers to the first host species data and 2 refers to the second host species data):

```
sbatch sbatch_bcin_genemodel.sh \
/path/to/counts1.csv /path/to/counts2.csv \
/path/to/sampleIDs1.csv /path/to/sampleIDs2.csv \
/path/to/batch1.csv /path/to/batch2.csv \
/path/to/outputdir/
```
