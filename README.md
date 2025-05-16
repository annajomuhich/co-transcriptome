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

## Step 3a: model_means.R - Run Generalized Linear Mixed Model with Negative Binomial on host reads

This script inputs counts, sample IDs, and batch information for one host species and outputs for each gene:

- anova results and variances for each model term
- emmeans and standard errors
- DEG data for the infected term, including log2FCs and p values

Run your sbatch command like this (you don't need to hard-code input files into the script)

`sbatch sbatch_modelmeans.sh /path/to/counts.csv /path/to/sampleIDs.csv /path/to/batch.csv /path/to/output_dir/`

## Step 3b: bcin_genemodel.R - Run Generalized Linear Mixed Model with Negatie Binomial on botrytis reads

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

## Step 4: mr2mods_dataset_setup.R

Prepare dataset for mr2mods pipeline. 

input:
- host_norm_counts_expressed.csv
- bcin_norm_counts_expressed.csv

output:
- normalized.matrix

## Step 5: mr2mods

Then put normalized.matrix through the mr2mods pipeline. (https://github.itap.purdue.edu/jwisecav/mr2mods)

Run your sbatch command like this:

`sbatch sbatch_mr2mods.sh path/to/normalized.matrix path/to/outputdir/`

## Step 6: mr2mods_network_reformat.R

Reformat the output from mr2mods. This puts decay rates together and joins with gene annotations.

output:
- all_clusters.csv

## Step 7a: network_reg.R

Combines with DEG data to assign regulation status to networks. Combines networks from the two datasets.

input:
- all_clusters.csv

output:
- network_reg.csv
- network_reg_sig.csv (filtered for significant networks)

## Step 7b: clustertopGO.R

Runs GO enrichment analysis on clusters. This uses a Fisher test to determine significant overrepresentation of particular GO terms within a network in comparison to the gene universe (i.e. GO terms of all the genes in the genome).

input:
- all_clusters.csv (all detected networks for a given host & botrytis dataset)
- host and botrytis annotation files

output:
- hostcluster_topGO.csv
- bcincluster_topGO.csv

```
sbatch sbatch_clustertopGO.sh \
path/to/host_annotation.csv \   # arg 1
path/to/bcin_annotation.csv \   # arg 2
path/to/all_clusters.csv \      # arg 3
path/to/output_dir/             # arg 4
```

## Step 7c: ext_PCC.R

Extracts PCC values for genes within a network/cluster. Can input a single file with multiple clusters, differentiated by 'Cluster' column. Extracted PCCs can be used for drawing/visualizing networks.

input:
- clusters file (.csv containing select networks you want PCCs for)
    - must have columns `Cluster`, `decay_rate`, and `gene`
- /expr.mr/ (directory from mr2mods, contains all PCCs for all combos of genes in the dataset)
    - note, this is a large directory. I can't store these locally.

output:
- e.g. d5_101_pccs.csv (will output a separate file for each unique cluster in the input file)

```
sbatch sbatch_extpcc.sh \
path/to/cluster_file.csv \
path/to/expr.mr/ \
path/to/output_dir/
```

## Step 8: lesion_expr_model.R

Generates anova tables with variance for models

`lesion_size ~ gene_expression * host`

input:
- cucfab_lsm.csv (lesion sizes)
- bcin_adjusted_emmeans.csv (bcin expression data)
- host_ortho_expressed.csv (expression data of 1:1 orthologs of the two hosts)

output:
- lesion_Bcexpr_anova.csv
- lesion_hostexpr_anova.csv

Note: this only takes about 15 minutes to run so can run it locally. sbatch script is included if needed.

Arguments as follows:

```
sbatch sbatch_lesion_expr.sh \
/path/to/lesionsizes.csv \	#arg 1
/path/to/bcin_counts.csv \	#arg 2
/path/to/host_counts.csv \		#arg 3
/path/to/outputdir/					#arg 4
```

## Step X: FDR.R and DEG_padj.R

FDR correction

Though note that I've since integrated this correction into the above so don't need these

- bcin_genemodel.R
- model_means.R
- lesion_expr_model.R