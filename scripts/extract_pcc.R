library(tidyverse)

#farm WD for parallel jobs
setwd("~/kliebengrp/ext_pcc/ext_pcc_1")

### MAKE CHANGES HERE:
#define path to cluster list
cluster_path <- "SM_clusters/pv_CHS_ofinterest.csv"
#define path to PCCs
pcc_path <- "~/kliebengrp/pvbc_coex/output/expr.mr"
#define SM tag for output (e.g. terp, bot, boa, etc)
SM <- "chs"

### ------------------- Extract PCCs -------------------------
#Get networks to extract PCCs from
networks <- read.csv(cluster_path)

#go into directory with PCCs
#farm wd
setwd(pcc_path)

#get list of clusters
clusts <- unique(networks$Cluster)

#Loop through each listed cluster
for (i in 1:length(clusts)) {
	#set cluster variable
	cluster <- clusts[i]
	#filter the networks by that cluster
	network_filt <- networks %>%
		filter(Cluster == cluster)
	#Join gene codes
	#genes <- left_join(network_filt, gene_codes, by = "gene")
	#get gene list (in gene codes)
	genes <- network_filt$gene
	#remove duplicates (any will break the code)
	genes <- unique(genes)
	#get decay rate of the cluster
	decay <- unique(network_filt$decay_rate)
	#make empty dataframe
	df <- data.frame(Source = character(), Target = character(), PCC = numeric(),
									 stringsAsFactors = FALSE)
	#Print status
	print(paste("Starting Cluster:", cluster, "Decay Rate:", decay, sep = " "))
	
	# Loop through each gene and create edges to all other genes
	# Took about an hour for ~400x400 genes
	for (i in 1:length(genes)) {
		source_gene <- genes[i]
		print(paste(i,"Getting PCC values for Source:",source_gene, sep = " "))
		# Generate targets for the current source gene
		target_genes <- genes[-i]  # Exclude the current source gene
		
		# Create a data frame for the current source gene and its targets
		source_targets <- data.frame(Source = rep(source_gene, length(target_genes)),
																 Target = target_genes,
																 stringsAsFactors = FALSE)
		# Read PCC values for each target gene and extract the PCC value for the current source-target pair
		for (j in 1:nrow(source_targets)) {
			target_gene <- source_targets$Target[j]
			# Read the PCC values file
			pcc_data <- read_tsv(source_gene, col_names = F, show_col_types = F)
			
			# Extract the PCC value for the current source-target pair
			pcc_value <- pcc_data$X3[pcc_data$X1 == target_gene]
			
			# Store the PCC value in the dataframe
			source_targets$PCC[j] <- pcc_value
		}
		# Append the source-target data frame to the main dataframe
		df <- rbind(df, source_targets)
	}
	#specify filename, path, and save
	filename <- paste0("d", decay, "_", cluster, "_", SM, ".csv")
	path <- paste0("~/kliebengrp/ext_pcc/ext_pcc_1/output/")
	write.csv(df, paste0(path,filename))
}


