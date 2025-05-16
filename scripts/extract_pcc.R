######### Extract PCCs from networks
######### May 2025 AJM

#Load packages
library(tidyverse)

#Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if correct number of arguments provided
if (length(args) != 3) {
	stop("Required arguments: 
		1) cluster data file path
		2) PCCs directory path
		3) output directory path")
}

# Check if output directory path ends with a slash, add if missing
if (!endsWith(args[3], "/")) {
  args[3] <- paste0(args[3], "/")
}

# Define input file paths
cluster_path <- args[1]
pcc_path <- args[2]
output_path <- args[3]

### ------------------- Extract PCCs -------------------------
#Get networks to extract PCCs from
networks <- read.csv(cluster_path)

#make output directory
dir.create(output_path)

#go into directory with PCCs
#farm wd
setwd(pcc_path)

#get list of clusters
clusts <- unique(networks$Cluster)

#Loop through each listed cluster
for (i in 1:length(clusts)) {
	message(paste("Starting Cluster:", clusts[i], sep = " "))
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
	filename <- paste0("d", decay, "_", cluster, "_pccs.csv")

	write.csv(df, paste0(output_path,filename))
}


