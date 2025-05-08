##### Network Reformatting from mr2mods output
##### May 2025 AJM

### =============== ADJUST INPUTS HERE ===============

input_path <- "./remove_7isos/raw_networks/Pv/" #dir where your indiv decay rate csvs are
output_path <- "./remove_7isos/reformatted_networks/Pv/"

hostref_path <- "~/UCDavis/Klieb_Lab/Projects/Fabaceae/Fab_Manuscript/data/gene_descriptions/Pvulgaris_442_v2.1.annotation_info.txt"
bcinref_path <- "~/UCDavis/Klieb_Lab/Projects/Fabaceae/Fab_Manuscript/data/gene_descriptions/Bcin_Annotations_Full_transcript.csv"

#define column names to retain from host gene descriptions
host_geneID <- "locusName" #for Pv
host_genename <- "arabi-symbol" #for Pv
host_genedesc <- "arabi-defline" #for Pv

### =============== Load libraries and data ===============
library(tidyverse)

d5 <- read.csv(paste0(input_path, "expr.mr_005.modules.csv"))
d10 <- read.csv(paste0(input_path, "expr.mr_010.modules.csv"))
d25 <- read.csv(paste0(input_path, "expr.mr_025.modules.csv"))
d50 <- read.csv(paste0(input_path, "expr.mr_050.modules.csv"))
d100 <- read.csv(paste0(input_path, "expr.mr_100.modules.csv"))

### =============== Loop to reformat by decay rate ===============

# Create a vector of dataframe names
dataframe_names <- c("d5", "d10", "d25", "d50", "d100")

for (name in dataframe_names) {
	# Get the current dataframe
	current_df <- get(name)
	# Lengthen the df so that each member of each cluster has its own row
	current_df <- current_df %>%
		mutate(gene = strsplit(Members, " ")) %>%
		unnest(gene) %>%
		select(-Members)
	# Convert gene code columns to character class
	current_df$gene <- as.character(current_df$gene)
	# Update the current dataframe in the global environment
	assign(name, current_df)
}

#Add columns for decay rates and bind dfs together
d5$decay_rate <- "5"
d10$decay_rate <- "10"
d25$decay_rate <- "25"
d50$decay_rate <- "50"
d100$decay_rate <- "100"
df <- rbind(d5, d10, d25, d50, d100)

### =============== Assign annotations ===============
#Got Pv gene descriptions from Phytozome
#Got the gene descriptions for Bcin from our google drive https://docs.google.com/spreadsheets/d/1em56HZLFOhC0-J_i1PARy9haOlQbv-M9/edit?usp=drive_link&ouid=113501740728630706254&rtpof=true&sd=true
#read them in:
host_desc <- read_tsv(hostref_path)
bcin_desc <- read.csv(bcinref_path)
#select columns we need
### NOTE: host annotation also has best hits for clamy and rice! For now just looking at arabidopsis
host_desc <- host_desc %>% select(host_geneID, host_genename, host_genedesc)
bcin_desc <- bcin_desc %>% select(X.Gene.ID., X.Gene.Name.or.Symbol., X.PFam.Description.)
#rename columns
colnames(host_desc) <- c("gene", "gene_name", "gene_desc")
colnames(bcin_desc) <- c("gene", "gene_name", "gene_desc")
#remove duplicate genes in the reference
#This just keeps the first time the gene shows up in the annotation
bcin_desc <- bcin_desc[!duplicated(bcin_desc$gene), ]
host_desc <- host_desc[!duplicated(host_desc$gene), ]
#bind host and bcin gene descriptions together
gene_desc <- rbind(host_desc, bcin_desc)
#need to remove the decimals from Bcin gene codes before joining
	remove_decimal <- function(string) {
		sub("\\.\\d+$", "", string)
	}
df$gene <- sapply(df$gene, remove_decimal)
#join gene descriptions onto cluster info
df <- left_join(df, gene_desc, by = "gene")

### =============== Reformat and save ===============

#make cluster_full column and rearrange
df <- df %>%
	mutate(cluster_full = paste0("d", decay_rate, "_", Cluster)) %>%
	select(cluster_full, everything())

# # Filtering to include only clusters that contain both cowpea and botrytis genes
# filtered_df <- df %>%
# 	group_by(cluster_full) %>%
# 	filter(any(startsWith(gene, "V")) & any(startsWith(gene, "B"))) %>%
# 	ungroup()

#Saving a copy of both the full clusters and the cross kingdom clusters
dir.create(output_path, recursive = T)
write.csv(df, paste0(output_path,"all_clusters.csv"), row.names = FALSE)
#write.csv(filtered_df, paste0(path, "reformatted_networks/all_hostbc_clusters.csv"), row.names = F)


