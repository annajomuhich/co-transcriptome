##### Network Reformatting
##### IT - Vigna unguiculata
##### December 2024 AJM

# Load required libraries
library(tidyverse)

path <- "./combined_analysis/mr2mods/vubc_networks/"

d5 <- read.csv(paste0(path, "expr.mr_005.modules.csv"))
d10 <- read.csv(paste0(path, "expr.mr_010.modules.csv"))
d25 <- read.csv(paste0(path, "expr.mr_025.modules.csv"))
d50 <- read.csv(paste0(path, "expr.mr_050.modules.csv"))
d100 <- read.csv(paste0(path, "expr.mr_100.modules.csv"))

###LOOP TO REFORMAT SEVERAL DATAFRAMES, 5-100
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

##### Let's assign gene descriptions
#Got Pv gene descriptions from Phytozome
#Got the gene descriptions for Bcin from our google drive https://docs.google.com/spreadsheets/d/1em56HZLFOhC0-J_i1PARy9haOlQbv-M9/edit?usp=drive_link&ouid=113501740728630706254&rtpof=true&sd=true
#read them in:
vu_desc <- read_tsv("annotation/Vunguiculata_540_v1.2.annotation_info.txt")
bcin_desc <- read.csv("annotation/Bcin_Annotations_Full_transcript.csv")
#select columns we need
### NOTE: Vu annotation also has best hits for clamy and rice! For now just looking at arabidopsis
vu_desc <- vu_desc %>% select(locusName, `Best-hit-arabi-name`, `Best-hit-arabi-defline`)
bcin_desc <- bcin_desc %>% select(X.Gene.ID., X.Gene.Name.or.Symbol., X.PFam.Description.)
#rename columns
colnames(vu_desc) <- c("gene", "gene_name", "gene_desc")
colnames(bcin_desc) <- c("gene", "gene_name", "gene_desc")
#remove duplicate genes in the reference
#This just keeps the first time the gene shows up in the annotation
bcin_desc <- bcin_desc[!duplicated(bcin_desc$gene), ]
vu_desc <- vu_desc[!duplicated(vu_desc$gene), ]
#bind vu and bcin gene descriptions together
gene_desc <- rbind(vu_desc, bcin_desc)
#need to remove the decimals from Bcin gene codes before joining
remove_decimal <- function(string) {
	sub("\\.\\d+$", "", string)
}
df$gene <- sapply(df$gene, remove_decimal)

#join gene descriptions onto cluster info
df <- left_join(df, gene_desc, by = "gene")

#make cluster_full column and rearrange
df <- df %>%
	mutate(cluster_full = paste0("d", decay_rate, "_", Cluster)) %>%
	select(cluster_full, everything())

# Filtering to include only clusters that contain both cowpea and botrytis genes
filtered_df <- df %>%
	group_by(cluster_full) %>%
	filter(any(startsWith(gene, "V")) & any(startsWith(gene, "B"))) %>%
	ungroup()

#Saving a copy of both the full clusters and the cross kingdom clusters
dir.create(paste0(path,"reformatted_networks/"))
write.csv(df, paste0(path,"reformatted_networks/all_clusters.csv"), row.names = FALSE)
write.csv(filtered_df, paste0(path, "reformatted_networks/all_vubc_clusters.csv"), row.names = F)


