######### clustertopGO - GO enrichment analysis on gene coexpression clusters
######### For plant host and Botrytis co-transcriptome
######### Phaseolus vulgaris (common bean)
######### February 2025
######### AJM & LDF

#Install packages ==============================================
library(topGO)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)

#Define paths (Make changes here) =============================
#input - reference annotation files
host_annot_path <- "input/Pvulgaris_442_v2.1.annotation_info.csv"
bcin_annot_path <- "input/Bcin_Annotations_Full_transcript.csv"
#input - cluster data
clust_path <- "input/all_clusters.csv"
#output - make output directory
dir.create("output")
#output - host topGO result
host_topGO_path <- "output/hostcluster_topGO.csv"
bcin_topGO_path <- "output/bcincluster_topGO.csv"

#GeneID2GO Creation - Host ==========================================
# Load annotation
df_annot <- read_csv(host_annot_path)
# Reformat annotation
#df_annot <- df_annot %>% rename(gene = locusName) #this was only req'd for cowpea annotation
df_annot <- df_annot %>%
  distinct(gene, .keep_all = TRUE) %>% # Removing duplicates based on 'gene'
  dplyr::select(gene, GO) # Selecting only the 'gene' and 'GO' columns
df_annot$GO <- gsub(",", ", ", df_annot$GO) #Separate the GO codes in the GO column by comma + space

#Convert this dataframe into a txt file so that it can be used for topGO
write.table(df_annot, file = "df_GO.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE) 
#load as geneID2GO format for topGO
geneID2GO <- readMappings(file = "df_GO.txt") 
#remove the file, we no longer need it
file.remove(file = "df_GO.txt")

#Gene Cluster setup - Host ======================================
# Load data
clust_all <- read.csv(clust_path)
# Filtering to remove Bcin genes
clust_host <- clust_all[-grep("^Bcin", clust_all$gene), ]

# Removing rows where the same gene shows up twice in one cluster
clust_host <- clust_host[!duplicated(clust_host[, c("cluster_full", "gene")]), ]

# Removing clusters that have <3 host genes
clust_host <- clust_host %>%
	group_by(cluster_full) %>%              # Group by the 'cluster_full' variable
	filter(n() >= 3) %>%                    # Keep only groups with 3 or more occurrences
	ungroup()                               # Ungroup the dataframe

# Selecting only the 'cluster_full' and 'gene' columns
clust_host <- clust_host[, c("cluster_full", "gene")]

# Setup for topGO loop ==========================================
# Get a list of unique cluster names in the "cluster_full" column
cluster_names <- unique(clust_host$cluster_full)

#Creates new empty df table to insert data
allRes_combined <- data.frame(GO.ID = character(),
                              Term = character(),
                              Annotated = numeric(),
                              Significant = numeric(),
                              Expected = numeric(),
                              `Rank in classicFisher` = numeric(),
                              classicFisher = character(),
                              cluster = character(),
                              stringsAsFactors = FALSE)

# Loop to run topGO on each cluster  =========================================================
for (cluster_name in cluster_names) {
  # Filter to only include the current cluster and remove duplicate gene rows
  df_cluster <- clust_host %>%
    filter(cluster_full == cluster_name) %>% distinct(gene, .keep_all = TRUE)
  # Create new geneList 
  geneNames <- names(geneID2GO) #get list of genes in the gene universe
  myInterestingGenes <- df_cluster$gene #get list of genes in the cluster, i.e. our interesting genes
  gene_list <- factor(as.integer(geneNames %in% myInterestingGenes)) #assign a binary integer factor (1 = cluster, 0 = not in cluster)
  names(gene_list) <- geneNames #populate gene list with all the gene names
  # topGO data object
  GOdata <- new("topGOdata",
                ontology = "MF", # Molecular Function (MF)
                allGenes = gene_list, # The binary vector 
                annot = annFUN.gene2GO, # Use the annotation function for gene-to-GO mapping
                gene2GO = geneID2GO) # Your gene-to-GO mapping
  # Enrichment Test
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher") # Fisher test
  #classical enrichment analysis that tests over representation of GO terms within the cluster
  resultFisher
  #Analysis of Results
  #Generates df containing the top topGO terms identified by the elim algorithm
  allRes <- GenTable(GOdata,
                     classicFisher = resultFisher)
  # Add the current cluster name as a new column in the allRes table
  allRes$cluster <- cluster_name
  # Append the current allRes (with cluster name) to the combined table
  allRes_combined <- rbind(allRes_combined, allRes)
  # Print a message confirming the cluster has been processed
  message(paste("Cluster", cluster_name, "processed and added to allRes_combined table"))
}
allRes_combined %>% write.csv(host_topGO_path)


#GeneID2GO Creation - Bcin ===================================================
# Load annotation
df_annot <- read.csv(bcin_annot_path, stringsAsFactors = FALSE)
# Reformat annotation
df_annot <- df_annot %>%
	distinct(X.Gene.ID., .keep_all = TRUE) %>% # Removing duplicates based on 'gene'
	dplyr::select(X.Gene.ID.,X.Computed.GO.Function.IDs.,X.Computed.GO.Process.IDs.) # Selecting only the 'gene' and 'GO' columns
# Rename the X.Gene.ID. column to gene
df_annot <- df_annot %>%
	rename(gene = X.Gene.ID.)
# Combine the two columns into one "GO" column, removing "N/A"
df_annot <- df_annot %>%
	mutate(
		GO = paste(
			ifelse(X.Computed.GO.Function.IDs. != "N/A", X.Computed.GO.Function.IDs., ""),
			ifelse(X.Computed.GO.Process.IDs. != "N/A", X.Computed.GO.Process.IDs., ""),
			sep = ",") %>%
			gsub("(^,|,$|,,)", "", .)) %>%  # Remove leading/trailing commas and duplicate commas
	select(gene, GO)  # Keep only the Gene ID and the new GO column
# Remove rows where the GO column is empty
df_annot <- df_annot %>%
	filter(GO != "")
#Separate the GO codes in the GO column by comma + space
df_annot$GO <- gsub(",", ", ", df_annot$GO) 
#Convert this dataframe into a txt file so that it can be used for topGO
write.table(df_annot, file = "df_GO.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE) 
#load as geneID2GO format for topGO
geneID2GO <- readMappings(file = "df_GO.txt") 
#remove the file, we no longer need it
file.remove(file = "df_GO.txt")

#Gene Cluster setup - Bcin ======================================
# Filtering to remove host genes
clust_bcin <- clust_all[grep("^Bcin", clust_all$gene), ]

# Removing rows where the same gene shows up twice in one cluster
clust_bcin <- clust_bcin[!duplicated(clust_bcin[, c("cluster_full", "gene")]), ]

# Removing clusters that have <3 bcin genes
clust_bcin <- clust_bcin %>%
	group_by(cluster_full) %>%              # Group by the 'cluster_full' variable
	filter(n() >= 3) %>%                    # Keep only groups with 3 or more occurrences
	ungroup()                               # Ungroup the dataframe

# Selecting only the 'cluster_full' and 'gene' columns
clust_bcin <- clust_bcin[, c("cluster_full", "gene")]

# Setup for topGO loop ==========================================
# Get a list of unique cluster names in the "cluster_full" column
cluster_names <- unique(clust_bcin$cluster_full)

#Creates new empty df table to insert data
allRes_combined <- data.frame(GO.ID = character(),
															Term = character(),
															Annotated = numeric(),
															Significant = numeric(),
															Expected = numeric(),
															`Rank in classicFisher` = numeric(),
															classicFisher = character(),
															cluster = character(),
															stringsAsFactors = FALSE)

# Loop to run topGO on each cluster  =========================================================
#Error tolerant loop
for (cluster_name in cluster_names) {
	message(paste("Processing cluster:", cluster_name))
	tryCatch({
		# Filter to only include the current cluster and remove duplicate gene rows
		df_cluster <- clust_bcin %>%
			filter(cluster_full == cluster_name) %>%
			distinct(gene, .keep_all = TRUE)
		# Create new geneList 
		geneNames <- names(geneID2GO) # Get list of genes in the gene universe
		myInterestingGenes <- df_cluster$gene # Get list of genes in the cluster
		gene_list <- factor(as.integer(geneNames %in% myInterestingGenes)) # Assign binary factor
		names(gene_list) <- geneNames # Populate gene list with gene names
		# topGO data object
		GOdata <- new("topGOdata",
									ontology = "MF", # Molecular Function
									allGenes = gene_list, # Binary vector
									annot = annFUN.gene2GO, # Annotation function
									gene2GO = geneID2GO) # Gene-to-GO mapping
		# Enrichment Test
		resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher") 
		# Analysis of Results
		allRes <- GenTable(GOdata, classicFisher = resultFisher)
		# Add the current cluster name as a new column in the allRes table
		allRes$cluster <- cluster_name
		# Append the current allRes (with cluster name) to the combined table
		allRes_combined <- rbind(allRes_combined, allRes)
		# Print a message confirming the cluster has been processed
		message(paste("Cluster", cluster_name, "processed successfully."))
	}, error = function(e) {
		# If an error occurs, print a message and skip this cluster_name
		message(paste("Error processing cluster", cluster_name, "- skipping. Error message:", e$message))
	})
}

allRes_combined %>% write.csv(bcin_topGO_path)
