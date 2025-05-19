######### clustertopGO - GO enrichment analysis on gene coexpression clusters
######### For plant host and Botrytis co-transcriptome
######### Phaseolus vulgaris (common bean)
######### May 2025
######### AJM & LDF

#Install packages ==============================================
#BiocManager::install("topGO")

#Load packages ==============================================
library(topGO)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)

#Define paths =============================
# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if correct number of arguments provided
if (length(args) != 4) {
  stop("Required arguments: 
       1) host annotation file path
       2) bcin annotation file path 
       3) cluster data file path
       4) output directory path")
}

# Check if output directory path ends with a slash, add if missing
if (!endsWith(args[4], "/")) {
  args[4] <- paste0(args[4], "/")
}

# Input - reference annotation files
host_annot_path <- args[1]
bcin_annot_path <- args[2]
# Input - cluster data
clust_path <- args[3]
output_path <- args[4]
# Output directory
dir.create(args[4])

#GeneID2GO Creation - Host ==========================================
# Load annotation
if (endsWith(host_annot_path, ".csv")) {
  df_annot <- read_csv(host_annot_path)
} else if (endsWith(host_annot_path, ".txt")) {
  df_annot <- read_delim(host_annot_path, delim="\t")
} else {
  stop("Host annotation file must end in .csv or .txt")
}
# Reformat annotation
# Check if 'locusName' exists and rename to 'gene' if it does
if("locusName" %in% colnames(df_annot) && !("gene" %in% colnames(df_annot))) {
  df_annot <- df_annot %>% rename(gene = locusName)
}
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

# Function to run topGO analysis on clusters
run_topGO_analysis <- function(clust_data, geneID2GO, output_prefix) {
  # Get unique cluster names
  cluster_names <- unique(clust_data$cluster_full)
  
  # #subset for testing REMOVE
  # cluster_names <- cluster_names[1:10]
  
  # Create empty dataframe for results
  allRes_combined <- data.frame(GO.ID = character(),
                              Term = character(),
                              Annotated = numeric(),
                              Significant = numeric(),
                              Expected = numeric(),
                              `Rank in classicFisher` = numeric(),
                              classicFisher = character(),
                              cluster = character(),
                              stringsAsFactors = FALSE)
  
  # Loop to run topGO on each cluster
  for (cluster_name in cluster_names) {
    message(paste("Processing cluster:", cluster_name))
    tryCatch({
      # Filter to only include the current cluster and remove duplicate gene rows
      df_cluster <- clust_data %>%
        filter(cluster_full == cluster_name) %>%
        distinct(gene, .keep_all = TRUE)
      
      # Create new geneList 
      geneNames <- names(geneID2GO)
      myInterestingGenes <- df_cluster$gene
      gene_list <- factor(as.integer(geneNames %in% myInterestingGenes))
      names(gene_list) <- geneNames
      
      # topGO data object
      suppressMessages(
        GOdata <- new("topGOdata",
                     ontology = "MF",
                     allGenes = gene_list,
                     annot = annFUN.gene2GO,
                     gene2GO = geneID2GO)
      )
      
      # Enrichment Test
      suppressMessages(
        resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
      )
      
      # Analysis of Results
      allRes <- GenTable(GOdata, classicFisher = resultFisher)
      allRes$cluster <- cluster_name
      allRes_combined <- rbind(allRes_combined, allRes)
      
      message(paste("Cluster", cluster_name, "processed successfully."))
    }, error = function(e) {
      message(paste("Error processing cluster", cluster_name, "- skipping. Error message:", e$message))
    })
  }
  
  # Write full results
  message(paste("Writing full results to", paste0(output_path, output_prefix, "_clustertopGO_all.csv")))
  allRes_combined %>% write.csv(paste0(output_path, output_prefix, "_clustertopGO_all.csv"), row.names = F)
  
  # Create and return summary
  summary <- allRes_combined %>%
    group_by(cluster) %>%
    slice_min(order_by = classicFisher, n = 1, with_ties = FALSE) %>%
    select(cluster, Term, classicFisher) %>%
    rename(!!paste0(output_prefix, "_GO_FisherTest") := classicFisher,
           !!paste0(output_prefix, "_GO_Term") := Term)
  
  return(summary)
}

# Run topGO analysis for host genes
host_summary <- run_topGO_analysis(clust_host, geneID2GO, "host")

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

# Run topGO analysis for Bcin genes
bcin_summary <- run_topGO_analysis(clust_bcin, geneID2GO, "bcin")

### Write combined host and bcin summary as one output 
summary <- full_join(host_summary, bcin_summary, by = "cluster") %>%
	arrange(cluster)
summary %>% write.csv(paste0(output_path, "clustertopGO_summary.csv"), row.names = F)
