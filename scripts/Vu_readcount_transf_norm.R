########### Readcount transformation and TMM normalization
########### 3 Bioreps of co-transcriptomic data
########### May 2025
########### Anna Jo Muhich

library(tidyverse)
library(ggplot2)
library(edgeR)

######## Load data ================================

#### Specify input paths - MAKE CHANGES HERE
hostrep1_path <- "co-transcriptome/raw_reads/badQC_removed/IT1/Host_readcounts.csv"
hostrep2_path <- "co-transcriptome/raw_reads/badQC_removed/IT2/Host_readcounts.csv"
hostrep3_path <- "co-transcriptome/raw_reads/badQC_removed/IT3/Host_readcounts.csv"
bcinrep1_path <- "co-transcriptome/raw_reads/badQC_removed/IT1/Bcin_readcounts.csv"
bcinrep2_path <- "co-transcriptome/raw_reads/badQC_removed/IT2/Bcin_readcounts.csv"
bcinrep3_path <- "co-transcriptome/raw_reads/badQC_removed/IT3/Bcin_readcounts.csv"
####
	
#read in data
hostrep1 <- read.csv(hostrep1_path)
hostrep2 <- read.csv(hostrep2_path)
hostrep3 <- read.csv(hostrep3_path)
bcinrep1 <- read.csv(bcinrep1_path)
bcinrep2 <- read.csv(bcinrep2_path)
bcinrep3 <- read.csv(bcinrep3_path)

######## Combine reps, total up exon reads, and reformat  ================================

#combine reps
#remove the first column of each successive rep so we don't have repeats of transcript names
host <- cbind(hostrep1, hostrep2[,-1], hostrep3[,-1])
bcin <- cbind(bcinrep1, bcinrep2[,-1], bcinrep3[,-1])

#clean up transcript names
bcin$transcript <- sub("^transcript:", "", bcin$transcript)

#remove noncoding RNA for botrytis
bcin <- bcin %>% filter(!grepl("^E", transcript))

#Collapse exon reads into total reads for each gene

#### MAKE CHANGES HERE FOR YOUR HOST SPECIES GENE ID
host <- host %>%    #Use this one for cowpea (e.g. Vigun10g095650.1.v1.2.exon.3 becomes Vigun10g095650)
	mutate(gene = sub("\\.[0-9]+\\.v1.2.exon.[0-9]+$", "", transcript)) %>%
	dplyr::select(transcript, gene, everything())
# host <- host %>%     #Use this one for common bean (e.g. Phvul.001G000600.2.v2.1.exon.20 becomes Phvul.001G000600)
#   mutate(gene = sub("\\.[0-9]+\\.v2.1.exon.[0-9]+$", "", transcript)) %>%
# 	dplyr::select(transcript, gene, everything())

bcin <- bcin %>%    #This one for Botrytis
	mutate(gene = sub("\\.[0-9]+$", "", transcript)) %>%
	dplyr::select(transcript, gene, everything())

#add together counts within each gene for each sample
host <- host %>%
	group_by(gene) %>%
	summarise(across(all_of(setdiff(names(host), c("gene", "transcript"))), sum))
bcin <- bcin %>%
  group_by(gene) %>%
	summarise(across(all_of(setdiff(names(host), c("gene", "transcript"))), sum))

######## Remove genes with no counts  ================================
#remove genes with no counts across any sample.

# Identify rows (genes) where all counts except the first column are zero
removed_genes_host <- host[rowSums(host[, -1] == 0) == ncol(host) - 1, 1]
# Then filter the dataframe to remove those rows
host <- host[rowSums(host[, -1] == 0) != ncol(host) - 1, ]

# Identify rows (genes) where all counts except the first column are zero
removed_genes_bcin <- bcin[rowSums(bcin[, -1] == 0) == ncol(bcin) - 1, 1]
# Then filter the dataframe to remove those rows
bcin <- bcin[rowSums(bcin[, -1] == 0) != ncol(bcin) - 1, ]

####### Make norm_counts_all for host ========================================
#has 0HAI mock and no further low readcount filtering

dir.create("co-transcriptome/normalized_counts")

### reformat as matrix where colnames are samples, rownames are genes, and cells contain numeric read counts
count_data <- host
#remove spaces from the readcounts
count_data[, -1] <- apply(count_data[, -1], 2, function(x) trimws(x))
#remove gene column
count_data <- count_data[,-1]
#convert readcounts to numeric values
count_data[!is.numeric(count_data)] <- lapply(count_data[!is.numeric(count_data)], as.numeric)
#convert to matrix
count_data <- as.matrix(count_data)
#put in gene rownames
row.names(count_data) <- host$gene
# Create a DGEList object
dge <- DGEList(counts = count_data)
#perform normalization using TMM
dge <- calcNormFactors(dge)
#access normalization factors
normalization_factors <- dge$samples$norm.factors
#extract normalized counts
normalized_counts <- cpm(dge, normalized.lib.sizes = T)
normalized_counts <- as.data.frame(normalized_counts)
normalized_counts <- rownames_to_column(normalized_counts,var = "gene")
#write csv
write.csv(normalized_counts, "co-transcriptome/normalized_counts/host_norm_counts_all.csv", row.names = F, col.names = T)

####### Make norm_counts_all for Botrytis ========================================
#has 0HAI mock and no further low readcount filtering

### reformat as matrix where colnames are samples, rownames are genes, and cells contain numeric read counts
count_data <- bcin
#remove spaces from the readcounts
count_data[, -1] <- apply(count_data[, -1], 2, function(x) trimws(x))
#remove gene column
count_data <- count_data[,-1]
#convert readcounts to numeric values
count_data[!is.numeric(count_data)] <- lapply(count_data[!is.numeric(count_data)], as.numeric)
#convert to matrix
count_data <- as.matrix(count_data)
#put in gene rownames
row.names(count_data) <- bcin$gene
# Create a DGEList object
dge <- DGEList(counts = count_data)
#perform normalization using TMM
dge <- calcNormFactors(dge)
#access normalization factors
normalization_factors <- dge$samples$norm.factors
#extract normalized counts
normalized_counts <- cpm(dge, normalized.lib.sizes = T)
normalized_counts <- as.data.frame(normalized_counts)
normalized_counts <- rownames_to_column(normalized_counts,var = "gene")
#write csv
write.csv(normalized_counts, "co-transcriptome/normalized_counts/bcin_norm_counts_all.csv", row.names = F, col.names = T)

####### Make norm_counts_expressed for host ========================================

#Remove 0HAI Mock samples
#For the purposes of building out the model. Exp design doesn't need the time element
#We can add these Mock reads back in later and compare directly to the 48 HAI mock (if needed)
host <- host %>%
	dplyr::select(!I74_1) %>%
	dplyr::select(!I74_2) %>%
	dplyr::select(!I74_3)

##### reformat as matrix where colnames are samples, rownames are genes, and cells contain numeric read counts
count_data <- host
#remove spaces from the readcounts
count_data[, -1] <- apply(count_data[, -1], 2, function(x) trimws(x))
#remove gene column
count_data <- count_data[,-1]
#convert readcounts to numeric values
count_data[!is.numeric(count_data)] <- lapply(count_data[!is.numeric(count_data)], as.numeric)
#convert to matrix
count_data <- as.matrix(count_data)
#put in gene rownames
row.names(count_data) <- host$gene
# Create a DGEList object
dge <- DGEList(counts = count_data)

#visualize cpm frequencies
cpms <- log2(cpm(dge) + 1)
cpms <- cpms %>% as.data.frame() %>%
	pivot_longer(cols = 1:length(colnames(cpms)),
							 names_to = "sample_ID",
							 values_to = "CPM")
ggplot(cpms, aes(x=CPM)) +
	geom_histogram(binwidth = 0.01, color = "black", fill = "blue", alpha = 0.7) +
	theme_minimal() +
	labs(title = "Histogram of CPM", x = "Log2 Counts per Million (CPM)", y = "Frequency")

#Filter out low-expressed genes
#retaining only genes that have >=1 CPM in ~20% of the samples (44 samples)
keep <- rowSums(cpm(dge) >= 1) >= 44
dge <- dge[keep, ]
# Identify the genes that will be removed
removed_genes_host2 <- rownames(dge)[!keep]

#get a csv of all the removed genes
removed_genes_host2 <- as.data.frame(removed_genes_host2)
removed_genes_host2$removal_reason <- "low readcount"
removed_genes_host2 <- removed_genes_host2 %>% rename(gene = removed_genes_host2)
removed_genes_host$removal_reason <- "no readcount"
removed_genes_host <- rbind(removed_genes_host, removed_genes_host2)
removed_genes_host %>% write.csv("co-transcriptome/normalized_counts/host_removed_genes.csv", row.names = F)

#visualize cpm frequencies
cpms <- log2(cpm(dge) + 1)
cpms <- cpms %>% as.data.frame() %>%
	pivot_longer(cols = 1:length(colnames(cpms)),
							 names_to = "sample_ID",
							 values_to = "CPM")
ggplot(cpms, aes(x=CPM)) +
	geom_histogram(binwidth = 0.01, color = "black", fill = "blue", alpha = 0.7) +
	theme_minimal() +
	labs(title = "Histogram of CPM", x = "Log2 Counts per Million (CPM)", y = "Frequency")

#perform normalization using TMM
dge <- calcNormFactors(dge)
#access normalization factors
normalization_factors <- dge$samples$norm.factors
#extract normalized counts
normalized_counts <- cpm(dge, normalized.lib.sizes = T)
normalized_counts <- as.data.frame(normalized_counts)
normalized_counts <- rownames_to_column(normalized_counts,var = "gene")

#write csv
write.csv(normalized_counts, "co-transcriptome/normalized_counts/host_norm_counts_expressed.csv", row.names = F, col.names = T)

####### Make norm_counts_expressed for Botrytis ========================================

#Remove 0HAI Mock samples
#For the purposes of building out the model. Exp design doesn't need the time element
#We can add these Mock reads back in later and compare directly to the 48 HAI mock (if needed)
bcin <- bcin %>%
	dplyr::select(!I74_1) %>%
	dplyr::select(!I74_2) %>%
	dplyr::select(!I74_3)

##### reformat as matrix where colnames are samples, rownames are genes, and cells contain numeric read counts
count_data <- bcin
#remove spaces from the readcounts
count_data[, -1] <- apply(count_data[, -1], 2, function(x) trimws(x))
#remove gene column
count_data <- count_data[,-1]
#convert readcounts to numeric values
count_data[!is.numeric(count_data)] <- lapply(count_data[!is.numeric(count_data)], as.numeric)
#convert to matrix
count_data <- as.matrix(count_data)
#put in gene rownames
row.names(count_data) <- bcin$gene
# Create a DGEList object
dge <- DGEList(counts = count_data)

#visualize cpm frequencies
cpms <- log2(cpm(dge) + 1)
cpms <- cpms %>% as.data.frame() %>%
	pivot_longer(cols = 1:length(colnames(cpms)),
							 names_to = "sample_ID",
							 values_to = "CPM")
ggplot(cpms, aes(x=CPM)) +
	geom_histogram(binwidth = 0.01, color = "black", fill = "blue", alpha = 0.7) +
	theme_minimal() +
	labs(title = "Histogram of CPM", x = "Log2 Counts per Million (CPM)", y = "Frequency")

#Filter out low-expressed genes
#retaining only genes that have >=1 CPM in ~20% of the samples (44 samples)
keep <- rowSums(cpm(dge) >= 1) >= 44
dge <- dge[keep, ]
# Identify the genes that will be removed
removed_genes_bcin2 <- rownames(dge)[!keep]

#get a csv of all the removed genes
removed_genes_bcin2 <- as.data.frame(removed_genes_bcin2)
removed_genes_bcin2$removal_reason <- "low readcount"
removed_genes_bcin2 <- removed_genes_bcin2 %>% rename(gene = removed_genes_bcin2)
removed_genes_bcin$removal_reason <- "no readcount"
removed_genes_bcin <- rbind(removed_genes_bcin, removed_genes_bcin2)
removed_genes_bcin %>% write.csv("co-transcriptome/normalized_counts/bcin_removed_genes.csv", row.names = F)

#visualize cpm frequencies
cpms <- log2(cpm(dge) + 1)
cpms <- cpms %>% as.data.frame() %>%
	pivot_longer(cols = 1:length(colnames(cpms)),
							 names_to = "sample_ID",
							 values_to = "CPM")
ggplot(cpms, aes(x=CPM)) +
	geom_histogram(binwidth = 0.01, color = "black", fill = "blue", alpha = 0.7) +
	theme_minimal() +
	labs(title = "Histogram of CPM", x = "Log2 Counts per Million (CPM)", y = "Frequency")

#perform normalization using TMM
dge <- calcNormFactors(dge)
#access normalization factors
normalization_factors <- dge$samples$norm.factors
#extract normalized counts
normalized_counts <- cpm(dge, normalized.lib.sizes = T)
normalized_counts <- as.data.frame(normalized_counts)
normalized_counts <- rownames_to_column(normalized_counts,var = "gene")

#write csv
write.csv(normalized_counts, "co-transcriptome/normalized_counts/bcin_norm_counts_expressed.csv", row.names = F, col.names = T)




