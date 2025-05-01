########### Readcount transformation and TMM normalization
########### 3 Bioreps of co-transcriptomic data
########### November 2024
########### Anna Jo Muhich

library(tidyverse)
library(ggplot2)
library(edgeR)

#read in data
hostrep1 <- read.csv("IT1/readcounts/Host_readcounts.csv")
hostrep2 <- read.csv("IT2/readcounts/Host_readcounts.csv")
hostrep3 <- read.csv("IT3/readcounts/Host_readcounts.csv")
bcinrep1 <- read.csv("IT1/readcounts/Bcin_readcounts.csv")
bcinrep2 <- read.csv("IT2/readcounts/Bcin_readcounts.csv")
bcinrep3 <- read.csv("IT3/readcounts/Bcin_readcounts.csv")

#We have 2 libraries for I27_1, because yield was bad on the first library.
#Investigating:
I27_1A <- hostrep1 %>% select(transcript, I27_1)
I27_1B <- hostrep2 %>% select(transcript, I27_1)
summarize(I27_1A, sum(I27_1))
summarize(I27_1B, sum(I27_1))
#The one from IT2 worked much better and has ~10X more readcounts. Removing I27_1 from rep 1.
hostrep1 <- hostrep1 %>% select(!I27_1)
bcinrep1 <- bcinrep1 %>% select(!I27_1)

#combine reps
#remove the first column of each successive rep so we don't have repeats of transcript names
host <- cbind(hostrep1, hostrep2[,-1], hostrep3[,-1])
bcin <- cbind(bcinrep1, bcinrep2[,-1], bcinrep3[,-1])

#clean up transcript names
bcin$transcript <- sub("^transcript:", "", bcin$transcript)

#remove noncoding RNA for botrytis
bcin <- bcin %>% filter(!grepl("^E", transcript))

##### Collapse exon reads into total reads for each gene
#make gene column
#Use this one for cowpea (e.g. Vigun10g095650.1.v1.2.exon.3 becomes Vigun10g095650)
host <- host %>%
	mutate(gene = sub("\\.[0-9]+\\.v1.2.exon.[0-9]+$", "", transcript)) %>%
	dplyr::select(transcript, gene, everything())
# #Use this one for common bean (e.g. Phvul.001G000600.2.v2.1.exon.20 becomes Phvul.001G000600)
# host <- host %>%
#   mutate(gene = sub("\\.[0-9]+\\.v2.1.exon.[0-9]+$", "", transcript)) %>%
# 	dplyr::select(transcript, gene, everything())
#This one for Botrytis
bcin <- bcin %>%
	mutate(gene = sub("\\.[0-9]+$", "", transcript)) %>%
	dplyr::select(transcript, gene, everything())

#add together counts within each gene for each sample
#sample IDs start with "I"
host <- host %>%
  group_by(gene) %>%
  summarise(across(starts_with("I"), sum))
bcin <- bcin %>%
  group_by(gene) %>%
  summarise(across(starts_with("I"), sum))

#check colnames match
colnames(host) == colnames(bcin)

#bind readcounts together
df <- rbind(host, bcin)

#Remove samples with bad qc/low readcounts
df <- df %>%
	dplyr::select(!I10_3) %>% #good seq count, but ton of overrepresented sequences
	dplyr::select(!I5_4) %>% #good seq count, but ton of overrepresented sequences
	dplyr::select(!I73_1) #somewhat low read count and does not cluster with others in PCA

#remove genes with no counts across any sample.
#For IT this removed ~1778 genes
df <- df[rowSums(df[, -1] ==0) != ncol(df) - 1, ]

######################
#make norm_counts_all
#has 0HAI mock and no further low readcount filtering

### reformat as matrix where colnames are samples, rownames are genes, and cells contain numeric read counts
count_data <- df
#remove spaces from the readcounts
count_data[, -1] <- apply(count_data[, -1], 2, function(x) trimws(x))
#remove gene column
count_data <- count_data[,-1]
#convert readcounts to numeric values
count_data[!is.numeric(count_data)] <- lapply(count_data[!is.numeric(count_data)], as.numeric)
#convert to matrix
count_data <- as.matrix(count_data)
#put in gene rownames
row.names(count_data) <- df$gene
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
write.csv(normalized_counts, "combined_analysis/data/norm_counts_all.csv", row.names = F, col.names = T)

###########################
#make norm_counts_expressed

#Remove 0HAI Mock samples
#For the purposes of building out the model. Exp design doesn't need the time element
#We can add these Mock reads back in later and compare directly to the 48 HAI mock (if needed)
df <- df %>%
	dplyr::select(!I74_1) %>%
	dplyr::select(!I74_2) %>%
	dplyr::select(!I74_3)

##### reformat as matrix where colnames are samples, rownames are genes, and cells contain numeric read counts
count_data <- df
#remove spaces from the readcounts
count_data[, -1] <- apply(count_data[, -1], 2, function(x) trimws(x))
#remove gene column
count_data <- count_data[,-1]
#convert readcounts to numeric values
count_data[!is.numeric(count_data)] <- lapply(count_data[!is.numeric(count_data)], as.numeric)
#convert to matrix
count_data <- as.matrix(count_data)
#put in gene rownames
row.names(count_data) <- df$gene
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
write.csv(normalized_counts, "combined_analysis/data/norm_counts_expressed.csv", row.names = F, col.names = T)


