###### Remove poor QC samples - Cowpea
###### May 2025 AJM

library(tidyverse)
library(ggplot2)
library(edgeR)

######## Load data ================================

#### Specify input paths - MAKE CHANGES HERE
hostrep1_path <- "raw_reads/all_samples/IT1/Host_readcounts.csv"
hostrep2_path <- "raw_reads/all_samples/IT2/Host_readcounts.csv"
hostrep3_path <- "raw_reads/all_samples/IT3/Host_readcounts.csv"
bcinrep1_path <- "raw_reads/all_samples/IT1/Bcin_readcounts.csv"
bcinrep2_path <- "raw_reads/all_samples/IT2/Bcin_readcounts.csv"
bcinrep3_path <- "raw_reads/all_samples/IT3/Bcin_readcounts.csv"
####

#read in data
hostrep1 <- read.csv(hostrep1_path)
hostrep2 <- read.csv(hostrep2_path)
hostrep3 <- read.csv(hostrep3_path)
bcinrep1 <- read.csv(bcinrep1_path)
bcinrep2 <- read.csv(bcinrep2_path)
bcinrep3 <- read.csv(bcinrep3_path)

######## Remove bad samples  ================================

### MAKE CHANGES HERE FOR YOUR DATASET BASED ON MULTIQC and PCAs

#remove bad QC samples
#We have 2 libraries for I27_1, because yield was bad on the first library.
#The one from IT2 worked much better and has ~10X more readcounts. Removing I27_1 from rep 1.
hostrep1 <- hostrep1 %>% select(!I27_1)
bcinrep1 <- bcinrep1 %>% select(!I27_1)
#Remove other samples with bad qc/low readcounts
hostrep1 <- hostrep1 %>%
	dplyr::select(!I73_1) #somewhat low read count and does not cluster with others in PCA
bcinrep1 <- bcinrep1 %>%
	dplyr::select(!I73_1) #somewhat low read count and does not cluster with others in PCA
hostrep3 <- hostrep3 %>%
	dplyr::select(!I10_3) %>% #good seq count, but ton of overrepresented sequences
	dplyr::select(!I5_4)  #good seq count, but ton of overrepresented sequences
bcinrep3 <- bcinrep3 %>%
	dplyr::select(!I10_3) %>% #good seq count, but ton of overrepresented sequences
	dplyr::select(!I5_4)  #good seq count, but ton of overrepresented sequences

#remove low readcount isolates
# List of samples prefixes to remove
prefixes <- c("I21_", "I72_", "I11_", "I47_", "I44_", "I34_", "I7_")
# Combine into a single regex pattern
pattern <- paste0("^(", paste(prefixes, collapse = "|"), ")")
# Function to remove columns whose names start with any of the prefixes
remove_unwanted_cols <- function(df) {
	df[ , !grepl(pattern, names(df))]
}
# Put dataframe in list
dfs <- list(bcinrep1, bcinrep2, bcinrep3,
						hostrep1, hostrep2, hostrep3)
# Apply the function to each dataframe
dfs <- lapply(dfs, remove_unwanted_cols)
#pull dfs back out
bcinrep1 <- dfs[[1]]
bcinrep2 <- dfs[[2]]
bcinrep3 <- dfs[[3]]
hostrep1 <- dfs[[4]]
hostrep2 <- dfs[[5]]
hostrep3 <- dfs[[6]]

######## Save new datasets ===================================

dir.create("raw_reads/badQC_removed")
setwd("raw_reads/badQC_removed/")

dir.create("IT1")
dir.create("IT2")
dir.create("IT3")

hostrep1 %>% write.csv("IT1/Host_readcounts.csv", row.names = F)
bcinrep1 %>% write.csv("IT1/Bcin_readcounts.csv", row.names = F)
hostrep2 %>% write.csv("IT2/Host_readcounts.csv", row.names = F)
bcinrep2 %>% write.csv("IT2/Bcin_readcounts.csv", row.names = F)
hostrep3 %>% write.csv("IT3/Host_readcounts.csv", row.names = F)
bcinrep3 %>% write.csv("IT3/Bcin_readcounts.csv", row.names = F)
