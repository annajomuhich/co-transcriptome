###### Remove poor QC samples - Cucumber
###### May 2025 AJM

library(tidyverse)

######## Load data ================================

#### Specify input paths - MAKE CHANGES HERE
hostrep1_path <- "co-transcriptome/raw_reads/all_samples/Bos1/Host_readcounts.csv"
hostrep2_path <- "co-transcriptome/raw_reads/all_samples/Bos2/Host_readcounts.csv"
hostrep3_path <- "co-transcriptome/raw_reads/all_samples/Bos3/Host_readcounts.csv"
bcinrep1_path <- "co-transcriptome/raw_reads/all_samples/Bos1/Bcin_readcounts.csv"
bcinrep2_path <- "co-transcriptome/raw_reads/all_samples/Bos2/Bcin_readcounts.csv"
bcinrep3_path <- "co-transcriptome/raw_reads/all_samples/Bos3/Bcin_readcounts.csv"
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

#Remove samples with bad qc/low readcounts
hostrep1 <- hostrep1 %>%
	dplyr::select(!B48_1) %>% #virtually no reads
	dplyr::select(!B7_1)  #low seq count, lots of overrepresented sequences, poor clustering
bcinrep1 <- bcinrep1 %>%
	dplyr::select(!B48_1) %>% #virtually no reads
	dplyr::select(!B7_1)  #low seq count, lots of overrepresented sequences, poor clustering
hostrep2 <- hostrep2 %>%
	dplyr::select(!B14_2) %>% #somewhat low read count and does not cluster with others in PCA
	dplyr::select(!B3_2)  #somewhat low read count and does not cluster with others in PCA
bcinrep2 <- bcinrep2 %>%
	dplyr::select(!B14_2) %>% #somewhat low read count and does not cluster with others in PCA
	dplyr::select(!B3_2)  #somewhat low read count and does not cluster with others in PCA

######## Save new datasets ===================================

dir.create("co-transcriptome/raw_reads/badQC_removed")
setwd("co-transcriptome/raw_reads/badQC_removed/")

dir.create("Bos1")
dir.create("Bos2")
dir.create("Bos3")

hostrep1 %>% write.csv("Bos1/Host_readcounts.csv", row.names = F)
bcinrep1 %>% write.csv("Bos1/Bcin_readcounts.csv", row.names = F)
hostrep2 %>% write.csv("Bos2/Host_readcounts.csv", row.names = F)
bcinrep2 %>% write.csv("Bos2/Bcin_readcounts.csv", row.names = F)
hostrep3 %>% write.csv("Bos3/Host_readcounts.csv", row.names = F)
bcinrep3 %>% write.csv("Bos3/Bcin_readcounts.csv", row.names = F)

setwd("../../..")
