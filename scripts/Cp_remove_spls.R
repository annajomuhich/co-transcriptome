###### Remove poor QC samples - Zucchini
###### May 2025 AJM

library(tidyverse)

######## Load data ================================

#### Specify input paths - MAKE CHANGES HERE
hostrep1_path <- "co-transcriptome/raw_reads/all_samples/ZE1/Host_readcounts.csv"
hostrep2_path <- "co-transcriptome/raw_reads/all_samples/ZE2/Host_readcounts.csv"
hostrep3_path <- "co-transcriptome/raw_reads/all_samples/ZE3/Host_readcounts.csv"
bcinrep1_path <- "co-transcriptome/raw_reads/all_samples/ZE1/Bcin_readcounts.csv"
bcinrep2_path <- "co-transcriptome/raw_reads/all_samples/ZE2/Bcin_readcounts.csv"
bcinrep3_path <- "co-transcriptome/raw_reads/all_samples/ZE3/Bcin_readcounts.csv"
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

# #Remove samples with bad qc/low readcounts
df <- df %>%
	dplyr::select(!Z27_2) %>% #lots of dups
	dplyr::select(!Z18_2) %>% #low seq count
	dplyr::select(!Z39_2) %>% #somewhat low read count and many overrepresented seqs
	dplyr::select(!Z71_2) %>% #low seq count
	dplyr::select(!Z28_2)  #lots of overrepresented seqs

#Remove samples with bad qc/low readcounts
hostrep2 <- hostrep2 %>%
	dplyr::select(!Z27_2) %>% #lots of dups
	dplyr::select(!Z18_2) %>% #low seq count
	dplyr::select(!Z39_2) %>% #somewhat low read count and many overrepresented seqs
	dplyr::select(!Z71_2) %>% #low seq count
	dplyr::select(!Z28_2)  #lots of overrepresented seqs
bcinrep2 <- bcinrep2 %>%
	dplyr::select(!Z27_2) %>% #lots of dups
	dplyr::select(!Z18_2) %>% #low seq count
	dplyr::select(!Z39_2) %>% #somewhat low read count and many overrepresented seqs
	dplyr::select(!Z71_2) %>% #low seq count
	dplyr::select(!Z28_2)  #lots of overrepresented seqs

######## Save new datasets ===================================

dir.create("co-transcriptome/raw_reads/badQC_removed")
setwd("co-transcriptome/raw_reads/badQC_removed/")

dir.create("ZE1")
dir.create("ZE2")
dir.create("ZE3")

hostrep1 %>% write.csv("ZE1/Host_readcounts.csv", row.names = F)
bcinrep1 %>% write.csv("ZE1/Bcin_readcounts.csv", row.names = F)
hostrep2 %>% write.csv("ZE2/Host_readcounts.csv", row.names = F)
bcinrep2 %>% write.csv("ZE2/Bcin_readcounts.csv", row.names = F)
hostrep3 %>% write.csv("ZE3/Host_readcounts.csv", row.names = F)
bcinrep3 %>% write.csv("ZE3/Bcin_readcounts.csv", row.names = F)

setwd("../../..")
