###### Remove poor QC samples - Common Bean
###### May 2025 AJM

library(tidyverse)
library(ggplot2)
library(edgeR)

######## Load data ================================

#### Specify input paths - MAKE CHANGES HERE
hostrep1_path <- "co-transcriptome/raw_reads/all_samples/UCC1/Host_readcounts.csv"
hostrep2_path <- "co-transcriptome/raw_reads/all_samples/UCC2/Host_readcounts.csv"
hostrep3_path <- "co-transcriptome/raw_reads/all_samples/UCC3/Host_readcounts.csv"
bcinrep1_path <- "co-transcriptome/raw_reads/all_samples/UCC1/Bcin_readcounts.csv"
bcinrep2_path <- "co-transcriptome/raw_reads/all_samples/UCC2/Bcin_readcounts.csv"
bcinrep3_path <- "co-transcriptome/raw_reads/all_samples/UCC3/Bcin_readcounts.csv"
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

hostrep2 <- hostrep2 %>%
	dplyr::select(!U11_2) %>%
	dplyr::select(!U57_2) %>%
	dplyr::select(!U64_2) %>%
	dplyr::select(!U52_2)
bcinrep2 <- bcinrep2 %>%
	dplyr::select(!U11_2) %>%
	dplyr::select(!U57_2) %>%
	dplyr::select(!U64_2) %>%
	dplyr::select(!U52_2)

######## Save new datasets ===================================

dir.create("co-transcriptome/raw_reads/badQC_removed")
setwd("co-transcriptome/raw_reads/badQC_removed/")

dir.create("UCC1")
dir.create("UCC2")
dir.create("UCC3")

hostrep1 %>% write.csv("UCC1/Host_readcounts.csv", row.names = F)
bcinrep1 %>% write.csv("UCC1/Bcin_readcounts.csv", row.names = F)
hostrep2 %>% write.csv("UCC2/Host_readcounts.csv", row.names = F)
bcinrep2 %>% write.csv("UCC2/Bcin_readcounts.csv", row.names = F)
hostrep3 %>% write.csv("UCC3/Host_readcounts.csv", row.names = F)
bcinrep3 %>% write.csv("UCC3/Bcin_readcounts.csv", row.names = F)

setwd("../../..")
