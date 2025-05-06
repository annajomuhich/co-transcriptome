###### Dataset setup for mr2mods from emmeans
###### May 2025

library(tidyverse)
library(stringr)

###### ADJUST PATH =========================
path <- "remove_7isos/normalized_counts/Vu/"

###### Make normalized.matrix ==================================
#load data
host <- read.csv(paste0(path, "host_norm_counts_expressed.csv"))
bcin <- read.csv(paste0(path, "bcin_norm_counts_expressed.csv"))

#bind dataframes together if colnames match
if (all(colnames(host) == colnames(bcin))) {
	df <- rbind(host, bcin)
}

#convert gene column to rownames
df <- column_to_rownames(df, var = "gene")

#need sample ID as column names and genes as rownames for mr2mods!

#write matrix csv
df %>% write.table(paste0(path,'normalized.matrix'), sep="\t", quote = FALSE)
