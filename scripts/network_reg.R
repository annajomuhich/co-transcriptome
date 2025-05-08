### Assigning regulation status to networks
### May 2025 AJM

library(tidyverse)

### ------- adjust paths ------------------------
pv_net_path <- "remove_7isos/reformatted_networks/Pv/all_clusters.csv"
vu_net_path <- "remove_7isos/reformatted_networks/Vu/all_clusters.csv"

pv_deg_path <- "remove_7isos/Pv_hostmodel/DEGs_infected.csv"
vu_deg_path <- "remove_7isos/Vu_hostmodel/DEGs_infected.csv"

pv_annot_path <- "~/UCDavis/Klieb_Lab/Projects/Fabaceae/Fab_Manuscript/data/gene_descriptions/Pvulgaris_442_v2.1.annotation_info.txt"
vu_annot_path <- "~/UCDavis/Klieb_Lab/Projects/Fabaceae/Fab_Manuscript/data/gene_descriptions/Vunguiculata_540_v1.2.annotation_info.txt"

output_path <- "remove_7isos/reformatted_networks/"

### ------- Prepare Pv Networks --------------
#read in network data
pv_networks <- read.csv(pv_net_path)
#read in data
deg_data <- read.csv(pv_deg_path)
# Add -log10(p-value) and significance column
deg_data <- deg_data %>%
	mutate(
		neg_log10_p = -log10(Pr..Chisq.),
		significance = case_when(
			Pr..Chisq. < 0.01 & log2FC > 1 ~ "Upregulated",
			Pr..Chisq. < 0.01 & log2FC < -1 ~ "Downregulated",
			TRUE ~ "Not Significant"
		)
	)

#prepare annotation file
annot <- read.delim(pv_annot_path)
# annotation file has multiple transcripts for some genes. Just keep the first transcript of each gene:
annot <- annot %>%
	distinct(locusName, .keep_all = TRUE)
#join
deg_data <- left_join(deg_data, annot, join_by(gene == locusName))

#select deg columns
deg_data <- deg_data %>% select(gene, significance)
#join
pv_netreg <- left_join(pv_networks, deg_data, by = "gene")
#make new column to differentiate host clusters
pv_netreg$clustfull_host <- paste0("pv_",pv_netreg$cluster_full)
pv_netreg <- pv_netreg %>% select(clustfull_host, everything())


### -------- Prepare Vu networks -------------
#read in network data
vu_networks <- read.csv(vu_net_path)
#read in data
deg_data <- read.csv(vu_deg_path)
# Add -log10(p-value) and significance column
deg_data <- deg_data %>%
	mutate(
		neg_log10_p = -log10(Pr..Chisq.),
		significance = case_when(
			Pr..Chisq. < 0.01 & log2FC > 1 ~ "Upregulated",
			Pr..Chisq. < 0.01 & log2FC < -1 ~ "Downregulated",
			TRUE ~ "Not Significant"
		)
	)

#prepare annotation file
annot <- read.delim(vu_annot_path)
# annotation file has multiple transcripts for some genes. Just keep the first transcript of each gene:
annot <- annot %>%
	distinct(locusName, .keep_all = TRUE)
#join
deg_data <- left_join(deg_data, annot, join_by(gene == locusName))

#select deg columns
deg_data <- deg_data %>% select(gene, significance)
#join
vu_netreg <- left_join(vu_networks, deg_data, by = "gene")
#make new column to differentiate host clusters
vu_netreg$clustfull_host <- paste0("vu_",vu_netreg$cluster_full)
vu_netreg <- vu_netreg %>% select(clustfull_host, everything())


### -------- make combined dataframes  -------------
#Make combined dataframe
netreg <- rbind(pv_netreg, vu_netreg)

#Change significance column name for clarity
netreg <- netreg %>% rename(infected_response = significance)

#save the dataframe
netreg %>% write.csv(paste0(output_path,"network_reg.csv"), row.names = F)

#Also would like a version that only has significant networks
netreg_sig <- netreg %>%
	filter(P.value < 0.05)
netreg_sig %>% write.csv(paste0(output_path,"network_reg_sig.csv"), row.names = F)

