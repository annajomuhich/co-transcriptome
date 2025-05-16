### Assigning regulation status to networks
### May 2025 AJM
### May 14 2025 AJM: adding DEG data for bcin between hosts

library(tidyverse)

### ------- adjust paths ------------------------
pv_net_path <- "remove_7isos/reformatted_networks/Pv/all_clusters.csv"
vu_net_path <- "remove_7isos/reformatted_networks/Vu/all_clusters.csv"

pv_deg_path <- "remove_7isos/Pv_hostmodel/DEGs_infected_padj.csv"
vu_deg_path <- "remove_7isos/Vu_hostmodel/DEGs_infected_padj.csv"
bc_deg_path <- "remove_7isos/bcin_genemodel_may2/bcin_DEGs_padj.csv"

pv_annot_path <- "~/UCDavis/Klieb_Lab/Projects/Fabaceae/Fab_Manuscript/data/gene_descriptions/Pvulgaris_442_v2.1.annotation_info.txt"
vu_annot_path <- "~/UCDavis/Klieb_Lab/Projects/Fabaceae/Fab_Manuscript/data/gene_descriptions/Vunguiculata_540_v1.2.annotation_info.txt"
bc_annot_path <- "~/UCDavis/Klieb_Lab/Projects/Fabaceae/Fab_Manuscript/data/gene_descriptions/Bcin_Annotations_Full_transcript.csv"

output_path <- "remove_7isos/reformatted_networks/"

### ------- Prepare Pv Networks --------------
#read in network data
pv_networks <- read.csv(pv_net_path)
#read in data
pv_deg_data <- read.csv(pv_deg_path)
bc_deg_data <- read.csv(bc_deg_path)
#remove unnecessary p values
pv_deg_data <- pv_deg_data %>% select(!c(ttest_pvalue, Pr..Chisq.))
bc_deg_data <- bc_deg_data %>% select(!p.value)
# Add -log10(p-value) and significance column
pv_deg_data <- pv_deg_data %>%
	mutate(
		neg_log10_p = -log10(p_adj),
		significance = case_when(
			p_adj < 0.05 & log2FC > 1 ~ "Inf Upregulated",
			p_adj < 0.05 & log2FC < -1 ~ "Inf Downregulated",
			TRUE ~ "Not Significant"
		)
	)
bc_deg_data <- bc_deg_data %>%
	mutate(
		neg_log10_p = -log10(p_adj),
		significance = case_when(
			p_adj < 0.05 & log2FC > 1 ~ "Pv Upregulated",
			p_adj < 0.05 & log2FC < -1 ~ "Vu Upregulated",
			TRUE ~ "Not Significant"
		)
	)

#prepare annotation file - host
pv_annot <- read.delim(pv_annot_path)
# annotation file has multiple transcripts for some genes. Just keep the first transcript of each gene:
pv_annot <- pv_annot %>%
	distinct(locusName, .keep_all = TRUE)
#select cols of interest
pv_annot <- pv_annot %>% select(locusName, arabi.symbol, arabi.defline)
#rename columns
colnames(pv_annot) <- c("gene", "gene_name", "gene_desc")
#join
pv_deg_data <- left_join(pv_deg_data, pv_annot, by = "gene")

#prepare annotation file - bcin
bc_annot <- read.csv(bc_annot_path)
# annotation file has multiple transcripts for some genes. Just keep the first transcript of each gene:
bc_annot <- bc_annot %>%
	distinct(X.Gene.ID., .keep_all = TRUE)
#select cols of interest
bc_annot <- bc_annot %>% select(X.Gene.ID., X.Gene.Name.or.Symbol., X.PFam.Description.)
#rename columns
colnames(bc_annot) <- c("gene", "gene_name", "gene_desc")
#join
bc_deg_data <- left_join(bc_deg_data, bc_annot, by = "gene")

#combine host and bcin together
deg_data <- rbind(pv_deg_data, bc_deg_data)

#select deg columns
deg_data <- deg_data %>% select(gene, significance)
#join
pv_netreg <- left_join(pv_networks, deg_data, by = "gene")
#make new column to differentiate host clusters
pv_netreg$clustfull_host <- paste0("pv_",pv_netreg$cluster_full)
pv_netreg <- pv_netreg %>% select(clustfull_host, everything())

# ### -------- Prepare Vu networks -------------
#read in network data
vu_networks <- read.csv(vu_net_path)
#read in data
vu_deg_data <- read.csv(vu_deg_path)
bc_deg_data <- read.csv(bc_deg_path)
#remove unnecessary p values
vu_deg_data <- vu_deg_data %>% select(!c(ttest_pvalue, Pr..Chisq.))
bc_deg_data <- bc_deg_data %>% select(!p.value)
# Add -log10(p-value) and significance column
vu_deg_data <- vu_deg_data %>%
	mutate(
		neg_log10_p = -log10(p_adj),
		significance = case_when(
			p_adj < 0.05 & log2FC > 1 ~ "Inf Upregulated",
			p_adj < 0.05 & log2FC < -1 ~ "Inf Downregulated",
			TRUE ~ "Not Significant"
		)
	)
bc_deg_data <- bc_deg_data %>%
	mutate(
		neg_log10_p = -log10(p_adj),
		significance = case_when(
			p_adj < 0.05 & log2FC > 1 ~ "Pv Upregulated",
			p_adj < 0.05 & log2FC < -1 ~ "Vu Upregulated",
			TRUE ~ "Not Significant"
		)
	)

#prepare annotation file - host
vu_annot <- read.delim(vu_annot_path)
# annotation file has multiple transcripts for some genes. Just keep the first transcript of each gene:
vu_annot <- vu_annot %>%
	distinct(locusName, .keep_all = TRUE)
#select cols of interest
vu_annot <- vu_annot %>% select(locusName, Best.hit.arabi.name, Best.hit.arabi.defline)
#rename columns
colnames(vu_annot) <- c("gene", "gene_name", "gene_desc")
#join
vu_deg_data <- left_join(vu_deg_data, vu_annot, by = "gene")

#prepare annotation file - bcin
bc_annot <- read.csv(bc_annot_path)
# annotation file has multiple transcripts for some genes. Just keep the first transcript of each gene:
bc_annot <- bc_annot %>%
	distinct(X.Gene.ID., .keep_all = TRUE)
#select cols of interest
bc_annot <- bc_annot %>% select(X.Gene.ID., X.Gene.Name.or.Symbol., X.PFam.Description.)
#rename columns
colnames(bc_annot) <- c("gene", "gene_name", "gene_desc")
#join
bc_deg_data <- left_join(bc_deg_data, bc_annot, by = "gene")

#combine host and bcin together
deg_data <- rbind(vu_deg_data, bc_deg_data)

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
netreg <- netreg %>% rename(DEG_status = significance)

#save the dataframe
netreg %>% write.csv(paste0(output_path,"network_reg.csv"), row.names = F)

#Also would like a version that only has significant networks
netreg_sig <- netreg %>%
	filter(P.value < 0.05)
netreg_sig %>% write.csv(paste0(output_path,"network_reg_sig.csv"), row.names = F)

