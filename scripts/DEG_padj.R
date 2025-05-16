### Quick script to join adjusted p values to DEGs
### In the future just put this in the model_means.R script

library(tidyverse)

#Define paths
pv_DEG_path <- "remove_7isos/Pv_hostmodel/DEGs_infected.csv"
vu_DEG_path <- "remove_7isos/Vu_hostmodel/DEGs_infected.csv"
pv_padj_path <- "remove_7isos/Pv_hostmodel/anova_FDR.csv"
vu_padj_path <- "remove_7isos/Vu_hostmodel/anova_FDR.csv"
bc_DEG_path <- "remove_7isos/bcin_genemodel_may2/bcin_DEGs.csv"
bc_padj_path <- "remove_7isos/bcin_genemodel_may2/bcin_anova_FDR.csv"

#read in data
pv_DEG <- read.csv(pv_DEG_path)
vu_DEG <- read.csv(vu_DEG_path)
pv_padj <- read.csv(pv_padj_path)
vu_padj <- read.csv(vu_padj_path)
bc_DEG <- read.csv(bc_DEG_path)
bc_padj <- read.csv(bc_padj_path)


### Pv
#filter adjusted p values to only include infected
pv_padj <- pv_padj %>%
	filter(variable == "infected") %>%
    select(gene, p_adj)

#join adjusted p values to DEGs
pv_DEG <- left_join(pv_DEG, pv_padj, by = "gene") %>%
	select(gene, everything())

#write out DEGs with adjusted p values
pv_DEG %>% write.csv("remove_7isos/Pv_hostmodel/DEGs_infected_padj.csv", row.names = F)

### Vu
#filter adjusted p values to only include infected
vu_padj <- vu_padj %>%
	filter(variable == "infected") %>%
    select(gene, p_adj)

#join adjusted p values to DEGs
vu_DEG <- left_join(vu_DEG, vu_padj, by = "gene") %>%
	select(gene, everything())

#write out DEGs with adjusted p values
vu_DEG %>% write.csv("remove_7isos/Vu_hostmodel/DEGs_infected_padj.csv", row.names = F)

### Bcin
#filter adjusted p values to only include infected
bc_padj <- bc_padj %>%
	filter(variable == "genotype") %>%
    select(gene, p_adj)

#join adjusted p values to DEGs
bc_DEG <- left_join(bc_DEG, bc_padj, by = "gene") %>%
	select(gene, everything())

#write out DEGs with adjusted p values
bc_DEG %>% write.csv("remove_7isos/bcin_genemodel_may2/bcin_DEGs_padj.csv", row.names = F)
