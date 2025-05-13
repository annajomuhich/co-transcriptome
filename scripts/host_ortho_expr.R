### Combine expression data for 1:1 orthologs between two hosts
### May 2025 AJM

### Define paths ================================
ortho_path <- "../../Fabaceae/fab_orthologs/reformatted_orthologs/PvVu_ortho_1to1.csv"

pv_expr_path <- "./remove_7isos/Pv_hostmodel/adjusted_emmeans.csv"
vu_expr_path <- "./remove_7isos/Vu_hostmodel/adjusted_emmeans.csv"

out_path <- "./remove_7isos/"

### Read in data ===============================
library(tidyverse)

ortho <- read.csv(ortho_path)

pv_expr <- read.csv(pv_expr_path)
vu_expr <- read.csv(vu_expr_path)

### Reformat data ==============================
ortho <- ortho %>%
    mutate(ortho_ID = paste0(Pv, "_", Vu)) %>%
    select(ortho_ID, Pv, Vu)

pv_expr_long <- pv_expr %>%
    select(!infected) %>%
    pivot_longer(cols = -1, names_to = "gene", values_to = "pv_expr")
vu_expr_long <- vu_expr %>%
    select(!infected) %>%
    pivot_longer(cols = -1, names_to = "gene", values_to = "vu_expr")

df <- pv_expr_long %>%
    left_join(ortho, by = c("gene" = "Pv")) %>%
    left_join(vu_expr_long, by = c("Vu" = "gene", "iso_name")) %>%
    rename(Pv = gene) %>%
    select(ortho_ID, Pv, Vu, iso_name, pv_expr, vu_expr) %>%
    filter(!is.na(ortho_ID)) %>%
    filter(!is.na(vu_expr)) %>%
    pivot_longer(cols = c(pv_expr, vu_expr), names_to = "host", values_to = "expr") %>%
    pivot_wider(
    id_cols = c(iso_name, host),
    names_from = ortho_ID,
    values_from = expr) %>%
    mutate(host = if_else(host == "pv_expr", "Pv", "Vu"))

df %>% write.csv(paste0(out_path, "host_ortho_expr.csv"), row.names = FALSE)

#Note that combined dataset has ~10000 genes removed because it only includes expressed genes in both hosts
#only about 8000 genes are left in either host


