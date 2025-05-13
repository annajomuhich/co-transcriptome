### Modeling transcript effect on lesion size
### May 2025 AJM

##### ============= Define paths ==============
args <- commandArgs(trailingOnly = TRUE)

#assign input files to specific variables
lesion_path <- args[1]		#e.g. cucfab_lsm.csv
bcin_expr_path <- args[2]   # e.g. bcin_adjusted_emmeans.csv
host_expr_path <- args[3] # e.g. host_ortho_expr.csv
output_path <- args[4] 		

# Ensure output directory ends with a slash (for safe concatenation)
if (!grepl("/$", output_path)) {
	output_path <- paste0(output_path, "/")
}

##### ============= Load data and libraries ===== 

library(tidyverse)
library(lme4)

message("loading lesion data from ", paste0(lesion_path))
lesion <- read.csv(lesion_path)
message("loading bcin expression data from ", paste0(bcin_expr_path))
bcin_expr <- read.csv(bcin_expr_path)
message("loading host expression data from ", paste0(host_expr_path))
host_expr <- read.csv(host_expr_path)

##### =============== Reformat lesion data ==================

lesion <- lesion %>% select(Isolate_name, UCC, IT)
lesion <- lesion %>% pivot_longer(cols = !Isolate_name,
																	names_to = "genotype",
																	values_to = "lesion_size")
lesion <- lesion %>% rename(iso_name = Isolate_name)

##### ============= Reformat expression data - Bcin ================

df <- left_join(bcin_expr, lesion, by = c("iso_name", "genotype")) %>%
	mutate(host = if_else(genotype == "UCC", "common bean", "cowpea")) %>%
	select(iso_name, host, lesion_size, everything()) %>%
	select(!genotype) %>%
	filter(iso_name != "Mock_48HAI")

df$host <- as.factor(df$host)

##### ============== Loop - Bcin model =================

#get list of transcripts to loop through
genes <- colnames(df)[grep("^Bcin", colnames(df))]

#set up dataframe with first gene in the list
gene <- genes[1]

message(paste("modeling", gene))

formula <- as.formula(paste("lesion_size ~", gene, "* host"))

model <- lm(formula, data = df) #works
#extract anova table
anova <- car::Anova(model)
anova <- as.data.frame(anova)
anova <- tibble::rownames_to_column(anova, var = "variable")
anova$gene <- gene

# Print diagnostic information for first gene
if(gene == genes[1]) {
	message("ANOVA structure for first gene:")
	print(str(anova))
	message("ANOVA column names:")
	print(names(anova))
}

# calculate variance percentage and clean variable names
anova <- anova %>%
	mutate(variance = `Sum Sq` / sum(`Sum Sq`) * 100) %>%
	mutate(variable = case_when(
		variable == gene ~ "gene_expression",
		variable == "host" ~ "host",
		variable == paste0(gene, ":host") ~ "gene_expression:host",
		TRUE ~ variable)) %>%
	mutate(Rsquared = summary(model)$r.squared)
anova <- anova %>% select(gene, everything())

# Print transformed structure for first gene
if(gene == genes[1]) {
	message("Transformed ANOVA structure:")
	print(str(anova))
	message("Final column names:")
	print(names(anova))
}

anova_all <- anova #initialize df for combined anovas

#define list for failed genes
failed_genes <- list()

#adjust gene list
genes <- genes[-1]

for (gene in genes) {
	message(paste("modeling", gene))
	tryCatch({
		formula <- as.formula(paste("lesion_size ~", gene, "* host"))
		
		model <- lm(formula, data = df)
		#extract anova table
		anova <- car::Anova(model)
		anova <- as.data.frame(anova)
		anova <- tibble::rownames_to_column(anova, var = "variable")
		anova$gene <- gene
		
		# calculate variance percentage and clean variable names
		anova <- anova %>%
			mutate(variance = `Sum Sq` / sum(`Sum Sq`) * 100) %>%
			mutate(variable = case_when(
				variable == gene ~ "gene_expression",
				variable == "host" ~ "host",
				variable == paste0(gene, ":host") ~ "gene_expression:host",
				TRUE ~ variable)) %>%
			mutate(Rsquared = summary(model)$r.squared)
		anova <- dplyr::select(anova, gene, everything())
		
		# append to growing ANOVA dataframe
		anova_all <- rbind(anova_all, anova)
		
	}, error = function(e) {
		message(paste("Error encountered for gene:", gene, "Skipping..."))
		failed_genes <<- c(failed_genes, gene)
	})
}

#Do FDR correction (BH)
# Split data by variable type
anova_split <- split(anova_all, anova_all$variable)
# Apply FDR correction to each variable group
anova_fdr <- lapply(anova_split, function(x) {
  p_values <- x$`Pr..F.`
  x$p_adj <- p.adjust(p_values, method = "BH")
  return(x)})
# Recombine into single dataframe
anova_corrected <- do.call(rbind, anova_fdr)
# Reset row names
rownames(anova_corrected) <- NULL

#write results
dir.create(output_path)

# failed_genes <- as.data.frame(failed_genes) %>% t()
# row.names(failed_genes) <- NULL
# colnames(failed_genes) <- "failed_genes"
failed_genes %>% write.csv(paste0(output_path, "failed_Bcgenes.csv"), row.names = F)

anova_corrected %>% write.csv(paste0(output_path, "lesion_Bcexpr_anova.csv"), row.names = F)

##### ================ Reformat data - Host ==========================

host_expr <- host_expr %>% filter(iso_name != "Mock")

lesion <- lesion %>%
	mutate(host = if_else(genotype == "UCC", "Pv", "Vu")) %>%
	select(!genotype)

df <- left_join(host_expr, lesion, by = c("iso_name", "host")) %>%
	select(iso_name, host, lesion_size, everything())


###### Loop - Pv/Vu 1:1 ortholog model ======================================

#get list of transcripts to loop through
genes <- colnames(df)[grep("^Phvul", colnames(df))]

#set up dataframe with first gene in the list
gene <- genes[1]

message(paste("modeling", gene))

formula <- as.formula(paste("lesion_size ~", gene, "* host"))

model <- lm(formula, data = df) #works
#extract anova table
anova <- car::Anova(model)
anova <- as.data.frame(anova)
anova <- tibble::rownames_to_column(anova, var = "variable")
anova$gene <- gene

# calculate variance percentage and clean variable names
anova <- anova %>%
	mutate(variance = `Sum Sq` / sum(`Sum Sq`) * 100) %>%
	mutate(variable = case_when(
		variable == gene ~ "gene_expression",
		variable == paste0(gene, ":host") ~ "gene_expression:host",
		TRUE ~ variable)) %>%
	mutate(Rsquared = summary(model)$r.squared)
anova <- anova %>% select(gene, everything())

anova_all <- anova #initialize df for combined anovas

#define list for failed genes
failed_genes <- list()

#adjust gene list
genes <- genes[-1]

for (gene in genes) {
	message(paste("modeling", gene))
	tryCatch({
		formula <- as.formula(paste("lesion_size ~", gene, "* host"))
		
		model <- lm(formula, data = df)
		#extract anova table
		anova <- car::Anova(model)
		anova <- as.data.frame(anova)
		anova <- tibble::rownames_to_column(anova, var = "variable")
		anova$gene <- gene
		
		# calculate variance percentage and clean variable names
		anova <- anova %>%
			mutate(variance = `Sum Sq` / sum(`Sum Sq`) * 100) %>%
			mutate(variable = case_when(
				variable == gene ~ "gene_expression",
				variable == "host" ~ "host",
				variable == paste0(gene, ":host") ~ "gene_expression:host",
				TRUE ~ variable)) %>%
			mutate(Rsquared = summary(model)$r.squared)
		anova <- dplyr::select(anova, gene, everything())
		
		# append to growing ANOVA dataframe
		anova_all <- rbind(anova_all, anova)
		
	}, error = function(e) {
		message(paste("Error encountered for gene:", gene, "Skipping..."))
		failed_genes <<- c(failed_genes, gene)
	})
}

#Do FDR correction (BH)
# Split data by variable type
anova_split <- split(anova_all, anova_all$variable)
# Apply FDR correction to each variable group
anova_fdr <- lapply(anova_split, function(x) {
  p_values <- x$`Pr..F.`
  x$p_adj <- p.adjust(p_values, method = "BH")
  return(x)})
# Recombine into single dataframe
anova_corrected <- do.call(rbind, anova_fdr)
# Reset row names
rownames(anova_corrected) <- NULL

failed_genes %>% write.csv(paste0(output_path, "failed_hostgenes.csv"), row.names = F)

anova_all %>% write.csv(paste0(output_path, "lesion_hostexpr_anova.csv"), row.names = F)
