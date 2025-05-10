### Modeling transcript effect on lesion size
### May 2025 AJM

##### ============= Define paths ==============
args <- commandArgs(trailingOnly = TRUE)

#assign input files to specific variables
lesion_path <- args[1]		#e.g. cucfab_lsm.csv
bcin_expr_path <- args[2]   # e.g. bcin_adjusted_emmeans.csv
pv_expr_path <- args[3] # e.g. adjusted_emmeans.csv
vu_expr_path <- args[4] # e.g. adjusted_emmeans.csv
output_path <- args[5] 		

# Ensure output directory ends with a slash (for safe concatenation)
if (!grepl("/$", output_dir)) {
	output_dir <- paste0(output_dir, "/")
}

##### ============= Load data and libraries ===== 

library(tidyverse)
library(lme4)

lesion <- read.csv(lesion_path)
bcin_expr <- read.csv(bcin_expr_path)
pv_expr <- read.csv(pv_expr_path)
vu_expr <- read.csv(vu_expr_path)

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

#subset for testing
genes <- genes[1:10]

#set up dataframe with first gene in the list
gene <- genes[1]

message(paste("modeling", gene))

#Note: this one wasn't working, I think because you can't assign continuous variables as random effects
#formula <- as.formula(lesion_size ~ (1|Bcin01g00010) + host + (Bcin01g00010 * host))

formula <- as.formula(paste("lesion_size ~", gene, "* host"))

model <- lm(formula, data = df) #works
#extract anova table
anova <- print(car::Anova(model)) #save to object. this also displays on output
anova <- rownames_to_column(anova, var = "variable") #get column for variable category

#gather variance data
fixed_var <- (diag(vcov(model)))
# fixed_var is a named numeric vector from diag(vcov(model))
percent_var <- 100 * fixed_var / sum(fixed_var)
# Put into a data frame for clarity
fixed_var_df <- data.frame(
	term = names(fixed_var),
	variance = percent_var)
fixed_var_df <- fixed_var_df %>%
	mutate(variable = gsub("hostcowpea", "host", term)) %>%
	select(!term)

#join to anova data
anova <- full_join(anova, fixed_var_df, by = "variable")
anova <- anova %>%
	mutate(variable = gsub(paste(gene), "gene_expression", variable))
anova$gene <- gene
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
		
		# extract ANOVA table
		anova <- car::Anova(model)  # don't wrap in print(), this returns the object
		anova <- tibble::rownames_to_column(as.data.frame(anova), var = "variable")
		
		# gather variance data
		fixed_var <- diag(vcov(model))
		percent_var <- 100 * fixed_var / sum(fixed_var)
		fixed_var_df <- data.frame(
			term = names(fixed_var),
			variance = percent_var
		) %>%
			mutate(variable = gsub("hostcowpea", "host", term)) %>%
			select(-term)
		
		# join to ANOVA data
		anova <- dplyr::full_join(anova, fixed_var_df, by = "variable") %>%
			mutate(variable = gsub(gene, "gene_expression", variable))
		
		anova$gene <- gene
		anova <- dplyr::select(anova, gene, everything())
		
		# append to growing ANOVA dataframe
		anova_all <- rbind(anova_all, anova)
		
	}, error = function(e) {
		message(paste("Error encountered for gene:", gene, "Skipping..."))
		failed_genes <<- c(failed_genes, gene)
	})
}

#write results
dir.create(output_path)

failed_genes <- as.data.frame(failed_genes) %>% t() %>% rename(failed_genes = V1)
row.names(failed_genes) <- NULL
colnames(failed_genes) <- "failed_genes"
failed_genes %>% write.csv(paste0(output_path, "failed_Bcgenes.csv"), row.names = F)

anova_all %>% write.csv(paste0(output_path, "lesion_Bcexpr_anova.csv"), row.names = F)

##### ================ Reformat data - Host ==========================

pv_expr <- pv_expr %>% select(-infected) %>% filter(iso_name != "Mock")
vu_expr <- vu_expr %>% select(-infected) %>% filter(iso_name != "Mock")

pv_lesion <- lesion %>% filter(genotype == "UCC")
vu_lesion <- lesion %>% filter(genotype == "IT")

pv_df <- left_join(pv_lesion, pv_expr, by = "iso_name") %>%
	mutate(host = if_else(genotype == "UCC", "common bean", NA)) %>%
	select(iso_name, host, lesion_size, everything()) %>%
	select(!genotype) %>%
	filter(iso_name != "Mock_48HAI")

vu_df <- left_join(vu_lesion, vu_expr, by = "iso_name") %>%
	mutate(host = if_else(genotype == "IT", "cowpea", NA)) %>%
	select(iso_name, host, lesion_size, everything()) %>%
	select(!genotype) %>%
	filter(iso_name != "Mock_48HAI")

###### Loop - Pv model ======================================

df <- pv_df

#get list of transcripts to loop through
genes <- colnames(df)[grep("^Phvul", colnames(df))]

#subset for testing
genes <- genes[1:10]

#set up dataframe with first gene in the list
gene <- genes[1]

message(paste("modeling", gene))

formula <- as.formula(paste("lesion_size ~", gene))

model <- lm(formula, data = df) #works
#extract anova table
anova <- print(car::Anova(model)) #save to object. this also displays on output
anova <- rownames_to_column(anova, var = "variable") #get column for variable category

#gather variance data
fixed_var <- (diag(vcov(model)))
# fixed_var is a named numeric vector from diag(vcov(model))
percent_var <- 100 * fixed_var / sum(fixed_var)
# Put into a data frame for clarity
fixed_var_df <- data.frame(
	variable = names(fixed_var),
	variance = percent_var)

#join to anova data
anova <- full_join(anova, fixed_var_df, by = "variable")
anova <- anova %>%
	mutate(variable = gsub(paste(gene), "gene_expression", variable))
anova$gene <- gene
anova <- anova %>% select(gene, everything())

anova_all <- anova #initialize df for combined anovas

#define list for failed genes
failed_genes <- list()

#adjust gene list
genes <- genes[-1]

for (gene in genes) {
	message(paste("modeling", gene))
	tryCatch({
		formula <- as.formula(paste("lesion_size ~", gene))
		
		model <- lm(formula, data = df) #works
		#extract anova table
		anova <- print(car::Anova(model)) #save to object. this also displays on output
		anova <- rownames_to_column(anova, var = "variable") #get column for variable category
		
		#gather variance data
		fixed_var <- (diag(vcov(model)))
		# fixed_var is a named numeric vector from diag(vcov(model))
		percent_var <- 100 * fixed_var / sum(fixed_var)
		# Put into a data frame for clarity
		fixed_var_df <- data.frame(
			variable = names(fixed_var),
			variance = percent_var)
		
		#join to anova data
		anova <- full_join(anova, fixed_var_df, by = "variable")
		anova <- anova %>%
			mutate(variable = gsub(paste(gene), "gene_expression", variable))
		anova$gene <- gene
		anova <- anova %>% select(gene, everything())
		
		# append to growing ANOVA dataframe
		anova_all <- rbind(anova_all, anova)
		
	}, error = function(e) {
		message(paste("Error encountered for gene:", gene, "Skipping..."))
		failed_genes <<- c(failed_genes, gene)
	})
}

failed_genes <- as.data.frame(failed_genes) %>% t() %>% rename(failed_genes = V1)
row.names(failed_genes) <- NULL
colnames(failed_genes) <- "failed_genes"
failed_genes %>% write.csv(paste0(output_path, "failed_Pvgenes.csv"), row.names = F)

anova_all %>% write.csv(paste0(output_path, "lesion_Pvexpr_anova.csv"), row.names = F)

###### Loop - Vu model ======================================

df <- vu_df

#get list of transcripts to loop through
genes <- colnames(df)[grep("^Vigun", colnames(df))]

#subset for testing
genes <- genes[1:10]

#set up dataframe with first gene in the list
gene <- genes[1]

message(paste("modeling", gene))

formula <- as.formula(paste("lesion_size ~", gene))

model <- lm(formula, data = df) #works
#extract anova table
anova <- print(car::Anova(model)) #save to object. this also displays on output
anova <- rownames_to_column(anova, var = "variable") #get column for variable category

#gather variance data
fixed_var <- (diag(vcov(model)))
# fixed_var is a named numeric vector from diag(vcov(model))
percent_var <- 100 * fixed_var / sum(fixed_var)
# Put into a data frame for clarity
fixed_var_df <- data.frame(
	variable = names(fixed_var),
	variance = percent_var)

#join to anova data
anova <- full_join(anova, fixed_var_df, by = "variable")
anova <- anova %>%
	mutate(variable = gsub(paste(gene), "gene_expression", variable))
anova$gene <- gene
anova <- anova %>% select(gene, everything())

anova_all <- anova #initialize df for combined anovas

#define list for failed genes
failed_genes <- list()

#adjust gene list
genes <- genes[-1]

for (gene in genes) {
	message(paste("modeling", gene))
	tryCatch({
		formula <- as.formula(paste("lesion_size ~", gene))
		
		model <- lm(formula, data = df) #works
		#extract anova table
		anova <- print(car::Anova(model)) #save to object. this also displays on output
		anova <- rownames_to_column(anova, var = "variable") #get column for variable category
		
		#gather variance data
		fixed_var <- (diag(vcov(model)))
		# fixed_var is a named numeric vector from diag(vcov(model))
		percent_var <- 100 * fixed_var / sum(fixed_var)
		# Put into a data frame for clarity
		fixed_var_df <- data.frame(
			variable = names(fixed_var),
			variance = percent_var)
		
		#join to anova data
		anova <- full_join(anova, fixed_var_df, by = "variable")
		anova <- anova %>%
			mutate(variable = gsub(paste(gene), "gene_expression", variable))
		anova$gene <- gene
		anova <- anova %>% select(gene, everything())
		
		# append to growing ANOVA dataframe
		anova_all <- rbind(anova_all, anova)
		
	}, error = function(e) {
		message(paste("Error encountered for gene:", gene, "Skipping..."))
		failed_genes <<- c(failed_genes, gene)
	})
}

failed_genes <- as.data.frame(failed_genes) %>% t() %>% rename(failed_genes = V1)
row.names(failed_genes) <- NULL
colnames(failed_genes) <- "failed_genes"
failed_genes %>% write.csv(paste0(output_path, "failed_Vugenes.csv"), row.names = F)

anova_all %>% write.csv(paste0(output_path, "lesion_Vuexpr_anova.csv"), row.names = F)

