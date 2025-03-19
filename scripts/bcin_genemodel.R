##### Bcin gene model - Fabales
##### For comparison of Bcin gene expression across multiple hosts
##### March 2025 AJM

### Load packages and data ===================
library(tidyverse)
library(glmmTMB)
library(emmeans)
library(car)

#load count file for each host
pv <- read.csv("~/bcin_genemodel/input/Pv/norm_counts_expressed.csv")
vu <- read.csv("~/bcin_genemodel/input/Vu/norm_counts_expressed.csv")

#load sample key for each host
pv_sample_key <- read.csv("~/bcin_genemodel/input/Pv/ucc_rnaseq_sampleIDs.csv")
vu_sample_key <- read.csv("~/bcin_genemodel/input/Vu/it_rnaseq2_sampleIDs.csv")

#load seq batch info for each host
pv_seq_batch <- read.csv("~/bcin_genemodel/input/Pv/full_sequenced_batches.csv")
vu_seq_batch <- read.csv("~/bcin_genemodel/input/Vu/full_sequenced_batches.csv")


### Prepare exp design info =======================
#reformat seq batch
#add a species identifier to seq batch number (bc seq batch 1 for each species is a separate batch)
pv_seq_batch <- pv_seq_batch %>%
	mutate(seq_batch = paste0("pv", seq_batch))
vu_seq_batch <- vu_seq_batch %>%
	mutate(seq_batch = paste0("vu", seq_batch))
colnames(pv_seq_batch) == colnames(vu_seq_batch) #check colnames match
seq_batch <- rbind(pv_seq_batch, vu_seq_batch)
#reformat sample key
colnames(pv_sample_key) == colnames(vu_sample_key) #check colnames match
sample_key <- rbind(pv_sample_key, vu_sample_key)
#join seqbatch into sample key
sample_key <- left_join(sample_key, seq_batch, by = "sample_ID")


### Reformat count data ======================
#remove host genes
pv <- pv %>% filter(startsWith(gene, "B"))
vu <- vu %>% filter(startsWith(gene, "B"))
#join dfs, but only retain genes that are present in both datasets
df <- inner_join(pv, vu, by = "gene") #this dropped ~350-400 genes for Fabales
#pivot longer
df_long <- pivot_longer(df,
												cols = !gene,
												names_to = "sample_ID",
												values_to = "count")
#check for NA values in count data (needs to be FALSE)
any(is.na(df_long$count))
#join with sample key
df_long <- left_join(df_long, sample_key, by = "sample_ID")
#check for NA values in seqbatch (needs to be FALSE)
any(is.na(df_long$seq_batch))
#make infected column? Or remove mocks?

#convert categorical variables to factor
head(df_long)
df_long$genotype <- as.factor(df_long$genotype)
df_long$iso_name <- as.factor(df_long$iso_name)
df_long$iso_number <- as.factor(df_long$iso_number)
df_long$tray <- as.factor(df_long$rep)
df_long$leaf <- as.factor(df_long$leaf)
df_long$inoc_position <- as.factor(df_long$inoc_position)
df_long$plant <- as.factor(df_long$plant)
df_long$seq_batch <- as.factor(df_long$seq_batch)
df_long$gene <- as.factor(df_long$gene)
#df_long$infected <- as.factor(df_long$infected)
head(df_long) 
#need one column for each gene, reformat:
df <- pivot_wider(df_long,
									names_from = gene,
									values_from = count)

### Run ANOVA ======================
#Below is an error-intolerant loop
#get list of genes to iterate through
genes <- unique(df_long$gene)
genes <- as.character(genes)
genes <- genes[320:330] #subset for testing

#set up dataframe with first gene in the list
gene <- genes[1]

print(paste(""))
print(paste("modeling", gene))
print(paste(""))

#write formula
formula <- as.formula(paste(gene, "~",
														"genotype +",
														"iso_name +",
														"genotype * iso_name +",
														"(1|tray) +",
														"(1|seq_batch)"))
#build model
model <- glmmTMB(formula,
								 data = df,
								 family = nbinom2) %>% suppressMessages()

#extract anova table
anova <- print(car::Anova(model)) #save to object. this also displays on output
anova <- rownames_to_column(anova, var = "variable") #get column for variable category
anova$gene <- paste(gene) #add gene column

#gather variance data
# Extract variance of fixed effects (diagonal of covariance matrix)
fixed_var <- as.data.frame(diag(vcov(model)$cond))
fixed_var <- rownames_to_column(fixed_var, var = "term")
colnames(fixed_var)[2] <- "variance" #change column name
var_sums <- fixed_var %>%
	summarise(genotype = sum(variance[grepl("^genotype", term) & !grepl(":", term)], na.rm = TRUE),
						iso_name = sum(variance[grepl("^iso_name", term)], na.rm = TRUE),
						`genotype:iso_name` = sum(variance[grepl("^genotype", term) & grepl(":", term)], na.rm = TRUE),
						intercept = sum(variance[grepl("Intercept", term)]),
						tot_var = sum(genotype, iso_name, `genotype:iso_name`, intercept))
#calculate percent variance for each term
var_sums <- var_sums %>%
	mutate(genotype = genotype/tot_var,
				 iso_name = iso_name/tot_var,
				 `genotype:iso_name` = `genotype:iso_name`/tot_var,
				 intercept = intercept/tot_var)
#add column for gene
var_sums <- var_sums %>%
	select(!tot_var)
#pivot longer
var_sums <- var_sums %>%
	pivot_longer(cols = everything(),
							 names_to = "variable",
							 values_to = "variance")
#join to anova data
anova <- full_join(anova, var_sums, by = "variable")
anova$gene <- gene
anova <- anova %>% select(gene, everything())

anova_all <- anova #initialize df for combined anovas

#calculate emmeans
emmresult <- emmeans(model, ~ iso_name + genotype)
#get emm summary
emmsummary <- summary(emmresult)

#set up dataframes for emmeans and SE
emm_df <- as.data.frame(emmsummary)
SE_df <- as.data.frame(emmsummary)
emm_df <- emm_df %>%
	dplyr::select(iso_name, genotype, emmean)
SE_df <- SE_df %>%
	dplyr::select(iso_name, genotype, SE)
colnames(emm_df)[3] <- gene
colnames(SE_df)[3] <- gene

#adjust gene list
genes <- genes[-1]

#modified error-tolerant loop:
# Initialize an empty list to track genes with errors
failed_genes <- list()

for (gene in genes) {
	print(paste(""))
	print(paste("modeling", gene))
	print(paste(""))
	
	# Try-catch block for the modeling process
	tryCatch({
		# Write formula
		formula <- as.formula(paste(gene, "~",
																"genotype +",
																"iso_name +",
																"genotype * iso_name +",
																"(1|tray) +",
																"(1|seq_batch)"))
		# Build model
		model <- glmmTMB(formula, 
										 data = df, 
										 family = nbinom2) %>% suppressMessages()

		#extract anova table
		anova <- print(car::Anova(model)) #save to object
		anova <- rownames_to_column(anova, var = "variable") #get column for variable category
		anova$gene <- paste(gene) #add gene column
		
		#gather variance data
		# Extract variance of fixed effects (diagonal of covariance matrix)
		fixed_var <- as.data.frame(diag(vcov(model)$cond))
		fixed_var <- rownames_to_column(fixed_var, var = "term")
		colnames(fixed_var)[2] <- "variance" #change column name
		var_sums <- fixed_var %>%
			summarise(genotype = sum(variance[grepl("^genotype", term) & !grepl(":", term)], na.rm = TRUE),
								iso_name = sum(variance[grepl("^iso_name", term)], na.rm = TRUE),
								`genotype:iso_name` = sum(variance[grepl("^genotype", term) & grepl(":", term)], na.rm = TRUE),
								intercept = sum(variance[grepl("Intercept", term)]),
								tot_var = sum(genotype, iso_name, `genotype:iso_name`, intercept))
		#calculate percent variance for each term
		var_sums <- var_sums %>%
			mutate(genotype = genotype/tot_var,
						 iso_name = iso_name/tot_var,
						 `genotype:iso_name` = `genotype:iso_name`/tot_var,
						 intercept = intercept/tot_var)
		#add column for gene
		var_sums <- var_sums %>%
			select(!tot_var)
		#pivot longer
		var_sums <- var_sums %>%
			pivot_longer(cols = everything(),
									 names_to = "variable",
									 values_to = "variance")
		#join to anova data
		anova <- full_join(anova, var_sums, by = "variable")
		anova$gene <- gene
		anova <- anova %>% select(gene, everything())
		
		anova_all <- rbind(anova_all, anova) #add this genes anova to the growing anova df
		
		# Calculate emmeans
		emmresult <- emmeans(model, ~ iso_name + genotype)
		# Get emm summary
		emmsummary <- summary(emmresult)
		# Join with existing dataframes
		emm_toadd <- emmsummary$emm %>% as.data.frame()
		colnames(emm_toadd) <- gene
		emm_df <- cbind(emm_df, emm_toadd)
		SE_toadd <- emmsummary$SE %>% as.data.frame()
		colnames(SE_toadd) <- gene
		SE_df <- cbind(SE_df, SE_toadd)
	}, error = function(e) {
		# Handle error: add gene to failed list and print a message
		print(paste("Error encountered for gene:", gene, "Skipping..."))
		failed_genes <<- c(failed_genes, gene) # Append to global list
	})
}

#write out results
dir.create("~/bcin_genemodel/output")
write.csv(emm_df, "~/bcin_genemodel/output/bcin_adjusted_emmeans.csv", row.names = F)
write.csv(SE_df, "~/bcin_genemodel/output/bcin_adjusted_SE.csv", row.names = F)
write.csv(anova_all, "~/bcin_genemodel/output/bcin_anova.csv", row.names = F)
write.csv(failed_genes, "~/bcin_genemodel/output/failed_genes.csv", row.names = FALSE)
