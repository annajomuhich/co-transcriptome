### Generalized Linear Mixed Model with Negative Binomial
### Calculation of model-corrected means
### ANOVA result extraction for each gene
### DEG collection for infected vs Mock
### AJM March 2025

# install.packages("TMB")
# install.packages("glmmTMB")

library(tidyverse)
library(glmmTMB)
library(emmeans)
library(car) #for glmmTMB anova

#Note about packages:
#The lsmeans package is being deprecated and further development will take place in its successor, emmeans.
#Users may use emmeans in almost exactly the same way as lsmeans, but a few function names and internal details are changed.

args <- commandArgs(trailingOnly = TRUE)

#assign input files to specific variables
counts_file <- args[1]		# Counts file (norm_counts_expressed.csv)
sampleIDs_file <- args[2] # Sample ID file (e.g. it_rnaseq2_sampleIDs.csv)
batch_file <- args[3] 		# batch list file (full_sequenced_batches.csv)
output_dir <- args[4]     # Output directory

# Ensure output directory ends with a slash (for safe concatenation)
if (!grepl("/$", output_dir)) {
	output_dir <- paste0(output_dir, "/")
}

#read in input files
df <- read.csv(counts_file, header = T)
sample_key <- read.csv(sampleIDs_file, header = T)
seq_batch <- read.csv(batch_file, header = T)

#reformat data so each row is a single observation
df_long <- pivot_longer(df,
												cols = !gene,
												names_to = "sample_ID",
												values_to = "count")

#check for NA values in count data (needs to be FALSE)
any(is.na(df_long$count))

##### join data with other categorical variables

#join data
df_long <- left_join(df_long, sample_key, by = "sample_ID")
df_long <- left_join(df_long, seq_batch, by = "sample_ID")

#make column for 'infected'
df_long <- df_long %>%
	mutate(infected = if_else(iso_name == "Mock_48HAI",
				 									 true = "no",
				 									 false = "yes"))

#now can collapse iso_names for mock samples so they are the same (just "Mock")
df_long <- df_long %>%
	mutate(iso_name = ifelse(iso_name %in% c("Mock_0HAI", "Mock_48HAI"),
													 "Mock",
													 iso_name))

###### finish formatting for model
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
df_long$infected <- as.factor(df_long$infected)
head(df_long) 

#need one column for each gene, reformat:
df <- pivot_wider(df_long,
						names_from = gene,
						values_from = count)

########
# Gene by gene ANOVA loop

#Below is an error-intolerant loop
#get list of genes to iterate through
genes <- unique(df_long$gene)
genes <- as.character(genes)
genes <- genes[14904:14909] #subset for testing

#set up dataframe with first gene in the list
gene <- genes[1]

print(paste(""))
print(paste("modeling", gene))
print(paste(""))

#write formula
formula <- as.formula(paste(gene, "~",
														"infected +",
														"infected/iso_name +",
														"tray +",
														"seq_batch"))
#build model
#note this runs a natural log link
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
	summarise(tray = sum(variance[grepl("^tray", term)]),
						infected = sum(variance[grepl("infectedyes$", term)]),
						seq_batch = sum(variance[grepl("^seq_batch", term)]),
						`infected:iso_name` = sum(variance[grepl("^infectedno:|^infectedyes:", term)], na.rm = T),
						intercept = sum(variance[grepl("Intercept", term)]),
						tot_var = sum(tray, infected, seq_batch, `infected:iso_name`, intercept))
#calculate percent variance for each term
var_sums <- var_sums %>%
	mutate(tray = tray/tot_var,
				 infected = infected/tot_var,
				 seq_batch = seq_batch/tot_var,
				 `infected:iso_name` = `infected:iso_name`/tot_var,
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
emmresult <- emmeans(model, specs = "iso_name")
#get emm summary
emmsummary <- summary(emmresult)
#convert natural log estimates to log2
# math: ln(x) = log2(x) / log2(e)
emmsummary$emmean <- emmsummary$emmean / log(2)
emmsummary$SE <- emmsummary$SE / log(2)

#set up dataframes for emmeans and SE
emm_df <- as.data.frame(emmsummary)
SE_df <- as.data.frame(emmsummary)
emm_df <- emm_df %>%
	dplyr::select(iso_name, infected, emmean)
SE_df <- SE_df %>%
	dplyr::select(iso_name, infected, SE)
colnames(emm_df)[3] <- gene
colnames(SE_df)[3] <- gene

#gather DEG infected data
#get infected emmeans. note this is natural log scale
inf_emmresult <- emmeans(model, specs = "infected")
#compute pairwise contrasts for infected DEGs
#revpairwise will give infectedyes - infectedno, so the change from mock to infected
inf_DEG <- contrast(inf_emmresult, method = "revpairwise")
inf_DEG <- summary(inf_DEG)
#convert to log2 scale
inf_DEG$estimate <- inf_DEG$estimate / log(2)
inf_DEG$SE <- inf_DEG$SE / log(2)
#change 'estimate' column to 'log2FC'
inf_DEG <- inf_DEG %>% rename(log2FC = estimate)
#add column for gene
inf_DEG$gene <- gene
inf_DEG <- inf_DEG %>% select(gene, everything())
#contrast() gives us a t-test p value. We will want to use the chisq p value instead.
#rename t-test p value to differentiate it
inf_DEG <- inf_DEG %>%
	rename(ttest_pvalue = p.value)
#join with the chisq p value
chisq_p <- anova %>% 
	filter(variable == "infected") %>%
	select(gene, `Pr(>Chisq)`)
inf_DEG <- left_join(inf_DEG, chisq_p, by = "gene")

#set up a new dataframe to collect looped results
inf_DEG_all <- inf_DEG

#adjust gene list
genes <- genes[-1]

#Error-tolerant loop:
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
																"infected +",
																"infected/iso_name +",
																"tray +",
																"seq_batch"))
		# Build model
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
			summarise(tray = sum(variance[grepl("^tray", term)]),
								infected = sum(variance[grepl("infectedyes$", term)]),
								seq_batch = sum(variance[grepl("^seq_batch", term)]),
								`infected:iso_name` = sum(variance[grepl("^infectedno:|^infectedyes:", term)], na.rm = T),
								intercept = sum(variance[grepl("Intercept", term)]),
								tot_var = sum(tray, infected, seq_batch, `infected:iso_name`, intercept))
		#calculate percent variance for each term
		var_sums <- var_sums %>%
			mutate(tray = tray/tot_var,
						 infected = infected/tot_var,
						 seq_batch = seq_batch/tot_var,
						 `infected:iso_name` = `infected:iso_name`/tot_var,
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
		emmresult <- emmeans(model, specs = "iso_name")
		# Get emm summary
		emmsummary <- summary(emmresult)
		# Join with existing dataframes
		emm_toadd <- emmsummary$emm %>% as.data.frame()
		colnames(emm_toadd) <- gene
		emm_df <- cbind(emm_df, emm_toadd)
		SE_toadd <- emmsummary$SE %>% as.data.frame()
		colnames(SE_toadd) <- gene
		SE_df <- cbind(SE_df, SE_toadd)
		
		#gather DEG infected data
		#get infected emmeans. note this is natural log scale
		inf_emmresult <- emmeans(model, specs = "infected")
		#compute pairwise contrasts for infected DEGs
		#revpairwise will give infectedyes - infectedno, so the change from mock to infected
		inf_DEG <- contrast(inf_emmresult, method = "revpairwise")
		inf_DEG <- summary(inf_DEG)
		#convert to log2 scale
		inf_DEG$estimate <- inf_DEG$estimate / log(2)
		inf_DEG$SE <- inf_DEG$SE / log(2)
		#change 'estimate' column to 'log2FC'
		inf_DEG <- inf_DEG %>% rename(log2FC = estimate)
		#add column for gene
		inf_DEG$gene <- gene
		inf_DEG <- inf_DEG %>% select(gene, everything())
		#contrast() gives us a t-test p value. We will want to use the chisq p value instead.
		#rename t-test p value to differentiate it
		inf_DEG <- inf_DEG %>%
			rename(ttest_pvalue = p.value)
		#join with the chisq p value
		chisq_p <- anova %>% 
			filter(variable == "infected") %>%
			select(gene, `Pr(>Chisq)`)
		inf_DEG <- left_join(inf_DEG, chisq_p, by = "gene")
		#join result with growing dataframe
		inf_DEG_all <- rbind(inf_DEG_all, inf_DEG)
		
	}, error = function(e) {
		# Handle error: add gene to failed list and print a message
		print(paste("Error encountered for gene:", gene, "Skipping..."))
		failed_genes <<- c(failed_genes, gene) # Append to global list
	})
}

#write out results
dir.create(output_dir)
write.csv(anova_all, paste0(output_dir, "anova.csv"), row.names = F)
write.csv(emm_df, paste0(output_dir, "adjusted_emmeans.csv"), row.names = F)
write.csv(SE_df, paste0(output_dir, "adjusted_SE.csv"), row.names = F)
write.csv(inf_DEG_all, paste0(output_dir, "DEGs_infected.csv"), row.names = F)

# Report the failed genes
print("The following genes caused errors:")
print(failed_genes)
# Convert the list of failed genes to a dataframe
failed_genes_df <- data.frame(gene = failed_genes)
failed_genes_df <- t(failed_genes_df) 
#colnames(failed_genes_df) <- "failed_gene"
# Save the dataframe as a CSV file
write.csv(failed_genes_df, paste0(output_dir, "failed_genes.csv"), row.names = FALSE)
