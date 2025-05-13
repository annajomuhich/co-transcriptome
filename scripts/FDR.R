### FDR correction for gene by gene p values
### May 2025 AJM

anova <- read.csv("./remove_7isos/Vu_hostmodel/anova.csv")
library(tidyverse)

# Split data by variable type
anova_split <- split(anova, anova$variable)
# Apply FDR correction to each variable group
anova_fdr <- lapply(anova_split, function(x) {
  p_values <- x$`Pr..Chisq.`
  x$p_adj <- p.adjust(p_values, method = "BH")
  return(x)})
# Recombine into single dataframe
anova_corrected <- do.call(rbind, anova_fdr)
# Reset row names
rownames(anova_corrected) <- NULL

# Write results
write.csv(anova_corrected, "./remove_7isos/Vu_hostmodel/anova_FDR.csv", row.names = FALSE)
