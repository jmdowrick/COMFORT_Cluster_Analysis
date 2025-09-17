library(tidyr)
config <- config::get()

# Load data ----
comfort_amino <- read.csv(paste0(config$data_processed, "comfort_amino_acids_cluster.csv"), row.names = 1)
comfort_bile <- read.csv(paste0(config$data_processed, "comfort_bile_acids_cluster.csv"), row.names = 1)
comfort_organic <- read.csv(paste0(config$data_processed, "comfort_organic_acids_cluster.csv"), row.names = 1)

names(comfort_amino) <- substring(names(comfort_amino), 2)
names(comfort_bile) <- substring(names(comfort_bile), 2)
names(comfort_organic) <- substring(names(comfort_organic), 2)

org_acids <- as.data.frame(t(comfort_organic))
colnames(org_acids) <- paste0('scfa_', colnames(org_acids))
org_acids <- org_acids %>% tibble::rownames_to_column(var = "rowname")
rm(comfort_organic)

amino_acids <- as.data.frame(t(comfort_amino))
colnames(amino_acids) <- paste0('amino_', colnames(amino_acids))
amino_acids <- amino_acids %>% tibble::rownames_to_column(var = "rowname")
rm(comfort_amino)

bile_acids <- as.data.frame(t(comfort_bile))
colnames(bile_acids) <- paste0('bile_', colnames(bile_acids))
bile_acids <- bile_acids %>% tibble::rownames_to_column(var = "rowname")
rm(comfort_bile)

metagenomics <- read.csv(paste0(config$data_processed, "comfort_metagenomics_diversityAndFBratio.csv"), row.names = 1)
metagenomics <- metagenomics %>% tibble::rownames_to_column(var = "rowname")

# Load symptom factor scores ----
factor_scores <- read.table(paste0(config$data_processed, "comfort_factor_scores.dat"), 
                            header=FALSE) 

factor_scores <- factor_scores %>%
  dplyr::select('V43', 'V45', 'V47', 'V49', 'V51', 'V53', 'V55') %>%
  dplyr::rename(., c("Depression" = "V43", "Anxiety" = "V45", 
                     "Diarrhea" = "V47", "Constipation" = "V49", 
                     "Pain/bloating" = "V51", "Upper GI" = "V53", 
                     "Nausea and vomiting" = "V55"))
ids = readr::read_csv(paste0(config$data_mplus, "comfort_PROs.csv"))
rownames(factor_scores) <- ids$id
factor_scores <- factor_scores %>% 
  janitor::clean_names() %>% 
  tibble::rownames_to_column(var = "rowname")

# Import diet, age, and exercise data ----
selected_diet <- c("id", "ENERC_sum", "PROT_sum", "FIBTG_sum", "FASAT_sum", "GLUS_sum", "SUGAR_sum")

comfort_diet <- readxl::read_excel(paste0(config$data_raw, "Diet diaries food aggregate ALL DAYS.xlsx"), sheet = "ALL DAYS AGGREGATED") %>%
  dplyr::select(all_of(selected_diet)) %>%
  as.data.frame() %>%
  dplyr::rename(rowname = id)

comfort_diet$rowname <- as.character(comfort_diet$rowname)

comfort_demographic <- readxl::read_excel(paste0(config$data_raw, config$file_PRO), sheet = "ALL_TOTAL") %>%
  dplyr::select("Participant #", AGE, METsum) %>%
  as.data.frame() %>%
  dplyr::rename(rowname = "Participant #")
comfort_demographic$rowname <- as.character(comfort_demographic$rowname)

# Import labels ----
merged_labels <- read.csv(paste0(config$data_cluster_assignments, "integrated_clusters.csv"), row.names = 1) %>%
  tibble::rownames_to_column(var = "rowname")

merged_labels$rowname <- as.character(merged_labels$rowname)
merged_labels$ref_clust <- as.factor(merged_labels$ref_clust)

# Combine all measurements with cluster labels----
combined <- merged_labels %>%
  full_join(bile_acids, by = "rowname") %>%
  full_join(org_acids, by = "rowname") %>%
  full_join(amino_acids, by = "rowname") %>%
  full_join(metagenomics, by = "rowname") %>%
  full_join(factor_scores, by = "rowname") %>%
  full_join(comfort_demographic, by = "rowname") %>%
  full_join(comfort_diet, by = "rowname") %>%
  drop_na(ref_clust)

combined <- combined %>% tibble::column_to_rownames(var = "rowname")

# Drop unstable clusters ----
combined <- combined %>%
  dplyr::filter(ref_clust %in% c(1, 2, 5, 6, 7, 8, 9, 11)) %>%
  droplevels()

combined$ref_clust <- as.factor(as.numeric(as.character(combined$ref_clust)))

# Perform normality test ----
listNames <- colnames(combined) 
listNames <- listNames[2:length(listNames)]

# Create an empty list to store the results
normality_results <- list()

# Perform Shapiro-Wilk test for each variable within each group
for (cluster in unique(combined$ref_clust)) {
  group_data <- combined[combined$ref_clust == cluster, ]
  
  # Apply Shapiro-Wilk test to each variable within this group
  group_results <- sapply(group_data[, listNames], function(x) shapiro.test(na.omit(x))$p.value)
  
  # Store the results in the list with the group name
  normality_results[[cluster]] <- group_results
}

# Combine all results into a data frame for easier manipulation
normality_results_df <- do.call(rbind, lapply(names(normality_results), function(group) {
  data.frame(Group = group, Variable = names(normality_results[[group]]), 
             P_Value = normality_results[[group]])
}))

# Apply multiple testing correction 
normality_results_df$P_Adjusted <- p.adjust(normality_results_df$P_Value, method = "bonferroni")
non_normal <- normality_results_df[normality_results_df$P_Adjusted < 0.05, ]

# Get list of unique names
non_normal_names <- unique(non_normal$Variable)

cld <- list()

# ---- 
# Statistical comparisons for non-normal variables
results <- list()
measures <- colnames(combined)
measures <- measures[measures != "ref_clust"]
measures <- measures[measures %in% non_normal_names] # select non-normal variables

results <- list()
p_values <- c()

# Perform Kruskal-Wallis test for each measure
for (i in 1:length(measures)) {
  measure_name <- measures[[i]]
  formula <- as.formula(paste(measure_name, "~ ref_clust"))
  kruskal_result <- kruskal.test(formula, data = combined)
  results[[measure_name]] <- kruskal_result
  p_values[measure_name] <- kruskal_result$p.value
}

# Adjust p-values for multiple comparisons
p_adjusted <- p.adjust(p_values, method = "bonferroni", n = length(listNames))
significant_measures_nonnormal <- names(which(p_adjusted < 0.05))

# Convert p-values to a data frame and save
p_adjusted_df <- as.data.frame(p_adjusted)
p_adjusted_df <- cbind(measure = rownames(p_adjusted_df), p_adjusted_df)
write.csv(p_adjusted_df, paste0(config$data_results, "Fig_4_kruskal_adj_p_values_list.csv"), row.names = FALSE)

# Perform Dunn's test for pairwise comparisons for significant measures
dunn_results <- list()
p_values_list <- list()

# Loop through each non-normal variable and perform Dunn's test
for (measure in significant_measures_nonnormal) {
  formula <- as.formula(paste(measure, "~ ref_clust")) 
  dunn_test_result <- FSA::dunnTest(formula, data = combined, method = "bh") 
  dunn_results[[measure]] <- dunn_test_result$res$P.adj
  names(dunn_results[[measure]]) <- gsub(" ", "", dunn_test_result$res$Comparison)
  cld[[measure]] <- multcompView::multcompLetters(dunn_results[[measure]])
}

# Combine the p-values into a data frame
dunn_p_values_df <- do.call(cbind, lapply(dunn_results, as.data.frame))

# Transpose the data frame to have measures as rows
dunn_p_values_df <- t(dunn_p_values_df)
rownames(dunn_p_values_df) <- names(dunn_results)

# Add variable names as a column
dunn_p_values_df <- as.data.frame(cbind(dunn_p_values_df, rownames(dunn_p_values_df)))
colnames(dunn_p_values_df)[ncol(dunn_p_values_df)] <- "sig_meas_names"
colnames(dunn_p_values_df) <- gsub(" ", "", colnames(dunn_p_values_df))

# ----
# Statistical comparisons for normal variables
results <- list()
measures <- colnames(combined)
measures <- measures[measures != "ref_clust"]
measures <- measures[!(measures %in% non_normal_names)] # filter out non-normal variables

for (i in 1:length(measures)) {
  measure_name <- measures[[i]]
  formula <- as.formula(paste(measure_name, "~ ref_clust"))
  fit <- lm(formula, data = combined, na.action=na.omit)
  results[[measure_name]] <- summary(aov(fit))
}

p_values <- sapply(results, function(x) x[[1]]$'Pr(>F)'[1])
p_adjusted <- p.adjust(p_values, method = "bonferroni", n = length(listNames))
significant_measures <- which(p_adjusted < 0.05)

p_adjusted_df <- as.data.frame(p_adjusted)
p_adjusted_df <- as.data.frame(cbind(p_adjusted_df, row.names(p_adjusted_df)))

# Find which groups differ from list of significant measures ----
sig_meas_names_normal <- names(significant_measures) 
results_tukeyhsd <- list()
p_values_list <-list()
for (measure_name in sig_meas_names_normal) {
  formula <- as.formula(paste(measure_name, "~ ref_clust"))
  fit <- lm(formula, data = combined, na.action=na.omit)
  results_tukeyhsd <- TukeyHSD(aov(fit))
  p_values_list[[measure_name]] <- results_tukeyhsd[[1]][, "p adj"]
  cld[[measure_name]] <- multcompView::multcompLetters4(aov(fit), TukeyHSD(aov(fit)))
}
tukey_p_values_df <- do.call(cbind, lapply(p_values_list, as.data.frame))

# transpose the data frame to have measures as rows
tukey_p_values_df <- t(tukey_p_values_df)
rownames(tukey_p_values_df) <- sig_meas_names_normal

tukey_p_values_df <- as.data.frame(cbind(tukey_p_values_df, sig_meas_names_normal))

# Combine non-normal and normal p-values for plotting
standardize_comparisons <- function(comparison) {
  parts <- strsplit(comparison, "-")[[1]]
  # Order the parts consistently (smaller number first)
  standardized <- paste(sort(parts), collapse = "-")
  return(standardized)
}

colnames(tukey_p_values_df) <- sapply(colnames(tukey_p_values_df), standardize_comparisons)
colnames(dunn_p_values_df) <- sapply(colnames(dunn_p_values_df), standardize_comparisons)

p_values_df <- full_join(tukey_p_values_df, dunn_p_values_df)

write.csv(p_values_df, paste0(config$data_results, "Fig_4_significant_merged_p_values.csv"), row.names = FALSE)

save(cld, combined, file=paste0(config$data_results, "Fig_4_comparison_figs.RData"))
