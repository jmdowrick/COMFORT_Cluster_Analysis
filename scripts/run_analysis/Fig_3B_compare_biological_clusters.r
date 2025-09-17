# Biological characterisation for NEMO clustering
config <- config::get()

# Import biological reference cluster labels
ref_clust <- read.csv(paste0(config$data_cluster_assignments, "biological_clusters.csv"), header=TRUE, row.names = 1) 

# Import factor scores
factor_scores <- read.table(paste0(config$data_results, "comfort_factor_scores.dat"), header=FALSE) 

factor_scores <- factor_scores %>%
  dplyr::select('V43', 'V45', 'V47', 'V49', 'V51', 'V53', 'V55') %>%
  dplyr::rename(., c("Depression" = "V43", "Anxiety" = "V45", 
                     "Diarrhea" = "V47", "Constipation" = "V49", 
                     "Pain/bloating" = "V51", "Upper GI" = "V53", 
                     "Nausea and vomiting" = "V55"))

ids = readr::read_csv(paste0(config$data_mplus, "comfort_PROs.csv"))
rownames(factor_scores) <- ids$id

# Combine all symptoms scores with cluster labels
factor_scores <- factor_scores %>% 
  janitor::clean_names() 

combined <- merge(factor_scores,
                       as.data.frame(ref_clust),
                       by='row.names')

row.names(combined) <- combined$Row.names
combined$Row.names <- NULL
combined$ref_clust <- as.factor(combined$ref_clust)

# Perform comparison between clusters
results <- list()
measures <- colnames(combined)
measures <- measures[measures != "ref_clust"]

for (i in 1:length(measures)) {
  measure_name <- measures[[i]]
  formula <- as.formula(paste(measure_name, "~ ref_clust"))
  fit <- lm(formula, data = combined, na.action=na.omit)
  results[[measure_name]] <- summary(aov(fit))
}

p_values <- sapply(results, function(x) x[[1]]$'Pr(>F)'[1])
p_adjusted <- p.adjust(p_values, method = "bonferroni")
significant_measures <- which(p_adjusted < 0.05)

p_adjusted_df <- as.data.frame(p_adjusted)
p_adjusted_df <- as.data.frame(cbind(p_adjusted_df, row.names(p_adjusted_df)))
write.csv(p_adjusted_df, file = paste0(config$data_results, "Fig_3B_bio_aov_p_values.csv"))
