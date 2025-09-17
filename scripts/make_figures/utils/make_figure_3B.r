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

symptom_differences <- combined %>%
  dplyr::select(depression, anxiety, diarrhea, constipation, pain_bloating, upper_gi, nausea_and_vomiting, ref_clust) 

symptom_diff.long <- symptom_differences %>%
  pivot_longer(depression:nausea_and_vomiting, names_to = 'variable', values_to = 'value')

g <- ggplot(symptom_diff.long, aes(x = ref_clust, y = value)) +
  geom_boxplot(outlier.size = 0.3,  fill = "#faa61b") +
  facet_wrap(~factor(variable, 
                     levels=c('depression', 'anxiety', 'diarrhea', 'constipation', 
                              'pain_bloating', 'upper_gi', 'nausea_and_vomiting')), 
             scales = 'free',
             labeller = as_labeller(
               c(depression = "Depression", anxiety = "Anxiety", diarrhea = "Diarrhea",
                 constipation = "Constipation", pain_bloating = "Pain & bloating", 
                 upper_gi = "Upper GI", nausea_and_vomiting = "Nausea & vomiting"))) + 
  theme_prism() +
  theme(text = element_text(family = "Source Sans Pro", face = "bold"),
        line = element_line(linewidth = 0.5)) + 
  xlab("Cluster (n)") +
  ylab("Factor score")

ggsave(filename=paste0(config$figures, "Fig_3B_NEMO_symptom_comparison.svg"), 
       plot=g,
       width = 7, 
       height = 7, 
       units = "in")
