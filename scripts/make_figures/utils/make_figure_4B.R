config <- config::get()
load(paste0(config$data_results, "Fig_4_comparison_figs.RData"))

symptom_differences <- combined %>%
  dplyr::select(depression, anxiety, diarrhea, constipation, pain_bloating, upper_gi, ref_clust) 

symptom_diff.long <- symptom_differences %>%
  pivot_longer(depression:upper_gi, names_to = 'variable', values_to = 'value')

symptom_significance <- symptom_differences %>%
  pivot_longer(depression:upper_gi, names_to = 'variable', values_to = 'value') %>%
  dplyr::group_by(ref_clust, variable) %>%
  dplyr::summarise(y.position = quantile(value, probs = 0.75, type = 7, na.rm = TRUE))

# Extract letters and convert to a dataframe
letters_df_depression <- data.frame(
  ref_clust = names(cld$depression$ref_clust$Letters),
  letter = cld$depression$ref_clust$Letters,
  variable = "depression"
)
letters_df_anxiety <- data.frame(
  ref_clust = names(cld$anxiety$ref_clust$Letters),
  letter = cld$anxiety$ref_clust$Letters,
  variable = "anxiety"
)
letters_df_diarrhea <- data.frame(
  ref_clust = names(cld$diarrhea$Letters),
  letter = cld$diarrhea$Letters,
  variable = "diarrhea"
)
letters_df_constipation <- data.frame(
  ref_clust = names(cld$constipation$ref_clust$Letters),
  letter = cld$constipation$ref_clust$Letters,
  variable = "constipation"
)
letters_df_pain_bloating <- data.frame(
  ref_clust = names(cld$pain_bloating$ref_clust$Letters),
  letter = cld$pain_bloating$ref_clust$Letters,
  variable = "pain_bloating"
)
letters_df_upper_gi <- data.frame(
  ref_clust = names(cld$upper_gi$ref_clust$Letters),
  letter = cld$upper_gi$ref_clust$Letters,
  variable = "upper_gi"
)

letters_df <- dplyr::full_join(letters_df_depression, letters_df_anxiety) %>%
  dplyr::full_join(letters_df_diarrhea) %>%
  dplyr::full_join(letters_df_constipation) %>%
  dplyr::full_join(letters_df_pain_bloating) %>%
  dplyr::full_join(letters_df_upper_gi)

symptom_significance <-  dplyr::left_join(symptom_significance, 
                                   letters_df, 
                                   by = c('ref_clust', 'variable'))

symptom_significance$ref_clust <- as.factor(symptom_significance$ref_clust)

g <- ggplot(symptom_diff.long, aes(x = ref_clust, y = value)) +
  geom_boxplot(outlier.size = 0.3, fill = "#faa61b") +
  facet_wrap(~factor(variable, 
                     levels=c('depression', 'anxiety', 'diarrhea', 'constipation', 
                              'pain_bloating', 'upper_gi', 'nausea_and_vomiting')), 
             scales = 'free',
             labeller = as_labeller(
               c(depression = "Depression", anxiety = "Anxiety", diarrhea = "Diarrhea",
                 constipation = "Constipation", pain_bloating = "Pain & bloating", 
                 upper_gi = "Upper GI", nausea_and_vomiting = "Nausea & vomiting"))) + 
  theme_prism() +
  theme(text = element_text(family = "Source Sans"),
        line = element_line(linewidth = 0.5)) + 
  xlab("Cluster (n)") +
  ylab("Factor score") +
  geom_text(data = symptom_significance, 
            aes(x = ref_clust, y = y.position, label = letter), 
            vjust = -0.5,
            nudge_x = 0.3
  )

ggsave(filename=paste0(config$figures, "Fig_4B_integrated_symptom_comparison.svg"), 
       width = 7, 
       height = 7, 
       units = "in")
