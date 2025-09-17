config <- config::get()
load(paste0(config$data_results, "Fig_4_comparison_figs.RData"))

amino_differences <- combined %>%
  dplyr::select(amino_valine, amino_isoleucine, amino_leucine, ref_clust)

amino_diff.long <- amino_differences %>%
  pivot_longer(amino_valine:amino_leucine, names_to = 'variable', values_to = 'value')

amino_significance <- amino_differences %>%
  pivot_longer(amino_valine:amino_leucine, names_to = 'variable', values_to = 'value') %>%
  dplyr::group_by(ref_clust, variable) %>%
  dplyr::summarise(y.position = quantile(value, probs = 0.75, type = 7, na.rm = TRUE))

# Extract letters and convert to a dataframe
letters_df_amino_valine <- data.frame(
  ref_clust = names(cld$amino_valine$ref_clust$Letters),
  letter = cld$amino_valine$ref_clust$Letters,
  variable = "amino_valine"
)
letters_df_amino_isoleucine <- data.frame(
  ref_clust = names(cld$amino_isoleucine$ref_clust$Letters),
  letter = cld$amino_isoleucine$ref_clust$Letters,
  variable = "amino_isoleucine"
)
letters_df_amino_leucine <- data.frame(
  ref_clust = names(cld$amino_leucine$ref_clust$Letters),
  letter = cld$amino_leucine$ref_clust$Letters,
  variable = "amino_leucine"
)

letters_df <- dplyr::full_join(letters_df_amino_valine, letters_df_amino_isoleucine) %>%
  dplyr::full_join(letters_df_amino_leucine)

amino_significance <-  dplyr::left_join(amino_significance, 
                                 letters_df, 
                                 by = c('ref_clust', 'variable'))

amino_significance$ref_clust <- as.factor(amino_significance$ref_clust)

g <- ggplot(amino_diff.long, aes(x = ref_clust, y = value)) +
  geom_boxplot(outlier.size = 0.3, fill = "#6baed6") +
  facet_wrap(~factor(variable, 
                     levels=c('amino_valine', 'amino_isoleucine', 'amino_leucine')), 
             scales = 'free',
             labeller = as_labeller(
               c(amino_valine = "Valine",
                 amino_isoleucine = "Isoleucine", 
                 amino_leucine = "Leucine"))) + 
  theme_prism() +
  theme(text = element_text(family = "Source Sans"),
        line = element_line(linewidth = 0.5)) + 
  xlab("Cluster (n)") +
  ylab("Concentration (µM)") +
  geom_text(data = amino_significance, 
            aes(x = ref_clust, y = y.position, label = letter), 
            vjust = -0.5,
            nudge_x = 0.3
  )

ggsave(filename=paste0(config$figures, "Fig_4E_integrated_amino_comparison.svg"), 
       width = 7, 
       height = 4, 
       units = "in")
