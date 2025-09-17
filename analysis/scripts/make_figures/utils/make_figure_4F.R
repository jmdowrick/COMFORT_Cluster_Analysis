config <- config::get()
load(paste0(config$data_results, "Fig_4_comparison_figs.RData"))

# Metagenomics (taxonomic abundance)
microbiome_differences <- combined %>%
  dplyr::select(firm_bac_ratio, average_richness, shannon, ref_clust) 

microbiome_diff.long <- microbiome_differences %>%
  pivot_longer(firm_bac_ratio:shannon, names_to = 'variable', values_to = 'value')

microbiome_diff.long$value[microbiome_diff.long$variable == "firm_bac_ratio"] <- log10(microbiome_diff.long$value[microbiome_diff.long$variable == "firm_bac_ratio"])
microbiome_significance <- microbiome_diff.long %>%
  dplyr::group_by(ref_clust, variable) %>%
  dplyr::summarise(y.position = quantile(value, probs = 0.75, type = 7, na.rm = TRUE))

# Extract letters and convert to a dataframe
letters_df_firm_bac_ratio <- data.frame(
  ref_clust = names(cld$firm_bac_ratio$Letters),
  letter = cld$firm_bac_ratio$Letters,
  variable = "firm_bac_ratio"
)
letters_df_shannon <- data.frame(
  ref_clust = names(cld$shannon$ref_clust$Letters),
  letter = cld$shannon$ref_clust$Letters,
  variable = "shannon"
)
letters_df_average_richness <- data.frame(
  ref_clust = names(cld$average_richness$ref_clust$Letters),
  letter = cld$average_richness$ref_clust$Letters,
  variable = "average_richness"
)

letters_df <- dplyr::full_join(letters_df_firm_bac_ratio, letters_df_shannon)%>%
  dplyr::full_join(letters_df_average_richness)

microbiome_significance <-  dplyr::left_join(microbiome_significance, 
                                      letters_df, 
                                      by = c('ref_clust', 'variable'))

microbiome_significance$ref_clust <- as.factor(microbiome_significance$ref_clust)

g <- ggplot(microbiome_diff.long, aes(x = ref_clust, y = value)) +
  geom_boxplot(outlier.size = 0.3, fill = "#9f9ac7") +
  facet_wrap(~factor(variable, 
                     levels=c('firm_bac_ratio', 'shannon', 'average_richness')), 
             scales = 'free',
             labeller = as_labeller(
               c(firm_bac_ratio = "Log(FB Ratio)", shannon = "Shannon",
                 average_richness = "Richness"))) + 
  theme_prism() +
  theme(text = element_text(family = "Source Sans"),
        line = element_line(linewidth = 0.5)) + 
  xlab("Cluster (n)") +
  ylab("Value") +
  geom_text(data = microbiome_significance, 
            aes(x = ref_clust, y = y.position, label = letter), 
            vjust = -0.5,
            nudge_x = 0.3
  )

ggsave(filename=paste0(config$figures, "Fig_4F_integrated_metagenomics_comparison.svg"), 
       width = 7, 
       height = 4, 
       units = "in")
