config <- config::get()
load(paste0(config$data_results, "Fig_4_comparison_figs.RData"))

# Bile acids ----
bile_differences <- combined %>%
  dplyr::select(bile_gdca, bile_hdca, bile_glca, bile_taurine, bile_cdca, 
                bile_ca, bile_tlca, bile_tcdca, bile_tdca, ref_clust) 

bile_differences$ref_clust <- as.factor(bile_differences$ref_clust)

bile_diff.long <- bile_differences %>%
  pivot_longer(bile_gdca:bile_tdca, names_to = 'variable', values_to = 'value') 

bile_significance <- bile_differences %>%
  pivot_longer(bile_gdca:bile_tdca, names_to = 'variable', values_to = 'value') %>%
  dplyr::group_by(ref_clust, variable) %>%
  dplyr::summarise(y.position = quantile(value, probs = 0.75, type = 7, na.rm = TRUE))

# Extract letters and convert to a dataframe

letters_df_gdca <- data.frame(
  ref_clust = names(cld$bile_gdca$Letters),
  letter = cld$bile_gdca$Letters,
  variable = "bile_gdca"
)
letters_df_hdca <- data.frame(
  ref_clust = names(cld$bile_hdca$Letters),
  letter = cld$bile_hdca$Letters,
  variable = "bile_hdca"
)
letters_df_glca <- data.frame(
  ref_clust = names(cld$bile_glca$Letters),
  letter = cld$bile_glca$Letters,
  variable = "bile_glca"
)
letters_df_taurine <- data.frame(
  ref_clust = names(cld$bile_taurine$Letters),
  letter = cld$bile_taurine$Letters,
  variable = "bile_taurine"
)
letters_df_cdca <- data.frame(
  ref_clust = names(cld$bile_cdca$Letters),
  letter = cld$bile_cdca$Letters,
  variable = "bile_cdca"
)
letters_df_ca <- data.frame(
  ref_clust = names(cld$bile_ca$Letters),
  letter = cld$bile_ca$Letters,
  variable = "bile_ca"
)
letters_df_tlca <- data.frame(
  ref_clust = names(cld$bile_tlca$Letters),
  letter = cld$bile_tlca$Letters,
  variable = "bile_tlca"
)
letters_df_tcdca <- data.frame(
  ref_clust = names(cld$bile_tcdca$Letters),
  letter = cld$bile_tcdca$Letters,
  variable = "bile_tcdca"
)
letters_df_tdca <- data.frame(
  ref_clust = names(cld$bile_tdca$Letters),
  letter = cld$bile_tdca$Letters,
  variable = "bile_tdca"
)

letters_df <- dplyr::full_join(letters_df_gdca, letters_df_hdca) %>%
  dplyr::full_join(letters_df_glca) %>%
  dplyr::full_join(letters_df_taurine) %>%
  dplyr::full_join(letters_df_cdca) %>%
  dplyr::full_join(letters_df_ca) %>%
  dplyr::full_join(letters_df_tlca) %>%
  dplyr::full_join(letters_df_tcdca) %>%
  dplyr::full_join(letters_df_tdca) 

letters_df$ref_clust <- as.factor(letters_df$ref_clust)

bile_significance <-  dplyr::left_join(bile_significance, 
                                letters_df, 
                                by = c('ref_clust', 'variable'))

# gdca - glycoheoxycholic acid 
# hdca - hyodeoxycholic acid
# ila - isolithocholic acid
# lca - lithocholic acid
# cdca - chenodeoxycholic acid
# tcdca - taurochenodeoxycholic acid
# tdca - taurodeoxycholic acid

g <- ggplot(bile_diff.long, aes(x = ref_clust, y = value)) +
  geom_boxplot(outlier.size = 0.3, fill = "#69badf") +
  scale_y_continuous(trans='log10') + 
  facet_wrap(~variable, scales = 'free',
             labeller = as_labeller(
               c(bile_gdca = "GDCA", bile_hdca = "HDCA", bile_glca = "GLCA",
                 bile_taurine = "Taurine", bile_cdca = "CDCA", bile_ca = "CA",
                 bile_tlca = "TLCA", bile_tcdca = "TCDCA", bile_tdca = "TDCA")),
             ncol = 3) + 
  ylab("Log_{10}(Concentration) (µg/g)") +
  theme_prism() +
  theme(text = element_text(family = "Source Sans"),
        line = element_line(linewidth = 0.5)) + 
  xlab("Cluster (n)") + 
  geom_text(data = bile_significance, 
            aes(x = ref_clust, y = y.position, label = letter), 
            vjust = -0.5,
            nudge_x = 0.3
  )

ggsave(filename=paste0(config$figures, "Fig_4D_integrated_bile_comparison.svg"), 
       width = 7, 
       height = 7, 
       units = "in")
