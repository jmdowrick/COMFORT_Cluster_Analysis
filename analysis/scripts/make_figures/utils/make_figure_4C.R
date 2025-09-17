config <- config::get()
load(paste0(config$data_results, "Fig_4_comparison_figs.RData"))

org_differences <- combined %>%
  dplyr::select(scfa_aa, scfa_pa, scfa_ba, scfa_va, scfa_suc_a, scfa_x4mva, ref_clust) 

org_diff.long <- org_differences %>%
  pivot_longer(scfa_aa:scfa_x4mva, names_to = 'variable', values_to = 'value')

org_significance <- org_differences %>%
  pivot_longer(scfa_aa:scfa_x4mva, names_to = 'variable', values_to = 'value') %>%
  dplyr::group_by(ref_clust, variable) %>%
  dplyr::summarise(y.position = quantile(value, probs = 0.75, type = 7, na.rm = TRUE))

# Extract letters and convert to a dataframe
letters_df_scfa_aa <- data.frame(
  ref_clust = names(cld$scfa_aa$Letters),
  letter = cld$scfa_aa$Letters,
  variable = "scfa_aa"
)
letters_df_scfa_pa <- data.frame(
  ref_clust = names(cld$scfa_pa$Letters),
  letter = cld$scfa_pa$Letters,
  variable = "scfa_pa"
)
letters_df_scfa_ba <- data.frame(
  ref_clust = names(cld$scfa_ba$ref_clust$Letters),
  letter = cld$scfa_ba$ref_clust$Letters,
  variable = "scfa_ba"
)
letters_df_scfa_va <- data.frame(
  ref_clust = names(cld$scfa_va$ref_clust$Letters),
  letter = cld$scfa_va$ref_clust$Letters,
  variable = "scfa_va"
)
letters_df_scfa_suc_a <- data.frame(
  ref_clust = names(cld$scfa_suc_a$Letters),
  letter = cld$scfa_suc_a$Letters,
  variable = "scfa_suc_a"
)
letters_df_scfa_x4mva <- data.frame(
  ref_clust = names(cld$scfa_x4mva$Letters),
  letter = cld$scfa_x4mva$Letters,
  variable = "scfa_x4mva"
)

letters_df <- dplyr::full_join(letters_df_scfa_aa, letters_df_scfa_pa) %>%
  dplyr::full_join(letters_df_scfa_ba) %>%
  dplyr::full_join(letters_df_scfa_va) %>%
  dplyr::full_join(letters_df_scfa_suc_a) %>%
  dplyr::full_join(letters_df_scfa_x4mva) 

org_significance <-  dplyr::left_join(org_significance, 
                               letters_df, 
                               by = c('ref_clust', 'variable'))

org_significance$ref_clust <- as.factor(org_significance$ref_clust)

g <- ggplot(org_diff.long, aes(x = ref_clust, y = value)) +
  geom_boxplot(outlier.size = 0.3, fill = "#6baed6") +
  scale_y_continuous(trans='log10') +
  facet_wrap(~factor(variable, 
                     levels=c('scfa_aa', 'scfa_ba', 'scfa_pa', 'scfa_va', 'scfa_suc_a', 'scfa_x4mva')), 
             scales = 'free',
             labeller = as_labeller(
               c(scfa_aa = "Acetic acid", 
                 scfa_pa = "Propionic acid", 
                 scfa_ba = "Butyric acid",
                 scfa_va = "Valeric acid", 
                 scfa_suc_a = "Succinic acid",
                 scfa_x4mva = "Isocaproic acid")
             )
  ) + 
  theme_prism() +
  theme(text = element_text(family = "Source Sans"),
        line = element_line(linewidth = 0.5)) + 
  xlab("Cluster (n)") +
  ylab("Log10(Concentration) (µmol/g)") +
  geom_text(data = org_significance, 
            aes(x = ref_clust, y = y.position, label = letter), 
            vjust = -0.5,
            nudge_x = 0.3
  )

ggsave(filename=paste0(config$figures, "Fig_4C_integrated_SCFA_comparison.svg"), 
       width = 7, 
       height = 7, 
       units = "in")
