load(paste0(config$data_results, "Fig_3A_biology_diagnosis_breakdown.rData"))

p <- ggplot(comb) +
  geom_bar(aes(x = `Cluster (n)`, y = `Count (n)`, fill = Diagnosis), position = "stack", stat = "identity", width = 1, color = "Black") +
  facet_grid(~ Class, switch = "x") +
  theme_minimal(base_size = 20) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(linewidth = 1, color = "#D9D9D9"),
    panel.grid.minor.y = element_line(linewidth = 0.5, color = "#D9D9D9"),
    plot.margin=unit(c(0,0,0,0), 'cm'),
    panel.spacing.x = unit(0, "null"),
    axis.text.x = element_blank()
  ) +
  scale_fill_discrete(type= c(
    '#cecece', # control
    '#a559aa', # FC
    '#EEB1F2', # FC low 
    '#082a54', # FD
    '#B1CEF2', # FD low 
    '#59a88c', # IBS-C
    '#B1F2E3', # IBS-C low 
    '#f0c571', # IBS-D
    '#FFE1AA', # IBS-D low 
    '#e02b35', # IBS-M
    '#FF9FA6', # IBS-M low  
    '#606060'  # Symptomatic control 
  ),
  breaks = c('Control', 'FC', 'FD', 'IBS-C', 'IBS-D', 'IBS-M')
  ) +
  scale_x_discrete(
  )

ggsave(filename=paste0(config$figures, "Fig_3A_biology_cluster_diagnoses.svg"), 
       width = 6, 
       height = 6.5, 
       units = "in")
