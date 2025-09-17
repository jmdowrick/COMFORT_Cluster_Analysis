symp_stability <- read_csv(paste0(config$data_results, "Fig_S1_symptom_stability.csv")) %>%
  dplyr::select(c(cluster, item_consensus)) 

symp_stability$significant <- ifelse(symp_stability$cluster == 2, "","*")
symp_stability$colour <- ifelse(symp_stability$cluster == 2, "#e8e8e8","#078281")

asteriskPos <- sapply(symp_stability$cluster, function(cl) {max(symp_stability$item_consensus[symp_stability$cluster==cl]+0.05)})
symp_stability$asteriskPos <- asteriskPos

p <- symp_stability %>% 
  ggplot(aes(x = factor(cluster), y = item_consensus, fill = colour)) +
  xlab("Cluster (n)") +
  ylab("Item consensus") +
  stat_halfeye(
    adjust = 0.5,
    justification = -0.2,
    .width = 0,
    point_colour = NA,
    show.legend = FALSE
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.color = "black",
    alpha = 0.7,
    show.legend = FALSE,
    size = 1
  ) + 
  geom_text(aes(label = significant, x = factor(cluster), y = asteriskPos), size = 7, colour= "black") +
  scale_fill_identity()+
  geom_hline(yintercept=0.5, linetype="dashed", color = "black") + 
  scale_y_continuous(n.breaks=6, limits = c(0,1)) + 
  theme_bw(base_size = 20) + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(text = element_text(size = 20, family = "Source Sans Pro", face = "bold"),
        axis.text.x = element_text(size = 20, family = "Source Sans Pro", face = "plain"), 
        axis.text.y = element_text(size = 20, family = "Source Sans Pro", face = "plain"), 
        axis.title = element_text(size = 20, "Source Sans Pro", face = "bold"))

ggsave(filename=paste0(config$sup_figures, "Fig_S1_symptom_stability.svg"), 
       width = 8, 
       height = 6, 
       units = "in")
