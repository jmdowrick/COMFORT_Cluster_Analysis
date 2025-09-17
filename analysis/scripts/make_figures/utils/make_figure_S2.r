
bio_stability <- read_csv(paste0(config$data_results, "Fig_S2_biology_stability.csv")) %>%
  dplyr::select(c(cluster, item_consensus)) 

bio_stability$significant <- "*"
bio_stability$colour <- "#078281"

asteriskPos <- sapply(bio_stability$cluster, function(cl) {max(bio_stability$item_consensus[bio_stability$cluster==cl]+0.05)})

asteriskPos[asteriskPos >= 1] <- 1
bio_stability$asteriskPos <- asteriskPos

p <- bio_stability %>% 
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
  coord_cartesian(clip = "off") +
  geom_text(aes(label = significant, x = factor(cluster), y = asteriskPos), size = 7, colour= "black") +
  scale_fill_identity()+
  geom_hline(yintercept=0.5, linetype="dashed", color = "black") + 
  scale_y_continuous(n.breaks=6, limits = c(0,1)) + 
  theme_bw(base_size = 20) + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(text = element_text(size = 20, family = "Source Sans Pro", face = "bold"),
        axis.text.x = element_text(size = 20, family = "Source Sans Pro", face = "plain"), 
        axis.text.y = element_text(size = 20, family = "Source Sans Pro", face = "plain"), 
        axis.title = element_text(size = 20, family = "Source Sans Pro", face = "bold"))

ggsave(filename=paste0(config$sup_figures, "Fig_S2_biology_stability.svg"), 
       width = 8, 
       height = 6, 
       units = "in")
