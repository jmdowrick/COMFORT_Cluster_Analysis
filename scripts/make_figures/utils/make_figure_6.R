load(paste0(config$data_results, "Fig_6_FEA_results.RData"))

# Figure 6A
enrichres_up_11vs6@result$significant <- NA
enrichres_up_11vs6@result$significant[enrichres_up_11vs6@result$`p.adjust`<0.1] <- '*'
enrichres_up_11vs6@result$GeneRatioCONT <- sapply(parse(text = enrichres_up_11vs6@result$GeneRatio), eval)

p <- enrichplot::dotplot(enrichres_up_11vs6) + 
  ggtitle("A. Upregulated") +
  scale_x_continuous(limits = c(0.047, 0.285), breaks = seq(0.1, 0.2, by = 0.1), minor_breaks = seq(0.05, 0.25, by = 0.1)) +
  geom_point(stroke = 0.9) +
  geom_text(aes(label = enrichres_up_11vs6@result$significant, x = enrichres_up_11vs6@result$GeneRatioCONT, y = enrichres_up_11vs6$Description), size = 7, colour= "black", nudge_y = 0.3) +
  theme(text = element_text(size = 14, family = "Source Sans Pro", face="bold"),
        axis.text.x = element_text(size = 14, family = "Source Sans Pro", face="plain"), 
        axis.text.y = element_text(size = 12, family = "Source Sans Pro", face="plain"), 
        axis.title = element_text(size = 14, family = "Source Sans Pro", face="plain"), 
        legend.title = element_text(size=14, family = "Source Sans Pro", face="bold"),
        legend.text = element_text(size = 10, family = "Source Sans Pro", face="plain"), 
        panel.border = element_rect(colour = "#323232", fill = NA, linewidth=1))

ggsave(paste0(config$figures, "Fig_6A_up_metagenomic_11vs6.svg"), 
       width = 6,
       height = 5, 
       units = "in")

# Figure 6B
enrichres_down_11vs6@result$significant <- NA
enrichres_down_11vs6@result$significant[enrichres_down_11vs6@result$`p.adjust`<0.1] <- '*'
enrichres_down_11vs6@result$GeneRatioCONT <- sapply(parse(text = enrichres_down_11vs6@result$GeneRatio), eval)

p <- enrichplot::dotplot(enrichres_down_11vs6) + 
  ggtitle("B. Downregulated") +
  scale_x_continuous(limits = c(0.04, 0.163), breaks = seq(0.08, 0.16, by = 0.04), minor_breaks = seq(0.04, 0.1, by = 0.04)) +
  geom_point(stroke = 0.9) +
  geom_text(aes(label = enrichres_down_11vs6@result$significant[1:10], x = enrichres_down_11vs6@result$GeneRatioCONT[1:10], y = enrichres_down_11vs6$Description[1:10]), size = 7, colour= "black", nudge_y = 0.35) +
  theme(text = element_text(size = 14, family = "Source Sans Pro", face="bold"),
        axis.text.x = element_text(size = 14, family = "Source Sans Pro", face="plain"), 
        axis.text.y = element_text(size = 12, family = "Source Sans Pro", face="plain"), 
        axis.title = element_text(size = 14, family = "Source Sans Pro", face="plain"), 
        legend.title = element_text(size=14, family = "Source Sans Pro", face="bold"),
        legend.text = element_text(size = 10, family = "Source Sans Pro", face="plain"), 
        panel.border = element_rect(colour = "#323232", fill = NA, linewidth=1))

ggsave(paste0(config$figures, "Fig_6B_down_metagenomic_11vs6.svg"), 
       width = 6,
       height = 5, 
       units = "in")

# Figure 6C
res_fmetab_up_11vs6$significant <- NA
res_fmetab_up_11vs6$significant[res_fmetab_up_11vs6$FDR<0.1] <- '*'

p <- ggplot(res_fmetab_up_11vs6, aes(x = minuslog10p, y = metaboliteSet)) +
  ggtitle("C. Overabundant (fecal)") +
  xlab("-log10(P)") +
  ylab("") +
  geom_point(aes(size = EnrichRatio, colour = FDR)) +
  geom_point(aes(size = EnrichRatio), fill = NA, shape = 1, stroke = 0.9, colour = "black") + 
  geom_text(aes(label = significant, x = minuslog10p, y = metaboliteSet), size = 7, colour= "black", nudge_y = 0.3) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) + 
  scale_x_continuous(limits = c(3.3, 4.4), breaks = seq(3.5, 4.5, by = 0.5), minor_breaks = seq(3.25, 4.25, by = 0.5)) +
  scale_color_gradient(low = "#357fb9", high = "#df6563", name = "q value", trans = 'reverse') +
  scale_size_continuous(name = "Enrich Ratio") +
  guides(color = guide_colorbar(order = 1), size = guide_legend(override.aes = list(fill = NA, shape = 1, colour = "black"), order = 2)) +
  theme_minimal() +
  theme(text = element_text(size = 14, family = "Source Sans Pro", face="bold", colour="#010101"),
        axis.text.x = element_text(size = 14, family = "Source Sans Pro", face="plain", colour="#010101"), 
        axis.text.y = element_text(size = 11, family = "Source Sans Pro", face="plain", colour="#010101"), 
        axis.title = element_text(size = 14, family = "Source Sans Pro", face="plain"), 
        legend.title = element_text(size=14, family = "Source Sans Pro", face="bold"),
        legend.text = element_text(size = 10, family = "Source Sans Pro", face="plain"), 
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), 
        axis.ticks = element_line())

ggsave(paste0(config$figures, "Fig_6C_up_fmetab_11vs6.svg"), 
       width = 7,
       height = 5, 
       units = "in")
