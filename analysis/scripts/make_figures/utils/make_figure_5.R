load(paste0(config$data_results, "Fig_5_FEA_results.RData"))

# Figure 5A
enrichres_up_1vs6@result$significant <- NA
enrichres_up_1vs6@result$significant[enrichres_up_1vs6@result$`p.adjust`<0.1] <- '*'
enrichres_up_1vs6@result$GeneRatioCONT <- sapply(parse(text = enrichres_up_1vs6@result$GeneRatio), eval)

p <- enrichplot::dotplot(enrichres_up_1vs6) + 
  ggtitle("A. Upregulated") + 
  scale_x_continuous(limits = c(0.046, 0.17), breaks = seq(0.05, 0.15, by = 0.05), minor_breaks = seq(0.075, 0.175, by = 0.05)) +
  geom_point(stroke = 0.9) +
  geom_text(aes(label = enrichres_up_1vs6@result$significant, x = enrichres_up_1vs6@result$GeneRatioCONT, y = enrichres_up_1vs6$Description), size = 7, colour= "black", nudge_y = 0.2) +
  theme(text = element_text(size = 14, family = "Source Sans Pro", face="bold"),
        axis.text.x = element_text(size = 14, family = "Source Sans Pro", face="plain"), 
        axis.text.y = element_text(size = 12, family = "Source Sans Pro", face="plain"), 
        axis.title = element_text(size = 14, family = "Source Sans Pro", face="plain"), 
        legend.title = element_text(size=14, family = "Source Sans Pro", face="bold"),
        legend.text = element_text(size = 10, family = "Source Sans Pro", face="plain"), 
        panel.border = element_rect(colour = "#323232", fill = NA, linewidth=1))

ggsave(paste0(config$figures, "Fig_5A_up_metagenomic_1vs6.svg"), 
       width = 6,
       height = 5, 
       units = "in")

# Figure 5B
enrichres_down_1vs6@result$significant <- NA
enrichres_down_1vs6@result$significant[enrichres_down_1vs6@result$`p.adjust`<0.1] <- '*'
enrichres_down_1vs6@result$GeneRatioCONT <- sapply(parse(text = enrichres_down_1vs6@result$GeneRatio), eval)

p <- enrichplot::dotplot(enrichres_down_1vs6) + 
  ggtitle("B. Downregulated") +
  scale_x_continuous(limits = c(0.04, 0.265), breaks = seq(0.1, 0.2, by = 0.1), minor_breaks = seq(0.05, 0.25, by = 0.1)) +
  geom_point(stroke = 0.9) +
  geom_text(aes(label = enrichres_down_1vs6@result$significant, x = enrichres_down_1vs6@result$GeneRatioCONT, y = enrichres_down_1vs6$Description), size = 7, colour= "black", nudge_y = 0.35) +
  theme(text = element_text(size = 14, family = "Source Sans Pro", face="bold"),
        axis.text.x = element_text(size = 14, family = "Source Sans Pro", face="plain"), 
        axis.text.y = element_text(size = 12, family = "Source Sans Pro", face="plain"), 
        axis.title = element_text(size = 14, family = "Source Sans Pro", face="plain"), 
        legend.title = element_text(size=14, family = "Source Sans Pro", face="bold"),
        legend.text = element_text(size = 10, family = "Source Sans Pro", face="plain"), 
        panel.border = element_rect(colour = "#323232", fill = NA, linewidth=1))

ggsave(paste0(config$figures, "Fig_5B_down_metagenomic_1vs6.svg"), 
       width = 6.5,
       height = 5, 
       units = "in")

# Figure 5C
res_fmetab_down_1vs6$significant <- NA
res_fmetab_down_1vs6$significant[res_fmetab_down_1vs6$FDR<0.1] <- '*'

p <- ggplot(res_fmetab_down_1vs6, aes(x = minuslog10p, y = metaboliteSet)) +
  ggtitle("C. Underabundant (fecal)") +
  xlab("-log10(P)") +
  ylab("") +
  geom_point(aes(size = EnrichRatio, colour = FDR)) +
  geom_point(aes(size = EnrichRatio), fill = NA, shape = 1, stroke = 0.9, colour = "black") + 
  geom_text(aes(label = significant, x = minuslog10p, y = metaboliteSet), size = 7, colour= "black", nudge_y = 0.25) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) + 
  scale_x_continuous(limits = c(2.45, 5.2), breaks = seq(3, 5, by = 1), minor_breaks = seq(2.5, 4.5, by = 1)) +
  scale_color_gradient(low = "#357fb9", high = "#df6563", name = "q value", trans = 'reverse') +
  scale_size_continuous(name = "Enrich Ratio") +
  guides(color = guide_colorbar(order = 1), size = guide_legend(override.aes = list(fill = NA, shape = 1, colour = "black"), order = 2)) +
  theme_minimal() +
  theme(text = element_text(size = 14, family = "Source Sans Pro", face="bold", colour="#010101"),
        axis.text.x = element_text(size = 14, family = "Source Sans Pro", face="plain", colour="#010101"), 
        axis.text.y = element_text(size = 12, family = "Source Sans Pro", face="plain", colour="#010101"), 
        axis.title = element_text(size = 14, family = "Source Sans Pro", face="plain"), 
        legend.title = element_text(size=14, family = "Source Sans Pro", face="bold"),
        legend.text = element_text(size = 10, family = "Source Sans Pro", face="plain"), 
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), 
        axis.ticks = element_line())

ggsave(paste0(config$figures, "Fig_5C_down_fmetab_1vs6.svg"), 
       width = 6.5,
       height = 5, 
       units = "in")

# Figure 5D
res_pmetab_down_1vs6$significant <- NA
res_pmetab_down_1vs6$significant[res_pmetab_down_1vs6$FDR<0.1] <- '*'

p <- ggplot(res_pmetab_down_1vs6, aes(x = minuslog10p, y = metaboliteSet)) +
  ggtitle("D. Underabundant (plasma)") +
  xlab("-log10(P)") +
  ylab("") +
  geom_point(aes(size = EnrichRatio, colour = FDR)) +
  geom_point(aes(size = EnrichRatio), fill = NA, shape = 1, stroke = 0.9, colour = "black") + # black border
  geom_text(aes(label = significant, x = minuslog10p, y = metaboliteSet), size = 7, colour= "black", nudge_y = 0.32) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) + 
  scale_x_continuous(limits = c(3.1, 4.25), breaks = seq(3.25, 4.25, by = 0.5), minor_breaks = seq(3.5, 4, by = 0.5)) +
  scale_color_gradient(low ="#357fb9" , high = "#df6563", name = "q value", trans="reverse") +
  scale_size_continuous(name = "Enrich Ratio") +
  guides(color = guide_colorbar(order = 1), size = guide_legend(override.aes = list(fill = NA, shape = 1, colour = "black"), order = 2)) +
  theme_minimal() + 
  theme(text = element_text(size = 14, family = "Source Sans Pro", face="bold", colour="#010101"),
        axis.text.x = element_text(size = 14, family = "Source Sans Pro", face="plain", colour="#010101"), 
        axis.text.y = element_text(size = 12, family = "Source Sans Pro", face="plain", colour="#010101"), 
        axis.title = element_text(size = 14, family = "Source Sans Pro", face="plain"), 
        legend.title = element_text(size=14, family = "Source Sans Pro", face="bold"),
        legend.text = element_text(size = 10, family = "Source Sans Pro", face="plain"), 
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), 
        axis.ticks = element_line())


ggsave(paste0(config$figures, "Fig_5D_down_pmetab_1vs6.svg"), 
       width = 6.5,
       height = 5, 
       units = "in")
