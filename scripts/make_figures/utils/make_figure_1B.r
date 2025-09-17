df_present <- read.csv(paste0(config$data_results, "Fig_1_data_availability.csv"), row.names = 1)
names(df_present) <- substring(names(df_present), 2)
df_present_numeric <- as.data.frame(lapply(df_present, as.numeric))

df_present$present <- apply(df_present_numeric, 1, function(x) length(which(x==1)))
df_present$partial <- apply(df_present_numeric, 1, function(x) length(which(x==0.7)))

df <- df_present %>%
  dplyr::select(c("present", "partial")) %>%
  rownames_to_column(var = "feature") %>%
  pivot_longer(cols=c("present", "partial")) 
    
df$feature <- factor(df$feature, levels=c("HADS", "PROMIS", "SAGIS", "gene", "taxa", "fmetab", "pmetab", "amino", "SCFA", "bile"))

p<-ggplot(df, aes(x = feature, y = value)) +
  geom_col(aes(fill=name)) +
  coord_flip() +
  theme_prism() +
  theme(panel.grid.major.x = element_line(color="gray"), panel.grid.minor.x = element_line(color="gray", linewidth = 0.5))

ggsave(paste0(config$figures, "Fig_1B_datacollection_count.svg"), 
plot = p,
width = 6,
height = 6, 
units = "in")