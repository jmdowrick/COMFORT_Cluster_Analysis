# Metagenomics ----
load(paste0(config$data_results, "Fig_5AB_diff_expressed_genes.RData"))

res <- lapply(names(deg_1vs6),
              function(x) enrichKO(gene = deg_1vs6[[x]]$rowname,
                                   minGSSize = 5,
                                   universe = bg_genes))

# Convert the enrichment result for each sample pattern to a dataframe with the pathways
names(res) <- names(deg_1vs6)
res_df <- lapply(names(res), function(x) rbind(res[[x]]@result))
names(res_df) <- names(res)
res_df <- do.call(rbind, res_df)

res_df <- res_df %>% mutate(minuslog10padj = -log10(p.adjust))

# Subset so pathways meet criteria
res_df <- dplyr::filter(res_df, res_df$p.adjust < 0.05) 

# Select only upregulated in Cluster 1
res_df_clust1 <- res_df %>%
  rownames_to_column() 

res_df_clust1 <- res_df_clust1[str_detect(res_df_clust1$rowname, "DGBI"), ] %>%
  dplyr::select(-rowname)

# Perform enrichment analysis
enrichres_up_1vs6 <- new("enrichResult",
                        readable = FALSE,
                        result = res_df_clust1,
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "BH",
                        qvalueCutoff = 0.1,
                        organism = "human",
                        ontology = "UNKNOWN",
                        gene = df$rowname,
                        keytype = "UNKNOWN",
                        universe = unique(bg_genes),
                        gene2Symbol = character(0),
                        geneSets = as.list(bg_genes))

# Select only downregulated in Cluster 1
res_df_clust6 <- res_df %>%
  rownames_to_column() 

res_df_clust6 <- res_df_clust6[str_detect(res_df_clust6$rowname, "Healthy"), ] %>%
  dplyr::select(-rowname)

# Perform enrichment analysis
enrichres_down_1vs6 <- new("enrichResult",
                        readable = FALSE,
                        result = res_df_clust6,
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "BH",
                        qvalueCutoff = 0.1,
                        organism = "human",
                        ontology = "UNKNOWN",
                        gene = df$rowname,
                        keytype = "UNKNOWN",
                        universe = unique(bg_genes),
                        gene2Symbol = character(0),
                        geneSets = as.list(bg_genes))

# Fecal metabolomics ----
# Functional enrichment analysis of metabolome data was performed using metaboanalyst. The following sections restructure the datasets for plotting. 
res_fmetab_down_1vs6 <- read_csv(paste0(config$data_FEA,'1vs6_fmetab_enrich_path_down_customrefset.csv')) %>%
  dplyr::rename("metaboliteSet" = "...1") %>%
  mutate(minuslog10padj = -log10(FDR), minuslog10p = -log10(`Raw p`), EnrichRatio = hits/expected) %>%
  dplyr::filter(`Raw p` < 0.05) %>%
  arrange(`Raw p`) %>%
  head(10)

res_fmetab_down_1vs6$metaboliteSet <- factor(res_fmetab_down_1vs6$metaboliteSet, levels = res_fmetab_down_1vs6$metaboliteSet[order(res_fmetab_down_1vs6$`Raw p`, decreasing = TRUE)])

# Plasma metabolomics ----
res_pmetab_down_1vs6 <- read_csv(paste0(config$data_FEA,'1vs6_pmetab_enrich_path_down_customrefset.csv')) %>%
  dplyr::rename("metaboliteSet" = "...1") %>%
  mutate(minuslog10padj = -log10(FDR), minuslog10p = -log10(`Raw p`), EnrichRatio = hits/expected) %>%
  dplyr::filter(`Raw p` < 0.05) %>%
  arrange(`Raw p`) %>%
  head(10) 

res_pmetab_down_1vs6$metaboliteSet <- factor(res_pmetab_down_1vs6$metaboliteSet, levels = res_pmetab_down_1vs6$metaboliteSet[order(res_pmetab_down_1vs6$`Raw p`, decreasing = TRUE)])

# Save outputs ----
save(enrichres_up_1vs6, enrichres_down_1vs6, res_fmetab_down_1vs6, res_pmetab_down_1vs6, file=paste0(config$data_results, "Fig_5_FEA_results.RData"))
