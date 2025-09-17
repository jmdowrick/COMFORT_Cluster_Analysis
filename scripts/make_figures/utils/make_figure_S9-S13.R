# import res files
load(file=paste0(config$data_results, "Fig_S9-S13_volcano_11vs6.RData"))

# Figure S9
EnhancedVolcano(res_pmetab_11vs6, pCutoff = 0.05,
                lab = rownames(res_pmetab_11vs6),
                x = 'log2FoldChange',
                y = 'pvalue')

ggsave(filename=paste0(config$sup_figures, "Fig_S9_pmetab_volcano_11vs6.svg"), 
       width = 8, 
       height = 6, 
       units = "in")

# Figure S10
EnhancedVolcano(res_amino_11vs6, pCutoff = 0.05,
                lab = rownames(res_amino_11vs6),
                x = 'log2FoldChange',
                y = 'pvalue')

ggsave(filename=paste0(config$sup_figures, "Fig_S10_amino_volcano_11vs6.svg"), 
       width = 8, 
       height = 6, 
       units = "in")

# Figure S11
EnhancedVolcano(res_fmetab_11vs6, pCutoff = 0.05,
                lab = rownames(res_fmetab_11vs6),
                x = 'log2FoldChange',
                y = 'pvalue')

ggsave(filename=paste0(config$sup_figures, "Fig_S11_fmetab_volcano_11vs6.svg"), 
       width = 8, 
       height = 6, 
       units = "in")

# Figure S12
EnhancedVolcano(res_bile_11vs6, pCutoff = 0.05,
                lab = rownames(res_bile_11vs6),
                x = 'log2FoldChange',
                y = 'pvalue')

ggsave(filename=paste0(config$sup_figures, "Fig_S12_bile_volcano_11vs6.svg"), 
       width = 8, 
       height = 6, 
       units = "in")

# Figure S13
EnhancedVolcano(res_scfa_11vs6, pCutoff = 0.05,
                lab = rownames(res_scfa_11vs6),
                x = 'log2FoldChange',
                y = 'pvalue')

ggsave(filename=paste0(config$sup_figures, "Fig_S13_scfa_volcano_11vs6.svg"), 
       width = 8, 
       height = 6, 
       units = "in")
