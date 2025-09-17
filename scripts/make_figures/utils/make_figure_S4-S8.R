# import res files
load(file=paste0(config$data_results, "Fig_S4-S8_volcano_1vs6.RData"))

# Figure S4
EnhancedVolcano(res_pmetab_1vs6, pCutoff = 0.05,
                lab = rownames(res_pmetab_1vs6),
                x = 'log2FoldChange',
                y = 'pvalue')

ggsave(filename=paste0(config$sup_figures, "Fig_S4_pmetab_volcano_1vs6.svg"), 
       width = 8, 
       height = 6, 
       units = "in")

# Figure S5
EnhancedVolcano(res_amino_1vs6, pCutoff = 0.05,
                lab = rownames(res_amino_1vs6),
                x = 'log2FoldChange',
                y = 'pvalue')

ggsave(filename=paste0(config$sup_figures, "Fig_S5_amino_volcano_1vs6.svg"), 
       width = 8, 
       height = 6, 
       units = "in")

# Figure S6
EnhancedVolcano(res_fmetab_1vs6, pCutoff = 0.05,
                lab = rownames(res_fmetab_1vs6),
                x = 'log2FoldChange',
                y = 'pvalue')

ggsave(filename=paste0(config$sup_figures, "Fig_S6_fmetab_volcano_1vs6.svg"), 
       width = 8, 
       height = 6, 
       units = "in")

# Figure S7
EnhancedVolcano(res_bile_1vs6, pCutoff = 0.05,
                lab = rownames(res_bile_1vs6),
                x = 'log2FoldChange',
                y = 'pvalue')

ggsave(filename=paste0(config$sup_figures, "Fig_S7_bile_volcano_1vs6.svg"), 
       width = 8, 
       height = 6, 
       units = "in")

# Figure S8
EnhancedVolcano(res_scfa_1vs6, pCutoff = 0.05,
                lab = rownames(res_scfa_1vs6),
                x = 'log2FoldChange',
                y = 'pvalue')

ggsave(filename=paste0(config$sup_figures, "Fig_S8_scfa_volcano_1vs6.svg"), 
       width = 8, 
       height = 6, 
       units = "in")
