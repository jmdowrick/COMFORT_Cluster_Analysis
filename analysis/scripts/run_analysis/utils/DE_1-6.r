# Differential expression analysis for Cluster 1 vs 6

toCompare <- c(1, 6) 
diseaseIndex <- 1
healthyIndex <- 2

combined <- merged_labels %>%
  full_join(L4.kegg.filtered, by = "sample_id") 

# Metagenomics ----
combined <- combined %>%
  dplyr::filter(ref_clust %in% toCompare) %>%  
  droplevels() %>%
  dplyr::rename(disease_state = ref_clust) %>%
  drop_na()

levels(combined$disease_state)[diseaseIndex] <- 'disease'
levels(combined$disease_state)[healthyIndex] <- 'healthy'
combined$disease_state <- as.character(combined$disease_state)

diagnosis <- combined$disease_state

combined$disease_state <- NULL
combined <- combined %>%
  column_to_rownames("sample_id")

combined <- as.data.frame(t(combined)) 

aldex.comp <- aldex(combined, diagnosis, mc.samples=128, test="t", effect=TRUE, 
                    include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE, gamma=NULL)

# Annotate genes with differential expression 
df <- aldex.comp %>%
  rownames_to_column() %>%
  dplyr::select(rowname, we.eBH, effect) %>%
  mutate(diffexpressed = case_when(
    we.eBH < 0.05 & effect > 0 ~ 'Healthy',
    we.eBH < 0.05 & effect < 0 ~ 'DGBI',
    we.eBH > 0 ~ 'NO'
  ))

df[,'rowname'] <- df[,'rowname'] %>%
  stringr::str_extract("^.{6}")

# Remove non-significant genes
df <- df[df$diffexpressed != 'NO',]

# Split dataframe into upregulated and downregulated genes
deg_1vs6 <- split(df, df$diffexpressed)

# Set up faecal biomarkers for differential analysis ----
# Faecal untargeted metabolomics 
fmetab <- merged_labels %>%
  full_join(df_fmetab, by = "sample_id") 

fmetab <- fmetab %>%
  dplyr::filter(ref_clust %in% toCompare) %>%  
  droplevels() %>%
  dplyr::rename(disease_state = ref_clust) %>%
  drop_na() 

fmetab <- fmetab[order(fmetab$disease_state),]
rownames(fmetab) <- NULL

levels(fmetab$disease_state)[diseaseIndex] <- 'disease'
levels(fmetab$disease_state)[healthyIndex] <- 'healthy'
fmetab$disease_state <- make.unique(as.character(fmetab$disease_state))

fmetab <- fmetab %>%
  column_to_rownames("disease_state")
fmetab$sample_id <- NULL
fmetab <- as.data.frame(t(fmetab)) 

fmetab <- fmetab %>%
  mutate_all(function(x) as.numeric(as.character(x))) %>%
  select_if(~ !any(is.na(.))) %>%
  select_if( is.numeric)%>%
  dplyr::select(sort(names(.), decreasing = TRUE))

# Bile acid 
bile <- merged_labels %>%
  full_join(comfort_bile, by = "sample_id") %>%
  dplyr::filter(ref_clust %in% toCompare) %>%
  droplevels() %>%
  dplyr::rename(disease_state = ref_clust) %>%
  drop_na()

bile <- bile[order(bile$disease_state),]
rownames(bile) <- NULL

levels(bile$disease_state)[diseaseIndex] <- 'disease'
levels(bile$disease_state)[healthyIndex] <- 'healthy'
bile$disease_state <- make.unique(as.character(bile$disease_state))

bile <- bile %>%
  column_to_rownames("disease_state")
bile$sample_id <- NULL
bile <- as.data.frame(t(bile)) 

bile <- bile %>%
  mutate_all(function(x) as.numeric(as.character(x))) %>%
  select_if(~ !any(is.na(.))) %>%
  select_if(is.numeric)%>%
  dplyr::select(sort(names(.), decreasing = TRUE))

# SCFA 
scfa <- merged_labels %>%
  full_join(comfort_organic, by = "sample_id") %>%
  dplyr::filter(ref_clust %in% toCompare) %>%
  droplevels() %>%
  dplyr::rename(disease_state = ref_clust) %>%
  drop_na()

scfa <- scfa[order(scfa$disease_state),]
rownames(scfa) <- NULL

levels(scfa$disease_state)[diseaseIndex] <- 'disease'
levels(scfa$disease_state)[healthyIndex] <- 'healthy'
scfa$disease_state <- make.unique(as.character(scfa$disease_state))

scfa <- scfa %>%
  column_to_rownames("disease_state")
scfa$sample_id <- NULL
scfa <- as.data.frame(t(scfa)) 

scfa <- scfa %>%
  mutate_all(function(x) as.numeric(as.character(x))) %>%
  select_if(~ !any(is.na(.))) %>%
  select_if(is.numeric) %>%
  dplyr::select(sort(names(.), decreasing = TRUE))

# Perform Robust Volcano plot comparison for faecal biomarkers ----

# Faecal untargeted metabolomics 
FC_fecal_1vs6 <- foldChngCalc(fmetab, 16, 35) 
p_fecal_1v6 <- NULL
for (i in 1:dim(fmetab)[1]){
  p_fecal_1v6[i]<-p.valcalc(as.numeric(fmetab[i,1:16]),as.numeric(fmetab[i,17:51]))
}

p_fecal_1v6 <- p.adjust(p_fecal_1v6, "BH", n = length(p_fecal_1v6) + dim(scfa)[1] + dim(bile)[1])

res_fmetab_1vs6 <- cbind(FC_fecal_1vs6, p_fecal_1v6) %>%
  as.data.frame() %>%
  dplyr::rename(log2FoldChange = FC_fecal_1vs6, pvalue = p_fecal_1v6)

res_fmetab_1vs6$log2FoldChange <- abs_log(res_fmetab_1vs6$log2FoldChange)
output_f_untargeted <- EnhancedVolcano(res_fmetab_1vs6, pCutoff = 0.05,
                                       lab = rownames(res_fmetab_1vs6),
                                       x = 'log2FoldChange',
                                       y = 'pvalue')

# Bile acid 
FC_bile_1vs6 <- foldChngCalc(bile, 16, 44) 

p_bile_1vs6 <- NULL
for (i in 1:dim(bile)[1]){
  p_bile_1vs6[i]<-p.valcalc(as.numeric(bile[i,1:16]),as.numeric(bile[i,17:60]))
}

p_bile_1vs6 <- p.adjust(p_bile_1vs6, "BH", n = length(p_bile_1vs6) + dim(scfa)[1] + dim(bile)[1])

res_bile_1vs6 <- cbind(FC_bile_1vs6, p_bile_1vs6) %>%
  as.data.frame() %>%
  dplyr::rename(log2FoldChange = FC_bile_1vs6, pvalue = p_bile_1vs6)


res_bile_1vs6$log2FoldChange <- abs_log(res_bile_1vs6$log2FoldChange)

output_f_bile <- EnhancedVolcano(res_bile_1vs6, pCutoff = 0.05,
                                 lab = rownames(res_bile_1vs6),
                                 x = 'log2FoldChange',
                                 y = 'pvalue')

# SCFA
FC_scfa_1vs6 <- foldChngCalc(scfa, 21, 52) 

p_scfa_1vs6 <- NULL
for (i in 1:dim(scfa)[1]){
  p_scfa_1vs6[i]<-p.valcalc(as.numeric(scfa[i,1:21]),as.numeric(scfa[i,22:73]))
}

p_scfa_1vs6 <- p.adjust(p_scfa_1vs6, "BH", n = length(p_bile_1vs6) + dim(scfa)[1] + dim(bile)[1])

res_scfa_1vs6 <- cbind(FC_scfa_1vs6, p_scfa_1vs6) %>%
  as.data.frame() %>%
  dplyr::rename(log2FoldChange = FC_scfa_1vs6, pvalue = p_scfa_1vs6)

res_scfa_1vs6$log2FoldChange <- abs_log(res_scfa_1vs6$log2FoldChange)

output_f_scfa <- EnhancedVolcano(res_scfa_1vs6, pCutoff = 0.05,
                                 lab = rownames(res_scfa_1vs6),
                                 x = 'log2FoldChange',
                                 y = 'pvalue')

# Export list of differentially expressed fecal metabolites ----

# Untargeted fecal metabolites
fmetab_foldchange <- output_f_untargeted$data$log2FoldChange[output_f_untargeted$data$Sig == "FC_P"] %>%
  as.data.frame() 
colnames(fmetab_foldchange) <- c("FC")
rownames(fmetab_foldchange) <-  output_f_untargeted$data$lab[output_f_untargeted$data$Sig == "FC_P"] 
fmetab_foldchange <- fmetab_foldchange %>%
  rownames_to_column("commonname")

# Bile
fmetab_foldchange_bile <- output_f_bile$data$log2FoldChange[output_f_bile$data$Sig == "FC_P"] %>%
  as.data.frame() 
colnames(fmetab_foldchange_bile) <- c("FC")
rownames(fmetab_foldchange_bile) <-  output_f_bile$data$lab[output_f_bile$data$Sig == "FC_P"] 
fmetab_foldchange_bile <- fmetab_foldchange_bile %>%
  rownames_to_column("commonname")

# SCFA
fmetab_foldchange_scfa <- output_f_scfa$data$log2FoldChange[output_f_scfa$data$Sig == "FC_P"] %>%
  as.data.frame() 
colnames(fmetab_foldchange_scfa) <- c("FC")
rownames(fmetab_foldchange_scfa) <-  output_f_scfa$data$lab[output_f_scfa$data$Sig == "FC_P"] 
fmetab_foldchange_scfa <- fmetab_foldchange_scfa %>%
  rownames_to_column("commonname")

fmetab_foldchange <- fmetab_foldchange %>%
  full_join(., fmetab_foldchange_scfa, by = c("commonname", "FC")) %>%
  full_join(., fmetab_foldchange_bile, by = c("commonname", "FC")) 

fmetab_foldchange <- fmetab_foldchange[!(fmetab_foldchange$commonname %in% colnames(df_f_lipid)), ]

# Export fecal metabolite list for metaboanalyst 
writexl::write_xlsx(fmetab_foldchange, paste0(config$data_results, "Fig_5_fmetab_1vs6_c18_hilic_targeted.xlsx"))

# Export fecal lipid metabolite list for metaboanalyst
fmetab_foldchange <- output_f_untargeted$data$log2FoldChange[output_f_untargeted$data$Sig == "FC_P"] %>%
  as.data.frame() 
colnames(fmetab_foldchange) <- c("FC")
rownames(fmetab_foldchange) <-  output_f_untargeted$data$lab[output_f_untargeted$data$Sig == "FC_P"] 
fmetab_foldchange <- fmetab_foldchange %>%
  rownames_to_column("commonname") 
fmetab_foldchange <- fmetab_foldchange[fmetab_foldchange$commonname %in% colnames(df_f_lipid), ]
fmetab_foldchange$commonname <- fmetab_foldchange$commonname %>%
  gsub("+HCOO", "", ., fixed = T) %>%
  gsub("+NH4", "", ., fixed = T) %>%
  gsub("+H", "", ., fixed = T) %>%
  gsub("-H", "", ., fixed = T) %>%
  gsub("+Na", "", ., fixed = T) %>%
  gsub("+O", "", ., fixed = T) 
writexl::write_xlsx(fmetab_foldchange,paste0(config$data_results, "Fig_5_fmetab_1vs6_lipids.xlsx"))

# Set up plasma biomarkers for differential analysis ----

# Plasma untargeted metabolomics 
pmetab_lab <- merged_labels %>%
  full_join(pmetab, by = "sample_id") %>%
  dplyr::filter(ref_clust %in% toCompare) %>%  
  droplevels() %>%
  dplyr::rename(disease_state = ref_clust) %>%
  drop_na()

pmetab_lab <- pmetab_lab[order(pmetab_lab$disease_state, decreasing = TRUE),]
rownames(pmetab_lab) <- NULL

levels(pmetab_lab$disease_state)[diseaseIndex] <- 'disease'
levels(pmetab_lab$disease_state)[healthyIndex] <- 'healthy'
pmetab_lab$disease_state <- make.unique(as.character(pmetab_lab$disease_state))

pmetab_lab <- pmetab_lab %>%
  column_to_rownames("disease_state")
pmetab_lab$sample_id <- NULL
pmetab_lab <- as.data.frame(t(pmetab_lab)) 

# Amino acids
amino_lab <- merged_labels %>%
  full_join(amino, by = "sample_id") %>%
  dplyr::filter(ref_clust %in% toCompare) %>%  
  droplevels() %>%
  dplyr::rename(disease_state = ref_clust) %>%
  drop_na()

amino_lab <- amino_lab[order(amino_lab$disease_state, decreasing = TRUE),]
rownames(amino_lab) <- NULL

levels(amino_lab$disease_state)[diseaseIndex] <- 'disease'
levels(amino_lab$disease_state)[healthyIndex] <- 'healthy'
amino_lab$disease_state <- make.unique(as.character(amino_lab$disease_state))

amino_lab <- amino_lab %>%
  column_to_rownames("disease_state")
amino_lab$sample_id <- NULL
amino_lab <- as.data.frame(t(amino_lab)) 

# Perform Robust Volcano plot comparison for plasma biomarkers ----

# Untargeted
FC_pmetab_1vs6 <- foldChngCalc(pmetab_lab, 15, 36) 

p_pmetab_1vs6 <- NULL
for (i in 1:dim(pmetab_lab)[1]){
  p_pmetab_1vs6[i]<-p.valcalc(as.numeric(pmetab_lab[i,1:15]),as.numeric(pmetab_lab[i,16:51]))
}

p_pmetab_1vs6 <- p.adjust(p_pmetab_1vs6, "BH", n = length(p_pmetab_1vs6) + dim(amino_lab)[1])
res_pmetab_1vs6 <- cbind(FC_pmetab_1vs6, p_pmetab_1vs6) %>%
  as.data.frame() %>%
  dplyr::rename(log2FoldChange = FC_pmetab_1vs6, pvalue = p_pmetab_1vs6)

res_pmetab_1vs6$log2FoldChange <- abs_log(res_pmetab_1vs6$log2FoldChange)

output_p_untargeted <- EnhancedVolcano(res_pmetab_1vs6, pCutoff = 0.05,
                                       lab = rownames(res_pmetab_1vs6),
                                       x = 'log2FoldChange',
                                       y = 'pvalue')

# Amino
FC_pmetab_1vs6_amino <- foldChngCalc(amino_lab, 11, 40) 

p_pmetab_1vs6_amino <- NULL
for (i in 1:dim(amino_lab)[1]){
  p_pmetab_1vs6_amino[i]<-p.valcalc(as.numeric(amino_lab[i,1:11]),as.numeric(amino_lab[i,12:51]))
}

p_pmetab_1vs6_amino <- p.adjust(p_pmetab_1vs6_amino, "BH", n = length(p_pmetab_1vs6_amino) + length(res_pmetab_1vs6$pvalue))

res_amino_1vs6 <- cbind(FC_pmetab_1vs6_amino, p_pmetab_1vs6_amino) %>%
  as.data.frame() %>%
  dplyr::rename(log2FoldChange = FC_pmetab_1vs6_amino, pvalue = p_pmetab_1vs6_amino)

res_amino_1vs6$log2FoldChange <- abs_log(res_amino_1vs6$log2FoldChange)

output_p_amino <- EnhancedVolcano(res_amino_1vs6, pCutoff = 0.05,
                                  lab = rownames(res_amino_1vs6),
                                  x = 'log2FoldChange',
                                  y = 'pvalue')

# Export differentially abundant plasma metabolites
pmetab_foldchange <- output_p_untargeted$data$log2FoldChange[output_p_untargeted$data$Sig == "FC_P"] %>%
  as.data.frame() 
colnames(pmetab_foldchange) <- c("FC")
rownames(pmetab_foldchange) <-  output_p_untargeted$data$lab[output_p_untargeted$data$Sig == "FC_P"] 
pmetab_foldchange <- pmetab_foldchange %>%
  rownames_to_column("commonname") 
pmetab_foldchange$commonname <- pmetab_foldchange$commonname %>%
  gsub("LP_", "", ., fixed = T) %>%
  gsub("CN_", "", ., fixed = T) %>%
  gsub("LN_", "", ., fixed = T) %>%
  gsub("HP_", "", ., fixed = T) %>%
  gsub("HN_", "", ., fixed = T) %>%
  gsub("+Na_MS2", "", ., fixed = T) %>%
  gsub("+NH4_MS2", "", ., fixed = T) %>%
  gsub("; [M+NH4]+", "", ., fixed = T) %>%
  gsub("+H_MS2", "", ., fixed = T) 

pmetab_foldchange_amino <- output_p_amino$data$log2FoldChange[output_p_amino$data$Sig == "FC_P"] %>%
  as.data.frame() 
colnames(pmetab_foldchange_amino) <- c("FC")
rownames(pmetab_foldchange_amino) <-  output_p_amino$data$lab[output_p_amino$data$Sig == "FC_P"] 
pmetab_foldchange_amino <- pmetab_foldchange_amino %>%
  rownames_to_column("commonname")

pmetab_foldchange <- pmetab_foldchange %>%
  full_join(., pmetab_foldchange_amino, by = c("commonname", "FC")) 

writexl::write_xlsx(pmetab_foldchange, paste0(config$data_results, "Fig_5_pmetab_1vs6.xlsx"))

# Save outputs to file ----
save(deg_1vs6, df, file=paste0(config$data_results, "Fig_5AB_diff_expressed_genes.RData"))
save(res_pmetab_1vs6, res_amino_1vs6, res_fmetab_1vs6, res_bile_1vs6, res_scfa_1vs6, file=paste0(config$data_results, "Fig_S4-S8_volcano_1vs6.RData"))
