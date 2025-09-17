library(plyr)
library(tidyverse)
library(ALDEx2)
library(MicrobiomeProfiler)
library(Rvolcano)
library(EnhancedVolcano)
library(janitor)

set.seed(123)

config <- config::get()

abs_log <- function(x){
  x[x==0] <- 1
  si <- sign(x)
  return(si * log2(si*x))
}

# Import metagenomic KEGG data ----
mapping <- data.frame(read.delim(file=paste0(config$data_raw, config$file_metagenomics_metadata), sep="\t",stringsAsFactors=F))
mapping <- mapping[order(mapping[,1]),]
mapping[,1] <- gsub(".daa", "", mapping[,1])

L4.kegg <- data.frame(read.delim(file=paste0(config$data_raw, config$file_metagenomics_gene_mgseq), sep="\t",stringsAsFactors=F))
rownames(L4.kegg) <- L4.kegg[,1]
colnames(L4.kegg) <- gsub("X","",colnames(L4.kegg)) 
hierarchy <- data.frame(read.delim(file=paste0(config$data_raw, config$file_metagenomics_gene_heir), sep="\t",stringsAsFactors=F))
L4.kegg <- L4.kegg[,order(colnames(L4.kegg))]
L4.kegg[,"kegg"] <- hierarchy$L4
rownames(L4.kegg) <- L4.kegg$kegg
L4.kegg <- subset(L4.kegg, select=-kegg)
L4.kegg <- L4.kegg[rowSums(L4.kegg) > 0,]
L4.kegg <- L4.kegg[,order(colnames(L4.kegg))]

bg_genes <- rownames(L4.kegg) %>%
  stringr::str_extract("^.{6}") 

save(bg_genes, file=paste0(config$data_processed, "FEA_bg_genes.RData"))

L4.kegg <- as.data.frame(L4.kegg[!(row.names(L4.kegg) %in% "Not assigned"),]) # drop anything containing others or not assigned
L4.kegg <- L4.kegg |> filter(!str_detect(rownames(L4.kegg), "Other"))

# Remove repeat samples
L4.kegg <- L4.kegg[,mapping$Repeat == "N"]
mapping <- mapping[mapping$Repeat == "N",]

# Change samples names (column names) to be the same as the participant names (e.g. "90001", "90002", etc)
colnames(L4.kegg) <- mapping$Participant

L4.kegg <- L4.kegg[rowSums(L4.kegg) > 0,]

# Isolate groups of interest and create indicator array
L4.kegg <- as.data.frame(t(L4.kegg))
L4.kegg.filtered <- L4.kegg %>% tibble::rownames_to_column(var = "sample_id")

# Import remaining biological data ----
df_f_hilic <- readxl::read_excel(paste0(config$data_raw, "COMFORT_Faecal_Metabolomics_All_data_Final_cleaned_names.xlsx"), sheet = config$sheet_fmetab_HILIC) %>%
  dplyr::select(c("Participant #", 12:ncol(.))) %>%
  as.data.frame()

df_f_c18 <- readxl::read_excel(paste0(config$data_raw, "COMFORT_Faecal_Metabolomics_All_data_Final_cleaned_names.xlsx"), sheet = config$sheet_fmetab_C18) %>%
  dplyr::select(c("Participant #", 12:ncol(.))) %>%
  as.data.frame()

df_f_lipid <- readxl::read_excel(paste0(config$data_raw, "COMFORT_Faecal_Metabolomics_All_data_Final_cleaned_names.xlsx"), sheet = config$sheet_fmetab_Lipid) %>%
  dplyr::select(c("Participant #", 12:ncol(.))) %>%
  as.data.frame()

comfort_f_metab <- merge(df_f_lipid, df_f_c18, by = "Participant #", all = TRUE) %>%
  merge(., df_f_hilic, by = "Participant #", all = TRUE) %>%
  dplyr::rename("id" = "Participant #")

rownames(comfort_f_metab) <- comfort_f_metab[, c(1)]
comfort_f_metab[, 1] <- NULL
comfort_f_metab <- as.data.frame(t(comfort_f_metab))

df_fmetab <- as.data.frame(t(comfort_f_metab))
df_fmetab <- df_fmetab %>% tibble::rownames_to_column(var = "sample_id")

comfort_bile <- readxl::read_excel(paste0(config$data_raw, config$file_bile_acid), sheet = config$sheet_bile_acid) %>%
  row_to_names(row_number = 1) %>%
  dplyr::select(c("Participant #", "GHDCA":"TCA")) %>%
  dplyr::rename("sample_id" = "Participant #") %>%
  as.data.frame() %>%
  mutate_all(function(x) as.numeric(as.character(x))) 
  
comfort_bile$sample_id <- as.character(comfort_bile$sample_id)

df <- readxl::read_excel(paste0(config$data_raw, config$file_organic_acid), sheet=config$sheet_organic_acid)

thresholds <- df[14, 13:26] %>%
  as.data.frame() %>%
  mutate_all(function(x) as.numeric(as.character(x)))

comfort_organic <- readxl::read_excel(paste0(config$data_raw, config$file_organic_acid), sheet = config$sheet_organic_acid) %>%
  row_to_names(row_number = 15) %>%
  dplyr::select(c("Sample Client ID", "FA":"HA")) %>%
  clean_names %>%
  dplyr::rename("sample_id" = "sample_client_id") %>%
  as.data.frame() %>%
  mutate_all(function(x) as.numeric(as.character(x)))

comfort_organic$sample_id <- as.character(comfort_organic$sample_id)

comfort_organic[, -c(1)] <- comfort_organic[, -c(1)] * (comfort_organic[, -c(1)] > thresholds[col(comfort_organic[, -c(1)])])

comfort_organic <- dplyr::select(comfort_organic, -c("x3mva")) %>%
  na.omit()

rm(thresholds, df)

pmetab <- readxl::read_excel(paste0(config$data_raw, config$file_pmetab), sheet = "All Data combined_knowns_final") %>%
  dplyr::select(c("Participant #", 19:ncol(.))) %>%
  as.data.frame() %>%
  dplyr::rename(., "sample_id" = "Participant #") %>%
  mutate_all(function(x) as.numeric(as.character(x))) %>%
  as.data.frame()

pmetab$sample_id <- as.character(pmetab$sample_id)

amino <- readxl::read_excel(paste0(config$data_raw, config$file_amino_acid), sheet = "Combined results") %>%
  row_to_names(row_number = 1) %>%
  dplyr::select(c("Participant #", "Aspartic Acid":"Trytophan")) %>%
  dplyr::rename("sample_id" = "Participant #") %>%
  mutate_all(function(x) as.numeric(as.character(x))) %>%
  as.data.frame()
amino$sample_id <- as.character(amino$sample_id)

# Import integrated clustering labels ----
merged_labels <- read.csv(paste0(config$data_cluster_assignments, "integrated_clusters.csv")) %>%
  dplyr::rename(sample_id = X)
  
merged_labels$sample_id <- as.character(merged_labels$sample_id)
merged_labels$ref_clust <- as.factor(merged_labels$ref_clust)

# Run differential expression analysis ----
source("analysis/scripts/run_analysis/utils/DE_1-6.r")
source("analysis/scripts/run_analysis/utils/DE_7-6.r")
source("analysis/scripts/run_analysis/utils/DE_8-6.r")
source("analysis/scripts/run_analysis/utils/DE_11-6.r")
