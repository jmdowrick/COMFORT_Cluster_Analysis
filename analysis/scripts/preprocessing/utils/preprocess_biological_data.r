library(plyr)
library(tidyverse)
library(janitor)
library(readr)
library(vegan)
library(lubridate)
library(caret)

config <- config::get()

# Amino acids ----
comfort_amino <- readxl::read_excel(paste0(config$data_raw, config$file_amino_acid), sheet = config$sheet_amino_acid) %>%
  row_to_names(row_number = 1) %>%
  dplyr::select(c("Participant #", "Aspartic Acid":"Trytophan")) %>%
  clean_names %>%
  dplyr::rename("id" = "participant_number") %>%
  as.data.frame() %>%
  dplyr::mutate_all(function(x) as.numeric(as.character(x)))

comfort_amino <- na.omit(comfort_amino)

rownames(comfort_amino) <- comfort_amino[, 1]
comfort_amino[, 1] <- NULL
comfort_amino <- as.data.frame(t(comfort_amino))

write.csv(comfort_amino, paste0(config$data_processed, "comfort_amino_acids_cluster.csv"), row.names = TRUE)

# Bile acids ----
comfort_bile <- readxl::read_excel(paste0(config$data_raw, config$file_bile_acid), sheet = config$sheet_bile_acid) %>%
  row_to_names(row_number = 1) %>%
  dplyr::select(c("Participant #", "GHDCA":"TCA")) %>%
  clean_names %>%
  dplyr::rename("id" = "participant_number") %>%
  as.data.frame() %>%
  dplyr::mutate_all(function(x) as.numeric(as.character(x)))

rownames(comfort_bile) <- comfort_bile[, 1]
comfort_bile[, 1] <- NULL
comfort_bile <- as.data.frame(t(comfort_bile))

write.csv(comfort_bile, paste0(config$data_processed, "comfort_bile_acids_cluster.csv"), row.names = TRUE)

# Fecal metabolomics ----
df_f_hilic <- readxl::read_excel(paste0(config$data_raw, config$file_fmetab), sheet = config$sheet_fmetab_HILIC) %>%
  dplyr::select(c("Participant #", 12:ncol(.))) %>%
  as.data.frame()

df_f_c18 <- readxl::read_excel(paste0(config$data_raw, config$file_fmetab), sheet = config$sheet_fmetab_C18) %>%
  dplyr::select(c("Participant #", 12:ncol(.))) %>%
  as.data.frame()

df_f_lipid <- readxl::read_excel(paste0(config$data_raw, config$file_fmetab), sheet = config$sheet_fmetab_Lipid) %>%
  dplyr::select(c("Participant #", 12:ncol(.))) %>%
  as.data.frame()

comfort_f_metab <- merge(df_f_lipid, df_f_c18, by = "Participant #", all = TRUE) %>%
  merge(., df_f_hilic, by = "Participant #", all = TRUE) %>%
  clean_names %>%
  dplyr::rename("id" = "participant_number")

rownames(comfort_f_metab) <- comfort_f_metab[, c(1)]
comfort_f_metab[, 1] <- NULL
comfort_f_metab <- as.data.frame(t(comfort_f_metab))

write.csv(comfort_f_metab, paste0(config$data_processed, "comfort_fmetab_cluster.csv"), row.names = TRUE)

# Metagenomics ----
zeroFilter <- function(x, threshold = 0.8) {
  percent0 <- apply(x, 1, function(x) {
    sum(x == 0) / length(x)
  }
  )
  x <- x[percent0 < threshold,]
  return(x)
}

# Taxonomic composition 
mapping <- data.frame(read.delim(file=paste0(config$data_raw, config$file_metagenomics_metadata), sep="\t",stringsAsFactors=F))
mapping <- mapping[order(mapping[,1]),]
mapping[,1] <- gsub(".daa", "", mapping[,1])

# Species level
L7 <- data.frame(read.delim(file=paste0(config$data_raw, config$file_metagenomics_taxa), sep="\t",stringsAsFactors=F))
rownames(L7) <- L7[,1]
L7 <- subset(L7, select=-X.Datasets)
colnames(L7) <- gsub("X","",colnames(L7)) 
L7 <- L7[rowSums(L7) > 0,]
L7 <- L7[,order(colnames(L7))]
L7 <- prop.table(as.matrix(L7),2)

# Remove repeat samples
L7 <- L7[,mapping$Repeat == "N"]

mapping <- mapping[mapping$Repeat == "N",]

# Change samples names (column names) to be the same as the participant names (e.g. "90001", "90002", etc)
colnames(L7) <- mapping$Participant

# Remove rows that only have zeros
L7 <- as.data.frame(L7[rowSums(L7) > 0,])

# Discard mean abundance <0.005% 
L7.filtered <- L7 %>% filter(rowMeans(.) > 0.00005)

# Discard taxa with >90% 0s
L7.filtered <- zeroFilter(L7.filtered, threshold = 0.9)

# Remove unassigned or contains environmental
to_remove = "root;Not assigned;"
L7.filtered <- as.data.frame(L7.filtered[!(row.names(L7.filtered) %in% to_remove),])
L7.filtered <- L7.filtered |> filter(!str_detect(rownames(L7.filtered), "environment"))

# Gene abundance 
mapping <- data.frame(read.delim(file=paste0(config$data_raw, config$file_metagenomics_metadata), sep="\t",stringsAsFactors=F))
mapping <- mapping[order(mapping[,1]),]
mapping[,1] <- gsub(".daa", "", mapping[,1])

# L3
L3.kegg <- data.frame(read.delim(file=paste0(config$data_raw, config$file_metagenomics_gene_mgseq), sep="\t",stringsAsFactors=F))
rownames(L3.kegg) <- L3.kegg[,1]
colnames(L3.kegg) <- gsub("X","",colnames(L3.kegg)) 
hierarchy <- data.frame(read.delim(file=paste0(config$data_raw, config$file_metagenomics_gene_heir), sep="\t",stringsAsFactors=F))
L3.kegg <- L3.kegg[,order(colnames(L3.kegg))]
L3.kegg[,"kegg"] <- hierarchy$L3
L3.kegg <- ddply(L3.kegg, "kegg", numcolwise(sum))
rownames(L3.kegg) <- L3.kegg$kegg
L3.kegg <- subset(L3.kegg, select=-kegg)
L3.kegg <- L3.kegg[rowSums(L3.kegg) > 0,]
L3.kegg <- L3.kegg[,order(colnames(L3.kegg))]
L3.kegg <- prop.table(as.matrix(L3.kegg),2)

L3.kegg <- as.data.frame(L3.kegg[!(row.names(L3.kegg) %in% "Not assigned"),]) # drop anything containing others or not assigned
L3.kegg <- L3.kegg |> filter(!str_detect(rownames(L3.kegg), "Other"))

# Remove repeat samples
L3.kegg <- L3.kegg[,mapping$Repeat == "N"]

mapping <- mapping[mapping$Repeat == "N",]

# Change samples names (column names) to be the same as the participant names (e.g. "90001", "90002", etc)
colnames(L3.kegg) <- mapping$Participant

# Remove rows that only have zeros
L3.kegg <- L3.kegg[rowSums(L3.kegg) > 0,]

# Discard mean abundance <0.005% 
L3.kegg.filtered <- L3.kegg %>% filter(rowMeans(.) > 0.00005)

# Discard pathways with >90% 0s
L3.kegg.filtered <- zeroFilter(L3.kegg.filtered, threshold = 0.9)

# Save output 
write.csv(L7.filtered, paste0(config$data_processed, "comfort_taxa_cluster.csv"), row.names = TRUE)
write.csv(L3.kegg.filtered, paste0(config$data_processed, "comfort_gene_cluster.csv"), row.names = TRUE)

# Organic acids ----
df <- readxl::read_excel(paste0(config$data_raw, config$file_organic_acid), sheet = config$sheet_organic_acid)

threshholds <- df[14, 13:26] %>%
  as.data.frame() %>%
  mutate_all(function(x) as.numeric(as.character(x)))

comfort_organic <- readxl::read_excel(paste0(config$data_raw, config$file_organic_acid), sheet = config$sheet_organic_acid) %>%
  row_to_names(row_number = 15) %>%
  dplyr::select(c("Sample Client ID", "FA":"HA")) %>%
  clean_names %>%
  dplyr::rename("id" = "sample_client_id") %>%
  as.data.frame() %>%
  mutate_all(function(x) as.numeric(as.character(x)))

comfort_organic[, -c(1)] <- comfort_organic[, -c(1)] * (comfort_organic[, -c(1)] > threshholds[col(comfort_organic[, -c(1)])])

comfort_organic <- dplyr::select(comfort_organic, -c("x3mva")) %>%
  na.omit()

rownames(comfort_organic) <- comfort_organic[, 1]
comfort_organic[, 1] <- NULL
comfort_organic <- as.data.frame(t(comfort_organic))

write.csv(comfort_organic, paste0(config$data_processed, "comfort_organic_acids_cluster.csv"), row.names = TRUE)

# Plasma metabolomics ---- 
df_p_hilic <- readxl::read_excel(paste0(config$data_raw, config$file_pmetab), sheet = config$sheet_pmetab_HILIC) %>%
  dplyr::select(c("Participant #", 19:ncol(.))) %>%
  as.data.frame() 
paste0(config$data_raw, config$file_bile_acid)

df_p_c18 <- readxl::read_excel(paste0(config$data_raw, config$file_pmetab), sheet = config$sheet_pmetab_C18) %>%
  dplyr::select(c("Participant #", 19:ncol(.))) %>%
  as.data.frame()

df_p_lipid <- readxl::read_excel(paste0(config$data_raw, config$file_pmetab), sheet = config$sheet_pmetab_Lipid) %>%
  dplyr::select(c("Participant #", 19:ncol(.))) %>%
  as.data.frame()
  
comfort_p_metab <- merge(df_p_lipid, df_p_c18, by = "Participant #", all = TRUE) %>%
  merge(., df_p_hilic, by = "Participant #", all = TRUE) %>%
  clean_names %>%
  dplyr::rename(., "id" = "participant_number") %>%
  mutate_all(function(x) as.numeric(as.character(x))) %>%
  select_if(~ !any(is.na(.)))

rownames(comfort_p_metab) <- comfort_p_metab[, c(1)]
comfort_p_metab[, 1] <- NULL
comfort_p_metab <- as.data.frame(t(comfort_p_metab))

write.csv(comfort_p_metab, paste0(config$data_processed, "comfort_pmetab_cluster.csv"), row.names = TRUE)
