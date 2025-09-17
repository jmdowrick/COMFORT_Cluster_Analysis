library(readr)
library(tidyverse)
library(vegan)
library(lubridate)
library(caret)
library(dplyr)

config <- config::get()

mapping <- data.frame(read.delim(file=paste0(config$data_raw, config$file_metagenomics_metadata), sep="\t",stringsAsFactors=F))
mapping <- mapping[order(mapping[,1]),]
mapping[,1] <- gsub(".daa", "", mapping[,1])

L7 <- data.frame(read.delim(file=paste0(config$data_raw, config$file_metagenomics_taxa), sep="\t",stringsAsFactors=F))
rownames(L7) <- L7[,1]
L7 <- subset(L7, select=-X.Datasets)
colnames(L7) <- gsub("X","",colnames(L7)) 
L7 <- L7[rowSums(L7) > 0,]
L7 <- L7[,order(colnames(L7))]

# Remove repeat samples
L7 <- L7[,mapping$Repeat == "N"]
mapping <- mapping[mapping$Repeat == "N",]

# Change samples names (column names) to be the same as the participant names (e.g. "90001", "90002", etc)
colnames(L7) <- mapping$Participant

# Remove rows that only have zeros
L7 <- as.data.frame(L7[rowSums(L7) > 0,])
L7 <- as.data.frame(t(L7))

# Rarefy to minimum sequencing depth
rarefaction_depth <- min(rowSums(L7))

# Number of repeated rarefactions
num_repeats <- 100

# Store each iteration in case R crashes
save_path <- paste0(config$data_processed, "alpha-diversity/")
fils <- list.files(save_path, pattern="RData$", full.names = TRUE, recursive = TRUE)
num_files <- tibble(dir = dirname(fils)) %>% 
  dplyr::count(dir)

# Perform repeated rarefactions
if (num_files$n < num_repeats) {
  for (i in (num_files$n + 1):num_repeats) {
    file_name <- paste0('shannon_', as.character(i),'.RData')
    rarefied_data <- rrarefy(L7, sample = rarefaction_depth)
    shannon_diversity <- diversity(rarefied_data, index = "shannon")
    save(shannon_diversity, file = paste0(save_path,file_name))
  }
}