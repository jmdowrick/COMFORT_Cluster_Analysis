library(readr)
library(tidyverse)
library(vegan)
library(lubridate)
library(caret)
library(dplyr)
library(rlist)

config <- config::get()

# Firmicutes-Bacteroidetes phyla ratio ----

zeroFilter <- function(x, threshold = 0.8) {
  percent0 <- apply(x, 1, function(x) {
    sum(x == 0) / length(x)
  }
  )
  x <- x[percent0 < threshold, ]
  return(x)
}


mapping <- data.frame(read.delim(file=paste0(config$data_raw, config$file_metagenomics_metadata), sep="\t",stringsAsFactors=F))
mapping <- mapping[order(mapping[,1]),]
mapping[,1] <- gsub(".daa", "", mapping[,1])

# Phylum level
L2 <- data.frame(read.delim(file=paste0(config$data_raw, "metagenomics/HVN_EG_B1234_TAXA_PHYLUM.txt"), sep="\t",stringsAsFactors=F))
rownames(L2) <- L2[,1]
L2 <- subset(L2, select=-X.Datasets)
colnames(L2) <- gsub("X","",colnames(L2)) 
L2 <- L2[rowSums(L2) > 0,]
L2 <- L2[,order(colnames(L2))]
L2 <- prop.table(as.matrix(L2),2)

# Remove repeat samples
L2 <- L2[,mapping$Repeat == "N"]
mapping <- mapping[mapping$Repeat == "N",]

# Change samples names (column names) to be the same as the participant names (e.g. "90001", "90002", etc)
colnames(L2) <- mapping$Participant

# Remove rows that only have zeros
L2 <- as.data.frame(L2[rowSums(L2) > 0,])

# Discard mean abundance <0.005% 
L2.filtered <- L2 %>% filter(rowMeans(.) > 0.00005)

# Discard taxa with >90% 0s
L2.filtered <- zeroFilter(L2.filtered, threshold = 0.9)

# Remove unassigned or contains environmental
to_remove = c("root;", "root;Not assigned;")
L2.filtered <- as.data.frame(L2.filtered[!(row.names(L2.filtered) %in% to_remove),])
L2.filtered <- L2.filtered |> filter(!str_detect(rownames(L2.filtered), "environment"))

L2.filtered <- as.data.frame(t(L2.filtered))

# Calculate firmicutes/bacteroidetes phyla ratio
firm_bac_ratio <- as.data.frame(L2.filtered$`root;cellular organisms;Bacteria;Terrabacteria group;Firmicutes;`/L2.filtered$`root;cellular organisms;Bacteria;FCB group;Bacteroidetes/Chlorobi group;Bacteroidetes;`)
colnames(firm_bac_ratio) <- "firm_bac_ratio"
row.names(firm_bac_ratio) <- row.names(L2.filtered)

rm(L2, L2.filtered, to_remove, mapping, zeroFilter)

# ----
# Calculate Alpha diversity

# Load shannon diversity calculations

path <- paste0(config$data_processed, "alpha-diversity/")
file_list <- list.files(path)

df <- list()

for (file in file_list){
  load(paste0(path,file))
  df <- cbind(df, shannon_diversity)
}

# Calculate the average Shannon diversity
df <- as.data.frame(df) 
row_names_shannon <- row.names(df)
df <- df %>%
  unnest(cols = colnames(df))
average_shannon <- rowMeans(df)
average_shannon <- as.data.frame(average_shannon)
colnames(average_shannon) <- "shannon"
row.names(average_shannon) <- row_names_shannon

# Calculate average richness
average_richness <- rarefy(L7, sample = rarefaction_depth, se = F)
average_richness <- as.data.frame(average_richness)
attr(average_richness$average_richness, "Subsample") <- NULL

# ----
# Merge all metagenomics data together

average_richness <- average_richness %>% tibble::rownames_to_column(var = "rowname")
average_shannon <- average_shannon %>% tibble::rownames_to_column(var = "rowname")
firm_bac_ratio <- firm_bac_ratio %>% tibble::rownames_to_column(var = "rowname")

metagenomics <- firm_bac_ratio %>%
  full_join(average_richness, by = "rowname") %>%
  full_join(average_shannon, by = "rowname")

metagenomics <- metagenomics %>% tibble::column_to_rownames(var = "rowname")

write.csv(metagenomics, paste0(config$data_processed, "comfort_metagenomics_diversityAndFBratio.csv"), row.names = TRUE)

