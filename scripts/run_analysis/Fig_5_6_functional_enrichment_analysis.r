library(MicrobiomeProfiler)
library(dplyr)
library(tidyverse)
library(readr)
config <- config::get()

load(paste0(config$data_results, "FEA_bg_genes.RData"))

# Cluster 1 vs Cluster 6
source("scripts/run_analysis/utils/FEA_1-6.r")

# Cluster 11 vs Cluster 6
source("scripts/run_analysis/utils/FEA_11-6.r")
