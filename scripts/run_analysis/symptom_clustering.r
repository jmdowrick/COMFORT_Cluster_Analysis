library(janitor)
library(plyr)
library(tidyr)
library(cowplot)
library(ggplot2)
library(dqrng)

set.seed(123)

config <- config::get()

factor_scores <- read.table(paste0(config$data_results, "comfort_factor_scores.dat"), header=FALSE) 

factor_scores <- factor_scores %>%
  dplyr::select('V43', 'V45', 'V47', 'V49', 'V51', 'V53', 'V55') %>%
  dplyr::rename(., c("Depression" = "V43", "Anxiety" = "V45", 
                     "Diarrhea" = "V47", "Constipation" = "V49", 
                     "Pain/bloating" = "V51", "Upper GI" = "V53", 
                     "Nausea and vomiting" = "V55"))

# Perform latent profile analysis ----
output <- factor_scores %>%
  tidyLPA::estimate_profiles(7,
                    variances = "equal",
                    covariances = "zero") 

ids <- read.csv(paste0(config$data_mplus, "comfort_PROs.csv"))

clust_labels_symptoms <- data.frame(ids$id, output$model_1_class_7$model$classification)
colnames(clust_labels_symptoms) <- c('id', 'cluster')

write.csv(clust_labels_symptoms, paste0(config$data_results, "symptom_clusters.csv"), row.names = FALSE)

# Calculate consensus matrix ----

# Transpose factor_scores and make the subject id the column name
row.names(factor_scores) <- clust_labels_symptoms$id 
factor_scores <- as.data.frame(t(factor_scores))

# Compute consensus matrix
num_resample <- 1000 # H
subjects <- clust_labels_symptoms$id

N <- length(subjects)

# original dataset has N subjects - dim(part)
resample_proportion <- 0.75

# create set of empty connectivity matrices (M and I)
I_k <- as.data.frame(matrix(0, N, N))
M_k <- as.data.frame(matrix(0, N, N))

colnames(I_k) <- subjects
colnames(M_k) <- subjects

rownames(I_k) <- subjects
rownames(M_k) <- subjects

for (h in 1:num_resample) {
  # bootstrap from the data (make D(h))
  D_h <- dqrng::dqsample(subjects, size = round(resample_proportion * N), replace = FALSE)

  # extract data associated with the bootstrapping
  factor_scores_h <- subset(factor_scores, select = as.character(D_h[as.character(D_h) %in% colnames(factor_scores)]))

  # get cluster assignments 
  temp_clust <- as.data.frame(t(factor_scores_h)) %>%
    tidyLPA::estimate_profiles(7,
                    variances = "equal",
                    covariances = "zero") 

  # create assignment with id as col names
  clustering <- temp_clust$model_1_class_7$model$classification
  
  same_cluster <- outer(clustering, clustering, "==")
  # compute M(h) - N x N matrix of whether subjects were in the same cluster
  M_k[names(clustering), names(clustering)] <- M_k[names(clustering), names(clustering)] + same_cluster

  # create I(h) - add 1 to participants present in this bootstrapping
  I_k[names(clustering), names(clustering)] <- I_k[names(clustering), names(clustering)] + 1
}
# M(K) <- compute consensus matrix from M - {M(1), ..., M(H)}
consensus_list_symptoms <- M_k / I_k

write.csv(consensus_list_symptoms, paste0(config$data_results, "symptom_consensus_matrix.csv"), row.names = TRUE)

# Export heat_map data ----
heat_map <- output$model_1_class_7$estimates

heat_map <- heat_map %>%
  dplyr::filter(Category == "Means") %>%
  dplyr::select(c(Parameter, Estimate, Class)) %>%
  tidyr::pivot_wider(names_from = Parameter, values_from = Estimate) %>%
  dplyr::select(-Class) %>%
  t() %>%
  as.data.frame()

colnames(heat_map) <- c("1", "2", "3", "4", "5", "6", "7")
row.names(heat_map) <- c("Depression", "Anxiety", "Diarrhea", "Constipation", "Pain/Bloating", "Upper GI", "Nausea/Vomiting")
write.csv(heat_map, paste0(config$data_results, "Fig_2B_LPA_heatmap.csv"), row.names = TRUE)
