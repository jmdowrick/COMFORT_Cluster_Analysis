library(SNFtool)
library(NEMO)
library(janitor)
library(cowplot)
library(ggplot2)
library(random)
library(dqrng)
library(data.table)

set.seed(123)

# dropUnsharedOmes drops participants with an insufficient number of datatypes.
# inputs:
# subjects = dataframe of all participant labels
# listofdf = list of dataframes (e.g., bile acid df and gene abundance df)
# threshold = minimum number of shared datatypes required to retain a subject
#
# outputs: (packaged as a list)
# subjects = dataframe of retained participant labels
# listofdf = list of dataframes containing retained participants only
dropUnsharedOmes <- function(subjects, listofdf, threshold = 2) {
  if (!is.data.frame(subjects)) {
    subjects <- as.data.frame(subjects)
    colnames(subjects) <- "Participant"
  }
  available_data <- subjects
  row.names(available_data) <- subjects[[1]]
  for (i in seq_along(listofdf)) {
    available_data[colnames(listofdf[[i]]), i + 1] <- 1
  }
  available_data[is.na(available_data)] <- 0
  numShared <- rowSums(dplyr::select(available_data, -c(Participant))) 
  
  # only retain subjects with enough shared data types
  subjects <- names(numShared[numShared >= threshold])
  for (i in seq_along(listofdf)) {
    listofdf[[i]] <- dplyr::select(as.data.frame(listofdf[[i]]), any_of(subjects))
  }
  
  return(list(subjects, listofdf))
}

config <- config::get()

# Load data ----
comfort_amino <- read.csv(paste0(config$data_processed, "comfort_amino_acids_cluster.csv"), row.names = 1)
comfort_bile <- read.csv(paste0(config$data_processed, "comfort_bile_acids_cluster.csv"), row.names = 1)
comfort_f_metab <- read.csv(paste0(config$data_processed, "comfort_fmetab_cluster.csv"), row.names = 1)
comfort_gene <- read.csv(paste0(config$data_processed, "comfort_gene_cluster.csv"), row.names = 1)
comfort_taxa <- read.csv(paste0(config$data_processed, "comfort_taxa_cluster.csv"), row.names = 1)
comfort_organic <- read.csv(paste0(config$data_processed, "comfort_organic_acids_cluster.csv"), row.names = 1)
comfort_p_metab <- read.csv(paste0(config$data_processed, "comfort_pmetab_cluster.csv"), row.names = 1)

names(comfort_amino) <- substring(names(comfort_amino), 2)
names(comfort_bile) <- substring(names(comfort_bile), 2)
names(comfort_f_metab) <- substring(names(comfort_f_metab), 2)
names(comfort_gene) <- substring(names(comfort_gene), 2)
names(comfort_taxa) <- substring(names(comfort_taxa), 2)
names(comfort_organic) <- substring(names(comfort_organic), 2)
names(comfort_p_metab) <- substring(names(comfort_p_metab), 2)

# Extract participant names
subjects <- read.table(paste0(config$data_processed, "participant_list.tsv"), header = TRUE)
subjects <- as.character(unlist(subjects))

# Ensure that participants are represented in at least two omics. If not, remove.
output <- dropUnsharedOmes(subjects,
  list(comfort_amino, comfort_bile, comfort_f_metab, comfort_gene,
       comfort_organic, comfort_p_metab, comfort_taxa)
)

subjects <- output[[1]]
listofdf <- output[[2]]

df_amino.aligned <- listofdf[[1]]
df_bile.aligned <- listofdf[[2]]
df_f_metab.aligned <- listofdf[[3]]
df_gene.aligned <- listofdf[[4]]
df_organic.aligned <- listofdf[[5]]
df_p_metab.aligned <- listofdf[[6]]
df_taxa.aligned <- listofdf[[7]]

# Create reference clustering ----
affinity.graph <- nemo.affinity.graph(listofdf)

# ask nemo to estimate the number of clusters.
num.clusters <- nemo.num.clusters(affinity.graph)

# clustering is the cluster assignment vector, ordered by the columns of affinity.graph.
ref_clust <- spectralClustering(affinity.graph, num.clusters)
names(ref_clust) <- colnames(affinity.graph)

clust_labels_biology <- as.data.frame(ref_clust) 
clust_labels_biology <- clust_labels_biology[order(rownames(clust_labels_biology)),,drop=FALSE]
write.csv(clust_labels_biology, paste0(config$data_results, "biological_clusters.csv"), row.names = TRUE)

# Compute consensus matrix ----
num_resample <- 1000 # H
N <- length(subjects)

# original dataset has N subjects - dim(part)
resample_proportion <- 0.75

# create set of empty connectivity matrices (M and I)
I_k <- as.data.frame(matrix(0, N, N))
M_k <- as.data.frame(matrix(0, N, N))

subjects <- subjects[order(subjects)]

colnames(I_k) <- subjects
colnames(M_k) <- subjects

rownames(I_k) <- subjects
rownames(M_k) <- subjects

for (h in 1:num_resample) {
  # bootstrap from the data (make D(h))
  D_h <- dqrng::dqsample(subjects, size = round(resample_proportion * N), replace = FALSE)

  # extract data associated with the bootstrapping
  df_amino_h <- subset(df_amino.aligned, select = as.character(D_h[as.character(D_h) %in% colnames(df_amino.aligned)]))
  df_bile_h <- subset(df_bile.aligned, select = as.character(D_h[as.character(D_h) %in% colnames(df_bile.aligned)]))
  df_f_metab_h <- subset(df_f_metab.aligned, select = as.character(D_h[as.character(D_h) %in% colnames(df_f_metab.aligned)]))
  df_gene_h<- subset(df_gene.aligned, select = as.character(D_h[as.character(D_h) %in% colnames(df_gene.aligned)]))
  df_organic_h <- subset(df_organic.aligned, select = as.character(D_h[as.character(D_h) %in% colnames(df_organic.aligned)]))
  df_p_metab_h <- subset(df_p_metab.aligned, select = as.character(D_h[as.character(D_h) %in% colnames(df_p_metab.aligned)]))
  df_taxa_h <- subset(df_taxa.aligned, select = as.character(D_h[as.character(D_h) %in% colnames(df_taxa.aligned)]))

  # drop subjects if no shared ome
  output <- dropUnsharedOmes(D_h, list(df_amino_h, df_bile_h, df_f_metab_h, 
                                       df_gene_h, df_organic_h,  df_p_metab_h, 
                                       df_taxa_h))
  D_h <- output[[1]]
  omics.list <- output[[2]]

  # cluster D(h) into k clusters
  affinity.graph <- nemo.affinity.graph(omics.list)
  clustering <- spectralClustering(affinity.graph, num.clusters)
  names(clustering) <- colnames(affinity.graph)

  # compute M(h) - N x N matrix of whether subjects were in the same cluster
  same_cluster <- outer(clustering, clustering, "==")
  M_k[colnames(affinity.graph), colnames(affinity.graph)] <- M_k[colnames(affinity.graph), colnames(affinity.graph)] + same_cluster
  
  # create I(h) - add 1 to participants present in this bootstrapping
  I_k[colnames(affinity.graph), colnames(affinity.graph)] <- I_k[colnames(affinity.graph), colnames(affinity.graph)] + 1
}

# M(K) <- compute consensus matrix from M - {M(1), ..., M(H)}
consensus_list_biology <- M_k / I_k

write.csv(consensus_list_biology, paste0(config$data_results, "biology_consensus_matrix.csv"), row.names = TRUE)
