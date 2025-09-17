# Cluster stability analysis processing
#
# This script consolidates the item-wise consensus scores for downstream statistical analysis in the file
# cluster_stability_analysis_stats.R

from utils.helper_functions import item_consensus, load_config
import pandas as pd

config_data = load_config('config.yml')

# Symptom clustering
df_consensus_symptom = pd.read_csv(config_data['default']['data_consensus_matrices']+"symptom_consensus_matrix.csv", index_col=0)
df_consensus_symptom.index = df_consensus_symptom.index.astype('str')
ref_clust_symptom = pd.read_csv(config_data['default']['data_cluster_assignments']+"symptom_clusters.csv", index_col=0)
ref_clust_symptom.index = ref_clust_symptom.index.astype('str')

symptom_stability = item_consensus(df_consensus_symptom, ref_clust_symptom)

symptom_stability.to_csv(config_data['default']['data_results']+"Fig_S1_symptom_stability.csv", header=True)

# Biological clustering
df_consensus_biology = pd.read_csv(config_data['default']['data_consensus_matrices']+"biology_consensus_matrix.csv", index_col=0)
df_consensus_biology.index = df_consensus_biology.index.astype('str')
ref_clust_biology = pd.read_csv(config_data['default']['data_cluster_assignments']+"biological_clusters.csv", index_col=0)
ref_clust_biology.index = ref_clust_biology.index.astype('str')

biology_stability = item_consensus(df_consensus_biology, ref_clust_biology)
biology_stability.to_csv(config_data['default']['data_results']+"Fig_S2_biology_stability.csv", header=True)

# Integrated clustering
df_consensus_integrated = pd.read_csv(config_data['default']['data_consensus_matrices']+"merged_consensus_matrix.csv", index_col=0)
df_consensus_integrated.index = df_consensus_integrated.index.astype('str')
ref_clust_integrated = pd.read_csv(config_data['default']['data_cluster_assignments']+"integrated_clusters.csv", index_col=0)
ref_clust_integrated.index = ref_clust_integrated.index.astype('str')

integrated_stability = item_consensus(df_consensus_integrated, ref_clust_integrated)

integrated_stability.to_csv(config_data['default']['data_results']+"Fig_S3_integrated_stability.csv", header=True)
