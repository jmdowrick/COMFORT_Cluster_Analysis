import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist
from utils.helper_functions import *

config_data = load_config('config.yml')

df_merged = pd.read_csv(config_data['default']['data_consensus_matrices']+"merged_consensus_matrix.csv", index_col=0)
df_merged.index = df_merged.index.astype('str')

# Number of clusters identified by running: scripts/make_figures/make_figures.py
model = AgglomerativeClustering(metric='precomputed', n_clusters = 11, linkage='complete').fit(1-df_merged)

unique_labels = np.unique(model.labels_)
centroids = np.array([df_merged[model.labels_ == i].mean(axis=0) for i in unique_labels])

linkage_matrix = linkage(pdist(centroids), method='complete')

leaf_order = leaves_list(linkage_matrix)
new_labels_map = {old_label: new_label+1 for new_label, old_label in enumerate(leaf_order)}

model.labels_ = np.array([new_labels_map[label] for label in model.labels_])
centroids_ordered = centroids[list(new_labels_map.keys())]
centroid_labels = [f"Cluster {i+1}" for i in list(new_labels_map.keys())]
linkage_matrix = linkage(pdist(centroids_ordered), method='complete')

df_linkage = pd.DataFrame(linkage_matrix)
df_linkage.to_csv(config_data['default']['data_results'] + 'Fig_4A_dendrogram.csv', header=True)

ref_clust_merged = pd.DataFrame(
    {'ref_clust': model.labels_},
    index= df_merged.index
)

ref_clust_merged.to_csv(config_data['default']['data_results']+'integrated_clusters.csv', header=True)