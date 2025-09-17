# Helper functions
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import yaml
import seaborn as sns
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import davies_bouldin_score
from matplotlib.ticker import FuncFormatter


def get_Hmeans_score(data, distance, link, center):  
    '''
    returns the  score regarding Davies Bouldin for points to centers
    INPUT:
        data - the dataset you want to fit Agglomerative to
        distance - the distance for AgglomerativeClustering
        link - the linkage method for AgglomerativeClustering
        center - the number of clusters you want (the k value)
    OUTPUT:
        score - the Davies Bouldin score for the Hierarchical model fit to the data
    '''
    # instantiate Hierarchical Clustering
    hmeans = AgglomerativeClustering(n_clusters=center, metric=distance, linkage=link)
    
    # fit the model to the data using the fit method
    model = hmeans.fit_predict(data)
    
    # Calculate Davies Bouldin score
    score = davies_bouldin_score(data, model)
    return score

def sort_consensus_matrix(consensus_matrix, ref_clust):
    df_consensus = pd.DataFrame(consensus_matrix)
    df_consensus = pd.concat([df_consensus, ref_clust], axis=1)
    df_consensus = df_consensus.T
    df_consensus = pd.concat([df_consensus, ref_clust], axis=1)

    df_consensus.sort_values(by = df_consensus.last_valid_index(), axis = 1, inplace = True) 
    df_consensus.sort_values(by = df_consensus.last_valid_index(), axis = 0, inplace = True) 
    
    return df_consensus

def item_consensus(consensus_matrix, ref_clust):
    if 'cluster' in ref_clust.columns:
        ref_clust = ref_clust.rename(columns={'cluster': 'ref_clust'})
      
    sorted_consensus_matrix = sort_consensus_matrix(consensus_matrix, ref_clust) 
      
    item_cons = pd.DataFrame(np.array([1,1,1]).reshape(1,-1), columns=['cluster', 'item_consensus', 'row'])
    
    for cluster in sorted_consensus_matrix['ref_clust'].unique():
        if ~np.isnan(cluster):
  
            cluster_items = ref_clust.loc[ref_clust['ref_clust'] == cluster].index
            current_cluster = sorted_consensus_matrix.loc[cluster_items, cluster_items]
            current_cluster = current_cluster.dropna(axis=1, how='all').dropna(axis=0, how='all')

            # Nk = number of items in cluster k
            N_k = current_cluster.shape[0]

            for row in current_cluster:
                if N_k > 1:
                    per_item = np.array([cluster, (1/(N_k-1)) * (sum(current_cluster[row]) - 1), row]) # subtract the diagonal
                    item_cons = pd.concat([item_cons, pd.DataFrame(per_item.reshape(1,-1), columns=item_cons.columns)], ignore_index=True)
                else:
                    per_item = np.array([cluster, 0, row])
                    item_cons = pd.concat([item_cons, pd.DataFrame(per_item.reshape(1,-1), columns=item_cons.columns)], ignore_index=True)
                
    item_cons.drop(0, inplace=True)
    item_cons['item_consensus'] = item_cons['item_consensus'].astype(float)

    return item_cons

def load_config(config_file_path):
    try:
        with open(config_file_path, 'r') as stream:
            config = yaml.safe_load(stream)
        return config
    except yaml.YAMLError as exc:
        print(f"Error loading YAML file: {exc}")
        return None
