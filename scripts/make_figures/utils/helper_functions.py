# Helper functions
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import yaml
from matplotlib import font_manager
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import davies_bouldin_score, silhouette_score
from scipy.cluster.hierarchy import dendrogram

def load_font(config_data):
    font_path = config_data['default']['font_path']
    source_sans = font_manager.FontProperties(fname=font_path)
    font_manager.fontManager.addfont(font_path)
    plt.rcParams["font.family"] = source_sans.get_name()
    plt.rcParams.update({'font.size': 14})
    return

def make_figure_1A(config_data):
    df_present = pd.read_csv(config_data['default']['data_results'] + "Fig_1_data_availability.csv", index_col=0)
    # Make plot
    idx = ['Fecal bile acids', 'Fecal organic acids', 'Plasma amino acids', 'Plasma metabolites', 'Fecal metabolites',
       'Fecal taxonomic abundance', 'Fecal gene abundance', 'SAGIS', 'PROMIS', 'HADS']

    cm = ['Blues', 'Blues', 'Blues', 'Reds', 'Reds', 'Purples', 'Purples', 'Oranges', 'Oranges', 'Oranges']
    
    f, axs = plt.subplots(10, 1, gridspec_kw={'hspace': 0})

    counter = 0
    for index, row in df_present.iterrows():
        sns.heatmap(np.array([row.values]), yticklabels=[idx[counter]], xticklabels=[], ax=axs[counter], cmap=cm[counter], cbar=False, vmax=2,vmin=0)
        counter += 1

    for ax in axs:
        ax.tick_params(axis='y', rotation=0)

    plt.savefig(config_data['default']['figures']+"Fig_1A_datacollection_heatmap.svg")
    plt.clf()
    return

def make_figure_2B(config_data):
    df_heatmap = pd.read_csv(config_data['default']['data_results']+"Fig_2B_LPA_heatmap.csv", index_col = 0)
    ax = sns.heatmap(df_heatmap, annot=True, fmt=".1f", linewidth=.5) 
    ax.set_facecolor('grey')
    plt.savefig(config_data['default']['figures']+"Fig_2B_symptom_heatmap.svg")
    plt.clf()
    return

def make_figure_4A(config_data):
    df_linkage = pd.read_csv(config_data['default']['data_results'] + 'Fig_4A_dendrogram.csv', index_col=0)
    linkage_matrix = df_linkage.values
    dendrogram(linkage_matrix, labels = [str(i+1) for i in range(11)], orientation = 'top', color_threshold=0, above_threshold_color='k')
    plt.axis('off')
    fig = plt.gcf()
    fig.set_size_inches(8,2)
    plt.savefig(config_data['default']['figures']+"Fig_4A_integrated_dendrogram.svg")
    plt.clf()
    return

def make_figure_S3A_S3B(config_data):
    df_merged = pd.read_csv(config_data['default']['data_consensus_matrices']+"merged_consensus_matrix.csv", index_col=0)

    # Figure S3A - Davies Bouldin Index
    centers = list(range(2, 15)) 
    avg_scores = []
    for center in centers:
        avg_scores.append(get_Hmeans_score((1-df_merged), 'precomputed', 'complete', center))

    fig, ax = plt.subplots(
        figsize=(6, 5)
    )

    plt.rcParams['svg.fonttype'] = 'none'
    plt.axvline(x = 11, color = "lightgray", label = 'axvline - full height', linestyle="--")
    plt.plot(centers, avg_scores, linestyle="-" , marker="o", color="#8D99AE", linewidth =2)
    plt.plot(centers[9], avg_scores[9], marker="o", color="#2B2F42")
    plt.xlabel('Number of clusters')
    plt.ylabel('Davies-Bouldin Index')
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.xticks(np.arange(2, 15, 1.0))
    plt.savefig(config_data['default']['sup_figures']+"Fig_S3A_DBI_merged.svg")
    plt.clf()
    #print(avg_scores[9])

    # Figure S3B - Silhouette Score
    silhouette_scores = []
    for n_clusters in range(2, 15):
        cluster_labels = AgglomerativeClustering(metric = 'precomputed', linkage='complete',  n_clusters=n_clusters).fit_predict(1-df_merged)
        silhouette_scores.append(silhouette_score(1-df_merged, cluster_labels))

    fig, ax = plt.subplots(
        figsize=(6, 5)
    )  

    plt.axvline(x = 11, color = "lightgray", label = 'axvline - full height', linestyle="--")
    plt.plot(range(2, 15), silhouette_scores, linestyle="-" , marker="o", color="#8D99AE", linewidth=2)
    plt.plot(11, silhouette_scores[9], marker="o", color="#2B2F42")
    plt.xlabel('Number of clusters')
    plt.ylabel('Silhouette score')
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.xticks(np.arange(2, 15, 1.0))
    plt.savefig(config_data['default']['sup_figures']+"Fig_S3B_silhouette_merged.svg")
    plt.clf()
    #print(silhouette_scores[9])

    return

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

def load_config(config_file_path):
    try:
        with open(config_file_path, 'r') as stream:
            config = yaml.safe_load(stream)
        return config
    except yaml.YAMLError as exc:
        print(f"Error loading YAML file: {exc}")
        return None