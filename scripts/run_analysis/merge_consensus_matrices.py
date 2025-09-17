import pandas as pd
import numpy as np
from utils.helper_functions import *

config_data = load_config('config.yml')

# Import symptom and biological clustering results
df_consensus_symptom = pd.read_csv(config_data['default']['data_consensus_matrices'] + "symptom_consensus_matrix.csv", index_col = 0)
df_consensus_symptom.index = df_consensus_symptom.index.astype('str')

df_consensus_biology = pd.read_csv(config_data['default']['data_consensus_matrices'] + "biology_consensus_matrix.csv", index_col = 0)
df_consensus_biology.index = df_consensus_biology.index.astype('str')

# Place consensus matrices in complete dataset
df_symptoms = pd.DataFrame(np.nan, index=np.char.mod('%d', np.arange(90001, 90350)), columns=np.char.mod('%d', np.arange(90001, 90350)))
df_symptoms.loc[df_consensus_symptom.index, df_consensus_symptom.columns] = df_consensus_symptom

df_biology = pd.DataFrame(np.nan, index=np.char.mod('%d', np.arange(90001, 90350)), columns=np.char.mod('%d', np.arange(90001, 90350)))
df_biology.loc[df_consensus_biology.index, df_consensus_biology.columns] = df_consensus_biology

# Calculate the average consensus matrix
df_merged = (df_biology.add(df_symptoms)).div(2)

# drop nan columns and rows
df_merged = df_merged.dropna(how='all')
df_merged = df_merged.dropna(axis=1, how='all')

# fill diagonal with 1
np.fill_diagonal(df_merged.values, 1)

# drop nans 
df_merged = df_merged.loc[df_merged.isna().sum() == 0, df_merged.isna().sum() == 0]
df_merged.to_csv(config_data['default']['data_consensus_matrices']+"merged_consensus_matrix.csv")