from utils.preprocessing_helper import *
import os
import subprocess

config_data = load_config('config.yml')

# Preprocess biological data
subprocess.call("Rscript analysis/scripts/preprocessing/utils/preprocess_biological_data.r", stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)

# Prepare patient reported outcome data for mplus factor analysis
df_mplus_PROs = prep_PROs_for_mplus(config_data)

df_mplus_PROs.to_csv(config_data['default']['data_mplus']+"comfort_PROs.csv", index=False)
with open(os.path.normpath(config_data['default']['data_mplus']+"item_names.txt"), 'w') as fp:
    fp.write(' '.join(df_mplus_PROs.columns.tolist()))

# Rarefy data for alpha diversity calculations
file_path = "./data/processed/alpha-diversity"
if os.path.isdir(file_path):
    for dirpath, dirnames, files in os.walk(file_path):
        if files:
            print('rarefaction already performed')
        if not files: 
            subprocess.call("Rscript analysis/scripts/preprocessing/utils/metagenomics_alpha_diversity.r", stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)
else:
    os.mkdir(file_path)
    subprocess.call("Rscript analysis/scripts/preprocessing/utils/metagenomics_alpha_diversity.r", stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)

# Extract diversity metrics and F-B ratio from metagenomic data
subprocess.call("Rscript analysis/scripts/preprocessing/utils/Fig_4_metagenomics_additional_analysis.r", stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)

# Figure 1 data availability (uses outputs from preprocess biological data)
df_present = fig_1_preprocessing(config_data)
df_present.to_csv(config_data['default']['data_results']+"Fig_1_data_availability.csv")