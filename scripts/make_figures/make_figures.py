from utils.helper_functions import *
import subprocess

config_data = load_config('config.yml')

load_font(config_data)

# Figure 1
make_figure_1A(config_data)

# Figure 2B
make_figure_2B(config_data)

# Figure 4A (dendrogram)
make_figure_4A(config_data)

# Supplementary Figures
make_figure_S3A_S3B(config_data)

# Make remaining figure panels using R 
subprocess.call("Rscript scripts/make_figures/utils/make_R_figures.r", stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)
