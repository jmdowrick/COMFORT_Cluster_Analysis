library(ggplot2)
library(cowplot)
library(ggprism)
library(tidyr)
library(tibble)
library(EnhancedVolcano)
library(showtext)
library(ggpubr)
library(readr)
library(tidyverse)
library(rstatix)
library(tidyquant)
library(ggdist)
library(ggthemes)
library(sysfonts)
library(showtext)

font_add_google("Source Sans Pro", regular.wt = 400, bold.wt = 700)
showtext_auto()

config <- config::get()

# Figure 1 ----
source("analysis/scripts/make_figures/utils/make_figure_1B.r")

# Figure 2 ----
source("analysis/scripts/make_figures/utils/make_figure_2A.r")

# Figure 3 ----
source("analysis/scripts/make_figures/utils/make_figure_3A.r")
source("analysis/scripts/make_figures/utils/make_figure_3B.r")

# Figure 4 ----
source("analysis/scripts/make_figures/utils/make_figure_4A.r")
source("analysis/scripts/make_figures/utils/make_figure_4B.r")
source("analysis/scripts/make_figures/utils/make_figure_4C.r")
source("analysis/scripts/make_figures/utils/make_figure_4D.r")
source("analysis/scripts/make_figures/utils/make_figure_4E.r")
source("analysis/scripts/make_figures/utils/make_figure_4F.r")

# Figure 5 ----
source("analysis/scripts/make_figures/utils/make_figure_5.r")

# Figure 6 ----
source("analysis/scripts/make_figures/utils/make_figure_6.r")

# Supplementary figures ----

# Stability analysis
source("analysis/scripts/make_figures/utils/make_figure_S1.r")
source("analysis/scripts/make_figures/utils/make_figure_S2.r")
source("analysis/scripts/make_figures/utils/make_figure_S3C.r")

# Differential analysis - 1 (DGBI-predominant 1) vs 6 (control reference)
source("analysis/scripts/make_figures/utils/make_figure_S4-S8.r")

# Differential analysis - 11 (DGBI-predominant 2) vs 6 (control reference)
source("analysis/scripts/make_figures/utils/make_figure_S9-S13.r")
