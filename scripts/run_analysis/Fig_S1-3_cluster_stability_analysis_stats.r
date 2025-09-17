# Significance of stability 
library(ggpubr)
library(readr)
library(tidyverse)
library(rstatix)
library(tidyquant)
library(ggdist)
library(ggthemes)
library(ggprism)

# Symptom ----

config <- config::get()

symp_stability <- read_csv(paste0(config$data_results, "Fig_S1_symptom_stability.csv")) %>%
  select(c(cluster, item_consensus)) 

clust_one <- symp_stability[symp_stability$cluster == 1,] 
ggqqplot(clust_one$item_consensus)

clust_two <- symp_stability[symp_stability$cluster == 2,] 
ggqqplot(clust_two$item_consensus)

clust_three <- symp_stability[symp_stability$cluster == 3,] 
ggqqplot(clust_three$item_consensus)

clust_four <- symp_stability[symp_stability$cluster == 4,] 
ggqqplot(clust_four$item_consensus)

clust_five <- symp_stability[symp_stability$cluster == 5,] 
ggqqplot(clust_five$item_consensus)

clust_six <- symp_stability[symp_stability$cluster == 6,] 
ggqqplot(clust_six$item_consensus)

clust_seven <- symp_stability[symp_stability$cluster == 7,] 
ggqqplot(clust_seven$item_consensus)

# cluster five has a normal distribution of item-wise consensus 

# Test against threshold 
clust_one_stability <- wilcox.test(clust_one$item_consensus, mu = 0.5, alternative = "greater")
clust_one_stability # stable

clust_two_stability <- wilcox.test(clust_two$item_consensus, mu = 0.5, alternative = "greater")
clust_two_stability # not stable

clust_three_stability <- wilcox.test(clust_three$item_consensus, mu = 0.5, alternative = "greater")
clust_three_stability # stable

clust_four_stability <- wilcox.test(clust_four$item_consensus, mu = 0.5, alternative = "greater")
clust_four_stability # stable

clust_five_stability <- t.test(clust_five$item_consensus, mu = 0.5, alternative = "greater")
clust_five_stability # stable

clust_six_stability <- wilcox.test(clust_six$item_consensus, mu = 0.5, alternative = "greater")
clust_six_stability # stable

clust_seven_stability <- wilcox.test(clust_seven$item_consensus, mu = 0.5, alternative = "greater")
clust_seven_stability # stable

# Biological ----
bio_stability <- read_csv(paste0(config$data_results, "Fig_S2_biology_stability.csv")) %>%
  select(c(cluster, item_consensus)) %>%
  group_by(cluster) %>%
  shapiro_test(item_consensus)

bio_stability

# View QQ plots
bio_stability <- read_csv(paste0(config$data_results, "Fig_S2_biology_stability.csv")) %>%
  select(c(cluster, item_consensus)) 

clust_one_bio <- bio_stability[bio_stability$cluster == 1,] 
ggqqplot(clust_one_bio$item_consensus)

clust_two_bio <- bio_stability[bio_stability$cluster == 2,] 
ggqqplot(clust_two_bio$item_consensus)

clust_three_bio <- bio_stability[bio_stability$cluster == 3,] 
ggqqplot(clust_three_bio$item_consensus)

clust_four_bio <- bio_stability[bio_stability$cluster == 4,] 
ggqqplot(clust_four_bio$item_consensus)

# Test against threshold 
clust_one_stab_bio <- wilcox.test(clust_one_bio$item_consensus, mu = 0.5, alternative = "greater")
clust_one_stab_bio # stable

clust_two_stab_bio <- wilcox.test(clust_two_bio$item_consensus, mu = 0.5, alternative = "greater")
clust_two_stab_bio # stable

clust_three_stab_bio <- wilcox.test(clust_three_bio$item_consensus, mu = 0.5, alternative = "greater")
clust_three_stab_bio # stable

clust_four_stab_bio <- wilcox.test(clust_four_bio$item_consensus, mu = 0.5, alternative = "greater")
clust_four_stab_bio # stable

# Integrated ----

# Test normality
stability <- read_csv(paste0(config$data_results, "Fig_S3_integrated_stability.csv")) %>%
  select(c(cluster, item_consensus)) %>%
  group_by(cluster) %>%
  shapiro_test(item_consensus)

stability

# View QQ plots
stability <- read_csv(paste0(config$data_results, "Fig_S3_integrated_stability.csv")) %>%
  select(c(cluster, item_consensus)) 

clust_one <- stability[stability$cluster == 1,] 
ggqqplot(clust_one$item_consensus)

clust_two <- stability[stability$cluster == 2,] 
ggqqplot(clust_two$item_consensus)

clust_three <- stability[stability$cluster == 3,] 
ggqqplot(clust_three$item_consensus)

clust_four <- stability[stability$cluster == 4,] 
ggqqplot(clust_four$item_consensus)

clust_five <- stability[stability$cluster == 5,] 
ggqqplot(clust_five$item_consensus)

clust_six <- stability[stability$cluster == 6,] 
ggqqplot(clust_six$item_consensus)

clust_seven <- stability[stability$cluster == 7,] 
ggqqplot(clust_seven$item_consensus)

clust_eight <- stability[stability$cluster == 8,] 
ggqqplot(clust_eight$item_consensus)

clust_nine <- stability[stability$cluster == 9,] 
ggqqplot(clust_nine$item_consensus)

clust_ten <- stability[stability$cluster == 10,] 
ggqqplot(clust_ten$item_consensus)

clust_eleven <- stability[stability$cluster == 11,] 
ggqqplot(clust_eleven$item_consensus)

# Test against threshold ----
clust_one_stability <- wilcox.test(clust_one$item_consensus, mu = 0.5, alternative = "greater")
clust_one_stability # stable

clust_two_stability <- wilcox.test(clust_two$item_consensus, mu = 0.5, alternative = "greater")
clust_two_stability # stable

clust_three_stability <- wilcox.test(clust_three$item_consensus, mu = 0.5, alternative = "greater")
clust_three_stability # not stable

clust_four_stability <- wilcox.test(clust_four$item_consensus, mu = 0.5, alternative = "greater")
clust_four_stability # not stable

clust_five_stability <- wilcox.test(clust_five$item_consensus, mu = 0.5, alternative = "greater")
clust_five_stability # stable

clust_six_stability <- wilcox.test(clust_six$item_consensus, mu = 0.5, alternative = "greater")
clust_six_stability # stable

clust_seven_stability <- wilcox.test(clust_seven$item_consensus, mu = 0.5, alternative = "greater")
clust_seven_stability # stable

clust_eight_stability <- wilcox.test(clust_eight$item_consensus, mu = 0.5, alternative = "greater")
clust_eight_stability # stable

clust_nine_stability <- wilcox.test(clust_nine$item_consensus, mu = 0.5, alternative = "greater")
clust_nine_stability # stable

clust_ten_stability <- wilcox.test(clust_ten$item_consensus, mu = 0.5, alternative = "greater")
clust_ten_stability # not stable

clust_eleven_stability <- wilcox.test(clust_eleven$item_consensus, mu = 0.5, alternative = "greater")
clust_eleven_stability # stable
