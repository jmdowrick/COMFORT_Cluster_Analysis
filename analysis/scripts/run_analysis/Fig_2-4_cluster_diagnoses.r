# This script restructures datasets for diagnosis breakdown figures
library(tidyr)
library(janitor)

diagnosis_restructure <- function(df, num_clust) {
    
    IBS_D <- table(factor(df[[2]][[1]]$cluster[df[[2]][[1]]$symptom > 0],levels=1:(num_clust)))
    IBS_D_low <- table(factor(df[[2]][[1]]$cluster[df[[2]][[1]]$symptom < 1],levels=1:(num_clust)))

    Control_symptoms <- table(factor(df[[2]][[2]]$cluster[df[[2]][[2]]$symptom > 0],levels=1:(num_clust)))
    Control <- table(factor(df[[2]][[2]]$cluster[df[[2]][[2]]$symptom < 1],levels=1:(num_clust)))

    IBS_C <- table(factor(df[[2]][[3]]$cluster[df[[2]][[3]]$symptom > 0],levels=1:(num_clust)))
    IBS_C_low <- table(factor(df[[2]][[3]]$cluster[df[[2]][[3]]$symptom < 1],levels=1:(num_clust)))

    IBS_M <- table(factor(df[[2]][[4]]$cluster[df[[2]][[4]]$symptom > 0],levels=1:(num_clust)))
    IBS_M_low <- table(factor(df[[2]][[4]]$cluster[df[[2]][[4]]$symptom < 1],levels=1:(num_clust)))

    FD <- table(factor(df[[2]][[5]]$cluster[df[[2]][[5]]$symptom > 0],levels=1:(num_clust)))
    FD_low <- table(factor(df[[2]][[5]]$cluster[df[[2]][[5]]$symptom < 1],levels=1:(num_clust)))

    FC <- table(factor(df[[2]][[6]]$cluster[df[[2]][[6]]$symptom > 0],levels=1:(num_clust)))
    FC_low <- table(factor(df[[2]][[6]]$cluster[df[[2]][[6]]$symptom < 1],levels=1:(num_clust)))

    df_cont_symp <- data.frame(Control_symptoms)
    df_cont <- data.frame(Control)

    df_ibsd <- data.frame(IBS_D)
    df_ibsd_low <- data.frame(IBS_D_low)

    df_ibsc <- data.frame(IBS_C)
    df_ibsc_low <- data.frame(IBS_C_low)

    df_ibsm <- data.frame(IBS_M)
    df_ibsm_low <- data.frame(IBS_M_low)

    df_fd <- data.frame(FD)
    df_fd_low <- data.frame(FD_low)

    df_fc <- data.frame(FC)
    df_fc_low <- data.frame(FC_low)

    df_comb <- data.frame(rbind(df_cont_symp, df_cont, 
                            df_ibsd, df_ibsd_low,
                            df_ibsc, df_ibsc_low,
                            df_ibsm, df_ibsm_low,
                            df_fd, df_fd_low,
                            df_fc, df_fc_low),
                    group = rep(c("Symp Cont.","Control", 
                                    "IBS-D", "IBS-D low",
                                    "IBS-C", "IBS-C low",
                                    "IBS-M", "IBS-M low",
                                    "FD", "FD low",
                                    "FC", "FC low"), each = num_clust))
    colnames(df_comb) <- c("Class", "Count (n)", "Diagnosis")

    df_comb$`Cluster (n)` <- with(df_comb, ifelse((Diagnosis == "Control" | Diagnosis == "Symp Cont."), "Healthy", "DGBI"))
    
    return(df_comb)
}

config <- config::get()

df <- readxl::read_excel(paste0(config$data_raw, config$file_PRO), sheet="METmin") %>%
  dplyr::select("Participant ID", "CASE_CONTROL") %>%
  dplyr::rename(id = 'Participant ID') %>%
  clean_names 

df_promis <- readxl::read_excel(paste0(config$data_raw, config$file_PRO),sheet = "PROMIS") %>%
  as.data.frame(.) %>%
  dplyr::select(c("Participant #", 6:ncol(.))) %>%
  dplyr::rename("id" = "Participant #") %>%
  clean_names

df$symptom <- +(apply(df_promis[2:7] > 55, 1, any, na.rm=TRUE)) 
df['id'] <- as.numeric(df$id)

# Symptom clustering ----
clust_labels_symptom <- read.csv(paste0(config$data_cluster_assignments, "symptom_clusters.csv"), header=TRUE) 
num_clust <- max(clust_labels_symptom$cluster)
df_symptom <- df %>% dplyr::inner_join(clust_labels_symptom, by=c("id")) %>% nest(data = -case_control)
comb <- diagnosis_restructure(df_symptom, num_clust)

save(comb,file=paste0(config$data_results, "Fig_2A_symptom_diagnosis_breakdown.rData"))

# Biological clustering ----
clust_labels_biology <- read.csv(paste0(config$data_cluster_assignments, "biological_clusters.csv"), header=TRUE) %>%
  dplyr::rename(c('cluster'= 'ref_clust','id' = 'X'))
num_clust <- max(clust_labels_biology$cluster)
df_biology <- df %>% dplyr::inner_join(clust_labels_biology, by=c("id")) %>% nest(data = -case_control)
comb <- diagnosis_restructure(df_biology, num_clust)

save(comb,file=paste0(config$data_results, "Fig_3A_biology_diagnosis_breakdown.rData"))

# Integrated clustering ----
clust_labels_integrated <- read.csv(paste0(config$data_cluster_assignments, "integrated_clusters.csv"), header=TRUE)  %>%
  dplyr::rename(c('cluster'= 'ref_clust','id' = 'X'))
clust_labels_integrated$cluster <- clust_labels_integrated$cluster
num_clust <- max(clust_labels_integrated$cluster) 
df_integrated <- df %>% 
  dplyr::inner_join(clust_labels_integrated, by=c("id")) %>% 
  nest(data = -case_control)

comb <- diagnosis_restructure(df_integrated, num_clust)

save(comb,file=paste0(config$data_results, "Fig_4A_integrated_diagnosis_breakdown.rData"))
