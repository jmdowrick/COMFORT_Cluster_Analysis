# COMFORT Cluster Analysis

This repository contains the code used for the cluster analysis of the COMFORT cohort (see [Heenan et al., 2020](https://doi.org/10.1159/000508160) for cohort details). The scripts are organised into categories: preprocessing, mplus factor analysis input files,  analysis, and figure generation. 

## preprocessing
Preprocessing was performed by running `preprocess_data.py`, which:
- preprocesses biological data,
- prepares patient reported outcome data for mplus factor analysis,
- performs rarefaction for alpha diversity calculations,
- extracts metagenomic diversity metrics and computes the F-B ratio,
- and creates the data frame underlying Figure 1.

## mplus_factor-analysis
The mplus input files are stored here. `240808_efa_comfort.inp` calculates the factor loading and `240808_esem_comfort_factor_scores.inp` calculates the factor scores.

## run_analysis
- The suitability of patient reported outcome data for factor analysis was performed using `pre-exploratory_factor_analysis_tests.r`. 
- The optimal number of profiles for latent profile analysis was determined in `num_symptom_profiles.r`
- Age, gender, and diagnostic composition of the cohort was performed in `cohort_characterisation.r`.
- Latent profile analysis was performed using `symptom_clustering.r`.
- Biological cluster analysis was performed using `biological_clustering.r`.
- Integrated cluster analysis was performed using `integrated_cluster.py`.
- Scripts for generating the dataframes underlying Figures 2-6 and S1-3 are named to match the corresponding figures. 

## make_figures
Panels for all figures in the paper were generated using `make_figures.py`.

## References
Heenan P, Creemers RH, Sharma S, Keenan J, Bayer S, Young W, Cooney J, Armstrong K, Fraser K, Skidmore PM, Talley NJ, Roy N, Gearry RB. *Cohort Profile: The Christchurch IBS cOhort to investigate Mechanisms FOr gut Relief and improved Transit (COMFORT).* Inflamm Intest Dis. 2020;5(3):132–143. doi:[10.1159/000508160](https://doi.org/10.1159/000508160).
