library(readxl)
library(tidyverse)

config <- config::get()

df_attributes <- readxl::read_excel(paste0(config$data_raw, config$file_PRO)) %>%
  select(c(Diagnosis_CaseControls, Gender_numerical, AGE, CASE_CONTROL))

# Age ----  
age_characteristics <- group_by(df_attributes, Diagnosis_CaseControls) %>%
  select(c(Diagnosis_CaseControls, AGE)) %>%
  drop_na() %>%
  summarise(avg = mean(AGE), std = sd(AGE))

age_differences <- select(df_attributes, c(Diagnosis_CaseControls, AGE))

control_age <- age_differences[age_differences$Diagnosis_CaseControls == 6, ]
case_age <- age_differences[age_differences$Diagnosis_CaseControls == 1, ]

age_ttest <- t.test(control_age$AGE, case_age$AGE)
age_ttest

# Gender ----
gender_characteristics <- group_by(df_attributes, Diagnosis_CaseControls) %>%
  select(c(Diagnosis_CaseControls, Gender_numerical)) %>%
  table() %>%
  prop.table(margin = 1)

gender_characteristics

gender_characteristics <- group_by(df_attributes, Diagnosis_CaseControls) %>%
  select(c(Diagnosis_CaseControls, Gender_numerical)) %>%
  table() 

chisq.test(gender_characteristics)

# Diagnosis breakdown ----
diagnosis_composition <- df_attributes %>%
  select(CASE_CONTROL) %>%
  table() 

diagnosis_composition
