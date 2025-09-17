require(EFA.dimensions)
require(data.table)
library(psych)
library(polycor)

config <- config::get()

df <- utils::read.csv(paste0(config$data_mplus, "comfort_PROs.csv"))

df[df == -999] <- NaN

# dropped items that were not loaded onto any items in mplus
df <- subset(df, select = -c(id, sagis_15, sagis_22))

cols <- c("E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10", "E11", 
          "E12", "E13", "E14", "sagis_1", "sagis_2", "sagis_3", "sagis_4",
          "sagis_5", "sagis_6", "sagis_7", "sagis_8", "sagis_9", "sagis_10",
          "sagis_11", "sagis_12", "sagis_13", "sagis_14", "sagis_16",
          "sagis_17", "sagis_18", "sagis_19", "sagis_20", "sagis_21")
df[cols] <- df[cols] + 1

# Screen for factor analysis data suitability
R <- polycor::hetcor(df, ML = "TRUE", parallel=TRUE)

kmo_out <- KMO(R$correlations)

bartlett_out <- cortest.bartlett(R$correlations, 315)

MAP(R$correlations, Ncases = 315, verbose=TRUE)

SCREE_PLOT(R$correlations, Ncases = 315, verbose = TRUE)
