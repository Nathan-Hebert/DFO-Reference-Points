#
# Generates two summary tables of results for the MSE sims
#
################################################################################

# Load libraries
library(tidyr); library(dplyr)

# Clear workspace
rm(list=ls())

# Set important values
rounding <- 2 # Number of decimal places to round tables to
nsims <- 100 # Needed for probability calculations

# Load the individual results for the scenarios and generate a combined dataframe
main_directory <- paste0(getwd(), '/Figs/MSE')
scenarios <- list.dirs(main_directory, full.names = FALSE, recursive = FALSE)
DF <- data.frame()
for (scenario in scenarios) {
  file_path <- paste0(main_directory, '/', scenario, '/DF.rda')
  if (file.exists(file_path)) {
    DF_new <- readRDS(file_path)
    DF_new$scenario <- scenario
    DF <- rbind(DF, DF_new)
  }
}

# Use combined dataframe to summarize long-term and short-term catch, 
# SSB/SSBMSYfirst5, SSB/SSBMSYentirets, and SSB/SSBMSYlast5
RES1 <- DF %>%
  filter(Yr <= 5 | (Yr > 25 & Yr <= 50)) %>%
  mutate(Period = case_when(Yr <= 5 ~ "ST", 
                            Yr > 25 & Yr <= 50 ~ "LT"),
         MP = factor(MP, levels = levels(as.factor(DF$MP)))) %>%
  group_by(scenario, Period, Sim, MP) %>%
  summarize(MedianCatch = median(C, na.rm = TRUE),
            MedianSSB_SSBMSYFirst5 = median(0.4*SSB/LRP_First5, na.rm = TRUE),
            MedianSSB_SSBMSYEntireTS = median(0.4*SSB/LRP_EntireTS, na.rm = TRUE),
            MedianSSB_SSBMSYLast5 = median(0.4*SSB/LRP_Last5, na.rm = TRUE)) %>%
  group_by(scenario, Period, MP) %>%
  summarize(catch = median(MedianCatch, na.rm = TRUE),
            SSB_SSBMSYFirst5 = median(MedianSSB_SSBMSYFirst5, na.rm = TRUE),
            SSB_SSBMSYEntireTS = median(MedianSSB_SSBMSYEntireTS, na.rm = TRUE),
            SSB_SSBMSYLast5 = median(MedianSSB_SSBMSYLast5, na.rm = TRUE))

# Use combined dataframe to summarize long-term and short-term P(SSB > first 5 
# years LRP), P(SSB > last 5 years LRP), P(SSB > entire ts LRP), and P(catch = 0 kt)
RES2 <- DF %>% mutate(I_Prob_Catch = if_else(C < 0.01, 1,0)) %>%
  filter(Yr <= 5 | (Yr > 25 & Yr <= 50)) %>%
  mutate(Period = case_when(Yr <= 5 ~ "ST", 
                          Yr > 25 & Yr <= 50 ~ "LT"),
         MP = factor(MP, levels = levels(as.factor(DF$MP)))) %>%
  group_by(scenario, Period, Yr, MP) %>%
  summarize(MedianLRPFirst5 = median(sum(I_Above_First5LRP)/nsims, na.rm = TRUE),
            MedianLRPEntireTS = median(sum(I_Above_EntireTSLRP)/nsims, na.rm = TRUE),
            MedianLRPLast5 = median(sum(I_Above_Last5LRP)/nsims, na.rm = TRUE),
            MedianProbCatch = median(sum(I_Prob_Catch)/nsims, na.rm = TRUE)) %>%
  group_by(scenario, Period, MP) %>%
  summarize(LRPFirst5 = median(MedianLRPFirst5, na.rm = TRUE),
            LRPEntireTS = median(MedianLRPEntireTS, na.rm = TRUE),
            LRPLast5 = median(MedianLRPLast5, na.rm = TRUE),
            ProbCatch = median(MedianProbCatch, na.rm = TRUE))

# Use summary values from above to generate results table... exclude "no-fishing"
RES <- RES1 %>% left_join(RES2) %>% pivot_wider(names_from = Period,
                                                values_from = c(SSB_SSBMSYFirst5, 
                                                                SSB_SSBMSYLast5,
                                                                SSB_SSBMSYEntireTS,
                                                                catch, 
                                                                LRPFirst5,
                                                                LRPLast5,
                                                                LRPEntireTS,
                                                                ProbCatch),
                                                names_sep = "_") %>% filter(MP != "NFref") 

# Calculate short-term and long-term average annual variability in yield (AAVY) 
# for each sim-HCR-scenario combo...
AAVY_df <- data.frame(Scenario = character(), MP = character(), Sim = integer(), 
                     ST_AAVY = numeric(), LT_AAVY = numeric())
Scenarios <- unique(DF$scenario)
MPs <- unique(DF$MP)
for (scenario in Scenarios) {
  for (mp in MPs) {
    for (i in 1:nsims) {
      catch <- DF$C[DF$scenario == scenario & DF$MP == mp & DF$Sim == i]
      AAVY_value_ST <- mean(abs(diff(catch[1:4])))
      AAVY_value_LT <- mean(abs(diff(catch[25:49])))
      AAVY_df <- rbind(AAVY_df, data.frame(Scenario = scenario, MP = mp, Sim = i, 
                                         ST_AAVY = AAVY_value_ST,
                                         LT_AAVY = AAVY_value_LT))
    }
  }
}
# ... and then take the medians across sims, filter out "no fishing", and append
# to the results table
AAVY_df_med <- AAVY_df %>% group_by(Scenario, MP) %>% 
  summarize(AAVY_ST = median(ST_AAVY), AAVY_LT = median(LT_AAVY)) %>% 
  filter(MP != "NFref") 
RES$AAVY_ST <- AAVY_df_med$AAVY_ST
RES$AAVY_LT <- AAVY_df_med$AAVY_LT

# Do some rounding, split the results table into 2, and save the two tables as csvs
RES[] <- lapply(RES, function(x) if(is.numeric(x)) round(x, rounding) else x)

RES[c("scenario","MP","LRPFirst5_ST","LRPFirst5_LT","SSB_SSBMSYFirst5_ST",
              "SSB_SSBMSYFirst5_LT","ProbCatch_ST","ProbCatch_LT", 
              "catch_ST","catch_LT","AAVY_ST","AAVY_LT")] %>% 
  mutate(MP = factor(MP, levels = c("MP_first5", "MP_entirets", "MP_last5",
                                    "MP_mixed","MP_mixed_BPA"))) %>% arrange(scenario, MP) %>%
  write.csv("mse_results_table1.csv", row.names = FALSE)

RES[c("scenario","MP","LRPFirst5_ST","LRPFirst5_LT","LRPEntireTS_ST","LRPEntireTS_LT",
              "LRPLast5_ST","LRPLast5_LT","SSB_SSBMSYFirst5_ST","SSB_SSBMSYFirst5_LT",
              "SSB_SSBMSYEntireTS_ST","SSB_SSBMSYEntireTS_LT","SSB_SSBMSYLast5_ST",
              "SSB_SSBMSYLast5_LT")] %>% 
  mutate(MP = factor(MP, levels = c("MP_first5", "MP_entirets", "MP_last5",
                                    "MP_mixed","MP_mixed_BPA"))) %>% arrange(scenario, MP) %>%
  write.csv("mse_results_table2.csv", row.names = FALSE)