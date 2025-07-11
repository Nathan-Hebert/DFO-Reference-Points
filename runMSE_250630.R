#
# Conduct closed-loop simulations for a given OM, and save the results
#
################################################################################

# Clear workspace
rm(list=ls())

# Load necessary libraries
library(MSEtool); library(SAMtool); library(tidyr)
source("MPs_250626.R")
source("R:/Science/Population Ecology Division/Shared/REFPT/REFPT_R_Project/0_ref_pt_helper_functions_231221.R")

# Choose the OM of interest and number of sims (<= OM@nsim)
OM_file_name <- "M_continued_decline.rds"
n_sims <- 100

###########################MSE Setup + Simulations##############################

# Grab the OM
OM <- readRDS(paste0("OM/",OM_file_name))

# Set important quantities
OM_name <- sub("\\.(rds)$", "", OM_file_name) # Determine OM name from file name
OM@nsim <- n_sims # Number of simulations
MPs <- c("MP_last5", "MP_first5", "MP_entirets", "MP_mixed", "MP_mixed_BPA", "NFref")
MP_num <- length(MPs)
init_yr <- OM@CurrentYr - (OM@nyears-1) # Initial year
term_yr <- OM@CurrentYr # Terminal year

# Set up storage for EMs, control points, and other quantities needed later
empty_matrix <- matrix(NA, nrow = OM@proyears, ncol = OM@nsim)
empty_list_matrix <- matrix(vector("list", OM@proyears/5*OM@nsim), 
                            nrow = OM@proyears/5, ncol = OM@nsim)
CP_elements <- c("F_MSY_First5", "F_MSY_Last5", "F_MSY_EntireTS", 
              "F_MSY_Last5Mixed","F_MSY_Last5_MixedBPA", 
              "SSB_MSY_First5", "SSB_MSY_Last5", "SSB_MSY_EntireTS", 
              "SSB_MSY_First5Mixed", "SSB_MSY_Last5Mixed", 
              "SSB_MSY_First5_MixedBPA",
              "q_First5", "q_Last5", "q_EntireTS", 
              "q_Mixed","q_MixedBPA",
              "SSB_Sigma_MixedBPA")
EM_elements <- c("First5","Last5","EntireTS","Mixed","MixedBPA")
CP_list <- setNames(replicate(length(CP_elements), empty_matrix, 
                              simplify = FALSE), CP_elements)
EM_list <- setNames(replicate(length(EM_elements), empty_list_matrix, 
                             simplify = FALSE), EM_elements)

# Initialize parallel processing and run MSE on the MPs
setup()
MSE_results <- lapply(seq_along(MPs), function(i) {
  runMSE(OM, MPs = MPs[i], silent = TRUE, extended = TRUE, checkMPs = FALSE)
})

# Adjust units within CP_list as required
for (name in c("SSB_MSY_First5", "SSB_MSY_Last5", "SSB_MSY_EntireTS", 
               "SSB_MSY_First5Mixed", "SSB_MSY_Last5Mixed",
               "SSB_MSY_First5_MixedBPA", "SSB_MSY_Last5_MixedBPA")) {
  if (name %in% names(CP_list)) {
    CP_list[[name]] <- CP_list[[name]]/1000
  }
}

# Store MSE results in DF
DF <- data.frame(SSB=unlist(lapply(MSE_results, function(mse) {as.vector(mse@SSB)}))/1000,
                 C=unlist(lapply(MSE_results, function(mse) {as.vector(mse@Catch)}))/1000,
                 FM = unlist(lapply(MSE_results, function(mse) {as.vector(mse@FM)})),
                 Rec = unlist(lapply(MSE_results, function(mse) {2*as.vector(mse@N[,1,,,1])}))/1000,
                 Sim=rep(rep(1:OM@nsim, OM@proyears), MP_num),
                 MP=rep(MPs, each = OM@proyears*OM@nsim),
                 Yr=rep(rep(1:OM@proyears, each=OM@nsim, MP_num)))
DF <- sim_LRP_SSB0(OM, MSE_results[[1]], DF, srr = "beverton-holt") # Append, to each row, the correct 
                                                                    # sim's LRPs, USRs, and first 5 years SSB0
for (col in c("First5", "Last5", "EntireTS")) {
  for (suffix in c("LRP", "USR")) {
    DF[[paste0("I_Above_", col, suffix)]] <- as.integer(DF$SSB > DF[[paste0(suffix, "_", col)]])
  }
}

# Grab historical sim data
hist <- data.frame(MSE_results[[1]]@SSB_hist)/1000
colnames(hist) <- init_yr:term_yr
hist <- gather(hist, key = "Year", value = "SSB")
hist$CB <- gather(data.frame(MSE_results[[1]]@CB_hist), value = "CB")$CB/1000
hist$FM <- gather(data.frame(MSE_results[[1]]@FM_hist), value = "FM")$FM
hist$Year <- as.numeric(hist$Year)
hist$Sim <- rep(1:OM@nsim, OM@nyears)

# Expand historical sim data to bring in the MPs
hist <- do.call(rbind, replicate(MP_num, hist, simplify = FALSE))
hist$MP <- rep(MPs, each = OM@nsim*OM@nyears)

# Save data
dir.create(file.path(paste0(getwd(),"/Figs/MSE/",OM_name)), showWarnings = FALSE)
saveRDS(DF,paste0(getwd(),"/Figs/MSE/",OM_name,"/DF.rda"))
saveRDS(hist, paste0(getwd(),"/Figs/MSE/",OM_name,"/hist.rda"))
saveRDS(CP_list, paste0(getwd(),"/Figs/MSE/",OM_name,"/CP_list.rda"))
saveRDS(EM_list, paste0(getwd(),"/Figs/MSE/",OM_name,"/EM_list.rda"))

# Generate figures
source("generate_MSE_figs_250630.R")