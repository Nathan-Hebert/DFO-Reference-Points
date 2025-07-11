#
# Setup one of the simulated stocks and compute reference points
#
################################################################################

# Load necessary libraries
library(ggplot2)
source('0_ref_pt_helper_functions_231221.R')

# Set seed
set.seed(666)

# Read in necessary inputs
D <- read.csv(paste0(getwd(), "/Data/Pons_inputs.csv"))
dat <- D[, stock]; names(dat) <- D$Parameter
assign(names(dat)[1], dat[1])

# Setup age information... age zero recruits
max_age <- dat['amax']
ages <- 0:max_age
nages <- length(ages)

# Set other important values
R0 <- 1 # Initial recruitment (scales the NAA) in millions
M <- rep(dat['M'], nages)
h <- dat['h']

# Calculate LAA and WAA
laa <- dat['Linf']*(1-exp(-dat['VBk']*(ages-dat['t0']))) # Units are cm
waa <- dat['LWa']*laa^dat['LWb'] # Units are kg

# Calculate selectivity and maturity curves based on laa
lmat_b1 <- log(9.5)/(dat['L95']-dat['L50'])
lmat_b0 <- log(9.5) - lmat_b1*dat['L95']
lmat <- function(x){1/(1+exp(-lmat_b0-lmat_b1*x))}
mat <- lmat(laa)
lsel_b1 <- log(9.5)/(dat['SL95']-dat['SL50'])
lsel_b0 <- log(9.5) - lsel_b1*dat['SL95']
lsel <- function(x){1/(1+exp(-lsel_b0-lsel_b1*x))}
sel <- lsel(laa)

# Fit the glms used to shift the selectivity and maturity curves by age... 
# then use their predictions for no shift as the selectivity and maturity vectors 
# from here on out
sel_mat_data <- as.data.frame(cbind(sel, mat, age = 0:max_age))
sel_fit <- glm(sel ~ age, data = sel_mat_data, family = "binomial")
mat_fit <- glm(mat ~ age, data = sel_mat_data, family = "binomial")
sel_age <- plogis(predict(sel_fit))
mat_age <- plogis(predict(mat_fit))

# Calculate base phi0, and the base alpha and beta SSR parameters
phi0 <- sum(survivorship(nages, M, sel_age, f=0)*waa*mat_age)
BHa <- 4*h/((1-h)*phi0)
BHb <- 1/R0*(BHa-1/phi0)

# Set up dataframe holding a grid of rec dev multipliers, M, mat, sel, waa
M_factors <- seq(0.6, 1.4, by = .2)
waa_factors <- seq(0.6, 1.4, by = .2)
sel_mat_shift <- seq(-1, 1, by = .5) # Left-shifted = -, right-shifted = +
recdev_factors <- seq(0.6, 1.4, by = .2)
DF <- expand.grid(stock = stock,
                  M_factor = M_factors,
                  recdev_factor = recdev_factors,
                  waa_factor = waa_factors,
                  sel_shift = sel_mat_shift,
                  mat_shift = sel_mat_shift,
                  alpha = BHa,
                  beta = BHb)
DF$M <- DF$M_factor*M[1]
# Adjust alpha and beta based on rec dev multipliers
for(i in 1:nrow(DF)){DF[i, c("alpha", "beta")] <- ab_adjust(DF$alpha[i],
                                                        DF$beta[i], 
                                                        DF$recdev_factor[i])}

# Loop through the grid and calculate RPs
for(i in 1:nrow(DF)){
  
  # Shift the selectivity and maturity curves as indicated
  mat_age <- 0:max_age-DF$mat_shift[i]
  sel_age <- 0:max_age-DF$sel_shift[i]
  maturity <- plogis(predict(mat_fit, 
                                newdata = data.frame(age = mat_age)))
  selectivity <- plogis(predict(sel_fit, 
                                newdata = data.frame(age = sel_age)))
  
  # Get new R0 from phi0, alpha, beta
  phi0 <- sum(survivorship(nages, DF$M[i], selectivity, f=0)*waa*DF$waa_factor[i]
              *maturity)
  R0 <- (DF$alpha[i]*phi0-1)/(DF$beta[i]*phi0)
  
  # Estimate RPs
  RP_DF <- as.data.frame(equilibrium_RPs(n_ages = nages, 
                             m = DF$M[i], 
                             selectivity = selectivity,
                             waa = waa*DF$waa_factor[i], mat = maturity,
                             srr = "beverton-holt", alpha = DF$alpha[i], beta = DF$beta[i],
                             SPR = 0.4))
  RP_DF <- RP_DF[c("phi_0", "SSB_0", "F_MSY", # Grab certain values
                   "SSB_MSY", "F_X.SPR", "F_MSY_SPR", "MSY")]
  DF[i, names(RP_DF)] <- as.vector(unlist(RP_DF[1,]))
  DF[i, "R0"] <- R0 # Store R0 along with the RPs
}

# Save DF so can visualize alongside those from the other stocks
saveRDS(DF, paste0(getwd(), "/simulated stock RPs/DF_stock", stock, ".rds"))