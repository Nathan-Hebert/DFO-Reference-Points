#
# Setups an OM for a simulated stock
#
################################################################################

# Clear workspace
rm(list = ls())

# Load necessary libraries
library(ggplot2)
library(RTMB)
source("R:/Science/Population Ecology Division/Shared/REFPT/REFPT_R_Project/0_ref_pt_helper_functions_231221.R")
source('OMftns_240820.R')

# Choose seed, stock, and number of sims to generate OM for
set.seed(666)
stock <- 'M'
sim_num <- 100

# Set productivity scenario
proj_scenario_name <- "strong_recovery" # Name for projection scenario
waa_change_hist <- 0.75 # For historical years... (final waa)/(initial waa) after 
                        # 50 years, with waa change following an exponential curve
waa_change_proj <- 1.25/0.75 # Same as above, but for projection years
mat_change_hist <- -0.02 # For historical years... fraction of an age to shift mat by
                         # each year, with - implying left shift
mat_change_proj <- +0.02 # Same as above, but for projection years
rec_timestep_factor <- c(1,1,1,1) # A length 4 vector of multipliers to set initial 
                                  # recruitment plus 3 step changes in recruitment (i.e., 
                                  # halfway through historical period, start of projections, 
                                  # and halfway through projections)
M_timestep_factor <- c(1,1,1,1) # A length 4 vector of (new M)/(Pons M), to set an 
                                # initial M plus 3 step changes in M (i.e., halfway 
                                # through historical period, start of projections, 
                                # and halfway through projections)

# Set historical F trajectory (linear climb to F_max for 2/3 of the time series, and 
# then a linear drop back down for the last third)
F_max <- 0.4 # Max historical F
F_initial <- 0.1 # Initial historical F
F_final <- 0.15 # Final historical F

#################################Setup inputs###################################

# Read in necessary inputs
D <- read.csv("R:/Science/Population Ecology Division/Shared/REFPT/REFPT_R_Project/data/Pons_inputs.csv")
pons_dat <- D[, stock]; names(pons_dat) <- D$Parameter
assign(names(pons_dat)[1], pons_dat[1])

# Setup age information... age zero recruits
max_age <- pons_dat['amax']
min_age <- 1
ages <- min_age:max_age
nages <- length(ages)

# Set other important values
nyears <- 50 # Number of historical years
initD <- 1 # Initial depletion
R0 <- 1000 # Initial recruitment (scales the NAA) in millions
myM <- c(rep(pons_dat['M']*M_timestep_factor[1], floor(nyears/2)),   # Historical M... 
         rep(pons_dat['M']*M_timestep_factor[2], ceiling(nyears/2))) # allows step change
h <- pons_dat['h'] # Steepness
fm <- c(seq(F_initial, F_max, length.out = ceiling(2/3*nyears)), # Historical F
        seq(F_max, F_final, length.out = floor(1/3*nyears)))

# Calculate initial LAA and WAA
laa <- pons_dat['Linf']*(1-exp(-pons_dat['VBk']*(ages-pons_dat['t0']))) # Units are cm
waa <- pons_dat['LWa']*laa^pons_dat['LWb'] # Units are kg

# Change WAA exponentially over time based on chosen factor
waa_matrix <- matrix(nrow = nages, ncol = nyears)
for (i in 1:nages) {
  waa_matrix[i, ] <- waa[i]*(1 + log(waa_change_hist)/50)^(0:(nyears-1))
}

# Calculate selectivity and maturity curves based on initial LAA
lmat_b1 <- log(9.5)/(pons_dat['L95']-pons_dat['L50'])
lmat_b0 <- log(9.5) - lmat_b1*pons_dat['L95']
lmat <- function(x){1/(1+exp(-lmat_b0-lmat_b1*x))}
mat <- lmat(laa)
lsel_b1 <- log(9.5)/(pons_dat['SL95']-pons_dat['SL50'])
lsel_b0 <- log(9.5) - lsel_b1*pons_dat['SL95']
lsel <- function(x){1/(1+exp(-lsel_b0-lsel_b1*x))}
sel <- lsel(laa)

# Fit the glms used to shift the selectivity and maturity curves by age
sel_mat_data <- as.data.frame(cbind(sel, mat, age = min_age:max_age))
sel_fit <- glm(sel ~ age, data = sel_mat_data, family = "binomial")
mat_fit <- glm(mat ~ age, data = sel_mat_data, family = "binomial")
sel_age <- plogis(predict(sel_fit))
mat_age <- plogis(predict(mat_fit))

# Shift mat over time based on chosen yearly amount, and calculate how the age at
# 50% mat changes
mat_matrix <- matrix(nrow = nages, ncol = nyears)
for (i in 1:nyears) {
  mat_matrix[, i] <- plogis(predict(mat_fit, 
                                    newdata = data.frame(age = min_age:max_age-
                                                          mat_change_hist*(i-1))))
}
amat50 <- apply(mat_matrix, 2, function(x) max(which(x < 0.5)))

# Calculate phi0 (in kg/rec), SSB0, and the alpha and beta SSR parameters
myphi0 <- sum(survivorship(nages, myM[1], sel_age, f=0)*waa_matrix[,1]*mat_matrix[,1])
SSB0 <- myphi0*R0 # Units are kt
BHa <- 4*h/((1-h)*myphi0)
BHb <- 1/R0*(BHa-1/myphi0)

# Calculate FAA as a combination of selectivity and fishing mortality
FAA <- as.data.frame(base::matrix(NA, nrow=nyears, ncol=nages))
for(i in 1:nrow(FAA)){
  for(j in 1:ncol(FAA)){
    FAA[i,j] <- exp(rnorm(1,log(fm[i]),pons_dat['sigF']))*sel_age[j]
  }
}

# Initialize the vector of SSB and the matrix of NAA... 
# NAA*waa is equal to SSB0 times initial depletion
c <- initD*SSB0/sum(R0*survivorship(nages, myM[1], sel_age, f=0)*waa_matrix[,1]*mat_matrix[,1])
NAA0 <- c*R0*survivorship(nages, myM[1], sel_age, f=0)
NAA <- base::array(0,dim=dim(FAA))
NAA[1,] <- NAA0
ssb <- numeric(nrow(NAA))
ssb[1] <- sum(NAA[1,]*mat_matrix[,1]*waa_matrix[,1])

# Fill out the rest of the NAA matrix and SSB vectors
for(i in 2:nrow(NAA)){
  for(j in 2:ncol(NAA)){
    NAA[i,j] <- NAA[i-1,j-1]*exp(-FAA[i-1,j-1]-myM[i-1])
    if(j==ncol(NAA)){
      NAA[i,j] <- NAA[i-1,j]*exp(-FAA[i-1,j]-myM[i-1]) + NAA[i,j]  # Sum of diag transition and plus group
    }  
  }
  # Recruitment
  NAA[i,1] <- exp(log(BHa*ssb[i-1]/(1+BHb*ssb[i-1])) + rnorm(1,0,pons_dat['sigrecdev']) )
  ssb[i] <- sum(NAA[i,]*mat_matrix[,i]*waa_matrix[,i])
}

# Obtain CAA in numbers (millions)
CAA <- FAA*NAA
catch <- numeric(nyears) # In kt
for(i in 1:length(catch)){catch[i] <- sum(CAA[i,]*waa_matrix[,i])}

# Plot some time series associated with the stock
DF <- data.frame(Year=1:nyears, SSB=ssb, REC=NAA[,1], 
                 Catch=catch)
ggplot(DF,aes(y=REC,x=SSB)) + geom_point() + theme_classic() + 
  geom_function(fun=function(x) BHa*x/(1+BHb*x), colour="blue",
                linetype=2, linewidth=1.2) +
  scale_x_continuous(limits = c(0, NA), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, NA), expand = c(0,0)) +
  xlab("SSB (kt)") + ylab("Recruitment (millions)") +
  geom_function(fun=function(x) 1/myphi0*x, colour="blue",linetype=2, linewidth=1.2) +
  theme(text = element_text(size = 16))
ggplot(DF,aes(y=Catch, x=Year)) + geom_path(linewidth = 1.2) + 
  theme_classic() + scale_x_continuous(limits = c(0, NA), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, NA), expand = c(0,0)) + ylab("Catch (kt)") +
  theme(text = element_text(size = 16))
ggplot(DF,aes(y=SSB, x=Year)) + geom_path(linewidth = 1.2) + 
  theme_classic() + scale_x_continuous(limits = c(0, NA), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, NA), expand = c(0,0)) +  
  ylab("SSB (kt)") + theme(text = element_text(size = 16))
ggplot(DF,aes(y=REC, x=Year)) + geom_path(linewidth = 1.2) + 
  theme_classic() + scale_x_continuous(limits = c(0, NA), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, NA), expand = c(0,0)) + 
  ylab("Recruitment (millions)") + theme(text = element_text(size = 16))

#####################Fit assessment model and generate OM#######################

# Function for SSB for one year in kt
annualssb <- function(N, SW, MO, y){
  ret <- 0 
  for(a in 2:nages){ret <- ret + SW[y,a]*MO[y,a]*N[y,a]}
  ret
}

# Setup some important quantities
minAge <- 1
nages <- length(minAge:max_age)
maxYear <- 50
minYear <- maxYear-nyears+1
obs_index_sd <- 0.25

# Grab the fishery NAA (in millions) and the SSB index with obs error
CAA[is.na(CAA) | CAA==0] <- NA
SSBO <- exp(log(ssb)+rnorm(length(ssb),0,obs_index_sd)) # kt

# Create the data list
dat <- list()
dat$logobs <- log(CAA) # Fishery NAA in millions
dat$logind <- log(SSBO) # Index of SSB in kt
dat$year <- minYear:maxYear
dat$age <-  minAge:max_age
dat$M <- array(NA,dim=c(nyears,nages))
for(i in 1:ncol(dat$M)){dat$M[,i] <- myM}
dat$SW <- t(waa_matrix) # Stock waa in kg
dat$CW <- dat$SW # Catch waa in kg
dat$MO <- t(mat_matrix) # Maturity at age
dat$SRR <- "Beverton-Holt"
dat$phi0 <- sum(survivorship(n_ages = dim(dat$M)[2], # Phi0 in kg/recruit 
                             m = colSums(dat$M[1:amat50[1],])/amat50[1], 
                             selectivity = rep(0,dim(dat$M)[2]),
                             f = 0)*
                  colSums(dat$SW[1:amat50[1],])/amat50[1]*
                  colSums(dat$MO[1:amat50[1],])/amat50[1])
dat$logobs <- dat$logobs + 
  array(rnorm(ncol(dat$logobs)*nrow(dat$logobs),0,0.1),dim=dim(dat$logobs)) # Add obs error

# Set up prior for h
hmu <- h*1.25-0.25
hsd <- 0.05 * 1.25
halpha <- hmu*((hmu*(1-hmu))/hsd^2-1)
hbeta <- (1-hmu)*((hmu*(1-hmu))/hsd^2-1)

# Set up parameter list
par <- list()
par$logN1Y <- numeric(nages) # Initial NAA
par$logrecs <- numeric(nyears-1) # Annual recruitment deviations
par$logsdO <- 0 # Observation error on CAA (applies to all NAA except recruitment)
par$logq <- 0  # Catchability coefficient
par$logFM <- rep(0,nyears) # Vector of annual fishing mortality
par$logsdS <- 0 # Observation error on survey index
par$trans_a50 <- 0; par$trans_k <- 0  # Fishery selectivity parameters needed for logistic                      
par$logR0 <- 0; par$trans_h <- 0; par$logsdR <- 0 # SRR parameters

# Objective function
f<-function(par){
  getAll(par,dat)

  # Sel parameters
  k <- exp(trans_k)
  a50 <- plogis(trans_a50)*max_age
  
  # SRR parameters
  R0 <- exp(logR0) # millions
  parh <- 0.8*plogis(trans_h)+0.2
  
  # Initialize joint negative log likelihood
  jnll <- 0
  
  # Fill in logN
  logN <- base::matrix(0, ncol=nages, nrow=nyears)
  # Initial year NAA (without recruitment)
  logN[1,] <- logN1Y
  for(y in 2:nyears) {
    logN[y,1] <- logN[1,1] + sum(logrecs[2:y-1])
  }
  
  # Transitions to fill in logN
  for(y in 2:nyears){
    Ffor <- exp(sum(logFM[1:(y-1)])) # Annual F from previous year used for transition
    for (a in 2:nages) { # a is column in logN [1 to max_age] and the age is the actual age
      selectivity <- calculate_selectivity(a-1, a50, k)
      logN[y, a] <- logN[y-1, a-1] - Ffor*selectivity - M[y-1, a-1] # NAA transition diagonals
      
      if (a == nages) { # nages-1 is the maximum age and goes into the selectivity function
        selectivity_plus <- calculate_selectivity(max_age, a50, k)
        logN[y, a] <- log(exp(logN[y, a])+
                            exp(logN[y-1,a] - Ffor*selectivity_plus-M[y-1,a])) # NAA transition plus group
      }
    }
  }
  for(y in 1:nyears){
    curF <- exp(sum(logFM[1:y]))
    for(a in 1:nages){ # Obs error on NAA for all years but 
      # exclude ages with no observations in CAA 
      FAA <- curF*calculate_selectivity(a, a50, k)
      jnll <- jnll - dnorm(logobs[y,a], 
                           log(FAA)-log(FAA+M[y,a])+log(1-exp(-FAA-M[y,a]))+logN[y,a], 
                           exp(logsdO), log = TRUE)
    }
  }
  for(y in 2:nyears){
    ssby <- sum(SW[y,]*MO[y,]*exp(logN[y-1,])) # In kt
    pred <- log(4*R0*parh*ssby/ ((1-parh)*R0*phi0+(5*parh-1)*ssby)) 
    jnll <- jnll - dnorm(logN[y,1], pred, exp(logsdR), log = TRUE) # Age 0's in NAA 
  }
  
  # Prior penalty for h
  pri_y <- (parh-0.2)/0.8
  jnll <- jnll - dbeta(pri_y, halpha, hbeta, log = T) - log(pri_y-pri_y^2)
  
  # Index contribution to likelihood
  for(y in 1:nyears){
    pred <- log(sum(exp(logN[y,])*MO[y,]*SW[y,])) # Index is SSB in kt
    jnll <- jnll - dnorm(log(exp(logq)*exp(logind[y])),pred,exp(logsdS),log=TRUE)
  }
  
  # Report important quantities
  ssbFUN <- function(N, SW, MO){
    nrow <- nrow(N)
    ncol <- ncol(N)
    ret <- numeric(nrow)
    for(y in 1:nrow){
      for(a in 1:ncol){
        ret[y] = ret[y]+SW[y,a]*MO[y,a]*N[y,a]
      }
    }
    ret
  }  
  ssb <- ssbFUN(N=exp(logN),SW,MO)
  ADREPORT(ssb)
  ADREPORT(log(ssb))
  REPORT(logN)
  
  jnll
}    

# Fit model
obj <- MakeADFun(f, par, silent=T)
opt <- nlminb(obj$par, obj$fn, obj$gr, control=list(eval.max=1000, iter.max=1000))
sdr <- RTMB::sdreport(obj)

# Store model outputs
OL <- list()
OL$recs <- exp(obj$report()$logN[,1]) # Rec in millions
OL$ssb <- sdr$value[1:nyears] # SSB in kt
OL$R0 <- exp(opt$par[names(opt$par)=="logR0"])
OL$h <- 0.8*plogis(opt$par[names(opt$par)=="trans_h"])+0.2
OL$alpha <- 4*OL$h/(dat$phi0-OL$h*dat$phi0)  # kg/rec
OL$beta <- 1/(OL$R0)*(OL$alpha - 1/(dat$phi0))
OL$N <- exp(obj$report()$logN)
OL$q <- exp(opt$par[names(opt$par)=="logq"])
OL$k <- exp(opt$par[names(opt$par)=="trans_k"])
OL$a50 <- plogis(opt$par[names(opt$par)=="trans_a50"])*max_age
OL$phi0 <- dat$phi0 # kg/rec or kt/million rec
OL$logSSB_terminal_se <- sdr$sd[names(sdr$value)=="log(ssb)"][sum(
  names(sdr$value)=="log(ssb)")]
OL$conv <- ifelse(opt$convergence == 0 & sdr$pdHess == TRUE & 
                    !(is.na(OL$logSSB_terminal_se)) == TRUE, 0, 1)

# Generate OM using model fit and set some parameters
OM <- SCA2OM(obj, dat, nsim = sim_num, myseed = 666)
OM@Misc <- list(stock = stock,
                proj_scenario_name = proj_scenario_name,
                waa_change_hist = waa_change_hist, 
                waa_change_proj = waa_change_proj,
                mat_change_hist = mat_change_hist,
                mat_change_proj = mat_change_proj,
                rec_timestep_factor = rec_timestep_factor,
                M_timestep_factor = M_timestep_factor,
                F_max = F_max,
                F_initial = F_initial,
                F_final = F_final)
OM <- obs_change(OM)
OM@Iobs <- rep(obs_index_sd,2)

### Alter projection scenario for OM as required###

OM <- cpars_pro5(OM) # Set baseline for productivity components 
                     # in projections (i.e., average of last 5 historical years)
OM <- waa_trend_proj(OM, annual_rate = 1+log(waa_change_proj)/50)  # Change waa 
                                                                   # exponentially 
                                                                   # in the projections 
                                                                   # based on specified 
                                                                   # rate
OM <- mat_shift_proj(OM, annual_shift = mat_change_proj) # Shift projection mat yearly
                                                         # by specified rate

# Specified step changes in recruitment (initial recruitment, halfway through 
# historical period, at start of projection years, and halfway through projection years)
Perr_y_original <- OM@cpars$Perr_y
OM@cpars$Perr_y[,(OM@maxage-2):(OM@maxage+OM@nyears+OM@proyears)] <- 
  rec_timestep_factor[1]*Perr_y_original[,(OM@maxage-2):(OM@maxage+OM@nyears+OM@proyears)]
OM@cpars$Perr_y[,(OM@maxage+OM@nyears-2-ceiling(OM@nyears/2)):(OM@maxage+OM@nyears+OM@proyears)] <- 
  rec_timestep_factor[2]*Perr_y_original[,(OM@maxage+OM@nyears-2-ceiling(OM@nyears/2)):(OM@maxage+OM@nyears+OM@proyears)]
OM@cpars$Perr_y[,(OM@maxage+OM@nyears-2):(OM@maxage+OM@nyears+OM@proyears)] <- 
  rec_timestep_factor[3]*Perr_y_original[,(OM@maxage+OM@nyears-2):(OM@maxage+OM@nyears+OM@proyears)]
OM@cpars$Perr_y[,(OM@maxage+OM@nyears-2+ceiling(OM@proyears/2)):(OM@maxage+OM@nyears+OM@proyears)] <- 
  rec_timestep_factor[4]*Perr_y_original[,(OM@maxage+OM@nyears-2+ceiling(OM@proyears/2)):(OM@maxage+OM@nyears+OM@proyears)]

# Specified step changes in M (at start of projection years, and halfway through)
for(i in 1:OM@nsim)
{
  OM@cpars$M_ageArray[i,2:(OM@maxage+1),
                      (OM@nyears+1):(OM@proyears+OM@nyears)] <- 
    M_timestep_factor[3]*pons_dat['M']
  OM@cpars$M_ageArray[i,2:(OM@maxage+1),
                      (OM@nyears+1+ceiling(OM@proyears/2)):(OM@nyears+OM@proyears)] <- 
    M_timestep_factor[4]*pons_dat['M']
}

# Set catch waa = stock waa and save OM
OM@cpars$Wt_age_C <- OM@cpars$Wt_age
saveRDS(OM,paste0(getwd(),"/OM/",stock,"_",proj_scenario_name,".rds"))