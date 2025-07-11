#
# SCA model-fitting function used by MPs
#
################################################################################

# Load library
library(RTMB)

SCA_fit <- function(x, Data){
  
  # SSB for one year in kt
  annualssb <- function(N, SW, MO, y){
    ret <- 0 
    for(a in 2:nages){ret <- ret + SW[y,a]*MO[y,a]*N[y,a]}
    ret
  }
  
  # Setup some important quantities
  nyears <- ncol(Data@Ind)
  minAge <- 1 # To remove age 0 from the OM for the EM fit
  maxAge <- Data@MaxAge
  nages <- maxAge-minAge+1
  maxYear <- Data@OM$CurrentYr[1]
  minYear <- OM@CurrentYr-nyears+1
  
  # Get catch NAA in millions
  pNAA <- Data@CAA[x,,minAge:maxAge+1]/rowSums(Data@CAA[x,,minAge:maxAge+1])
  pNAA_waa <- pNAA * t(Data@Misc$FleetPars$Wt_age_C[x,minAge:maxAge+1,1:nyears])
  pCAA <- pNAA_waa/rowSums(pNAA_waa)
  BAA <- pCAA*Data@Cat[x,]
  CAA_obs <- BAA/t(Data@Misc$FleetPars$Wt_age_C[x,minAge:maxAge+1,1:nyears])    
  CAA_obs[is.nan(CAA_obs) | CAA_obs==0] <- NA
  
  # Get the index of SSB with obs error
  theSSB <- colSums(Data@Misc$StockPars$SSB[x,,,1]*2)
  if(nyears > dim(Data@Misc$StockPars$SSB)[3]){
    theSSB <- c(theSSB, 
                colSums(Data@Misc$StockPars$SSB_P[x,,1:(1+nyears-dim(Data@Misc$StockPars$SSB)[3]),1]*2))
  }
  for(i in 1:length(theSSB)){ # For all years
    theSSB[i] <- theSSB[i]*exp(rnorm(1,0,Data@Obs$Isd[x]))
  }  
  SSB_obs <- theSSB
  
  # Create the data list
  dat <- list()
  dat$logobs <- log(CAA_obs) # Log fishery NAA for simulation x 
  dat$logind <- log(SSB_obs) # Index of SSB in kt for simulation x
  dat$year <- minYear:maxYear
  dat$age <-  minAge:maxAge
  dat$M <- t(Data@Misc$StockPars$M_ageArray[x,2:(maxAge+1),1:nyears]) # M at age; remove age 0
  dat$SW <- t(Data@Misc$StockPars$Wt_age[x,2:(maxAge+1),1:nyears])    # Stock waa in kg; remove age 0
  dat$CW <- t(Data@Misc$FleetPars$Wt_age_C[x,2:(maxAge+1),1:nyears])  # Catch waa in kg; remove age 0
  dat$MO <- t(Data@Misc$StockPars$Mat_age[x,2:(maxAge+1),1:nyears])   # Maturity at age; remove age 0
  dat$SRR <- "Beverton-Holt"
  dat$phi0 <- sum(survivorship(n_ages = dim(dat$M)[2], # Phi0 in kg/recruit
                               m = colSums(dat$M[1:5,])/5, selectivity = rep(0,dim(dat$M)[2]),
                               f = 0)*colSums(dat$SW[1:5,])/5*colSums(dat$MO[1:5,])/5)

  # Set up prior for h
  hmu <- Data@Misc$StockPars$hs[x]*1.25-0.25
  hsd <- 0.025 *1.25
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
  par$missing <- numeric(sum(is.na(dat$logobs))) # Zeros as missing data
  
  # Objective function
  f<-function(par){
    getAll(par,dat)
    
    # Set up missing observations as random effects
    logobs[is.na(logobs)] <- missing
    
    # Sel parameters
    k <- exp(trans_k)
    a50 <- plogis(trans_a50)*maxAge
    
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
      for (a in 2:nages) { # a is column in logN [1 to maxAge] and the age is the actual age
        selectivity <- calculate_selectivity(a-1, a50, k)
        logN[y, a] <- logN[y-1, a-1] - Ffor*selectivity - M[y-1, a-1] # NAA transition diagonals
        
        if (a == nages) { # nages-1 is the maximum age and goes into the selectivity function
          selectivity_plus <- calculate_selectivity(maxAge, a50, k)
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
  
  # Fit model (with fixed q = 1)
  obj <- MakeADFun(f, par, random = "missing", silent=T,
                   map = list(logq = factor(NA)))
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
  OL$q <- 1 # exp(opt$par[names(opt$par)=="logq"])
  OL$k <- exp(opt$par[names(opt$par)=="trans_k"])
  OL$a50 <- plogis(opt$par[names(opt$par)=="trans_a50"])*maxAge
  OL$phi0 <- dat$phi0 # kg/rec or kt/million rec
  OL$logSSB_terminal_se <- sdr$sd[names(sdr$value)=="log(ssb)"][sum(
    names(sdr$value)=="log(ssb)")]
  OL$conv <- ifelse(opt$convergence == 0 & sdr$pdHess == TRUE & 
                      !(is.na(OL$logSSB_terminal_se)) == TRUE, 0, 1)
  
  # Return model outputs
  return (OL)
}