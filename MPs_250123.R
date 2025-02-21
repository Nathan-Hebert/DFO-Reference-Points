#
# Setup MPs for use with MSEtool
#
################################################################################

# Load necessary libraries
source("SCA_function_240827.R")

# Uses MSY-based reference pts calculated with first 5 years
MP_first5 <- function (x, Data, reps = 1){
  # Determine year
  y <- max(Data@Year) - Data@LHYear + 1 

  # If new model needs to be fit (year 1, 6, 11, etc.)
  if(y%%5==1)
  {
    # Fit model until it converges
    compiler::enableJIT(0)
    mySCA <- try(SCA_fit(x, Data), silent = TRUE)
    while(inherits(mySCA, "try-error") || mySCA$conv != 0) {
      mySCA <- try(SCA_fit(x, Data), silent = TRUE)
    }
    
    # Calculate SSB_MSY, F_MSY using first 5 years
    yrstart <- 1
    yrend <- 5
    model_sel <- 1/(1+exp((mySCA$a50-c(0:Data@MaxAge))/mySCA$k))
    DATA <- data.frame(years = rep(yrstart:yrend, each = Data@MaxAge + 1),
                       ages = rep(0:Data@MaxAge, length(yrstart:yrend)),
                       waa = as.vector(Data@Misc$StockPars$Wt_age[x,,yrstart:yrend]), 
                       mat = as.vector(Data@Misc$StockPars$Mat_age[x,,yrstart:yrend]),
                       m = as.vector(Data@Misc$StockPars$M_ageArray[x,,yrstart:yrend]),
                       selectivity = rep(model_sel, length(yrstart:yrend)))
    RP_DF <- average_time_period_RPs(DATA = DATA, years = "years", 
                                     ages = "ages", m = "m", 
                                     selectivity = "selectivity", 
                                     waa = "waa", mat = "mat", 
                                     srr = "beverton-holt", 
                                     alpha = mySCA$alpha, 
                                     beta = mySCA$beta,
                                     SPR = 0.4, 
                                     average_yrs = list(yrstart:yrend))
 
    # Store SSB_MSY, F_MSY, and model-estimated q until the next model fit
    CP_list[["F_MSY_First5"]][y:(y+4),x] <<- RP_DF$F_MSY
    CP_list[["SSB_MSY_First5"]][y:(y+4),x] <<- RP_DF$SSB_MSY
    CP_list[["q_First5"]][y:(y+4),x] <<- mySCA$q
  }

  Rec <- new("Rec")
  # Load in stored values needed for HCR
  fcp <- CP_list[["F_MSY_First5"]][y,x]
  q <- CP_list[["q_First5"]][y,x]  
  lcp <- 0.4*CP_list[["SSB_MSY_First5"]][y,x]
  ucp <- 0.8*CP_list[["SSB_MSY_First5"]][y,x]    
    
  # Grab last index value (SSB in kt + obs error) and apply q
  indexEMunits <- q*sum(Data@Misc$StockPars$SSB_P[x,,y,1]*2)*exp(rnorm(1,0,Data@Obs$Isd[x]))
    
  # Implement HCR using above values
  slope <- fcp/(ucp-lcp) # Slope of the line on the HCR between (LCP,F=0) and 
                         # (UCP,F=Fmsy)
  if(indexEMunits < lcp){fmsy_adj <- 1e-6}
  if(indexEMunits > ucp){fmsy_adj <- fcp}
  if(indexEMunits >= lcp & indexEMunits <= ucp){fmsy_adj <- 
    (slope*indexEMunits-slope*lcp)}
  Rec@Effort <- fmsy_adj/Data@OM$FinF[x]
  
  print(x)
  print(y)
  print(Rec)
  return(Rec)
}  
class(MP_first5) <- "MP"

# Uses MSY-based reference pts calculated with last 5 years
MP_last5 <- function (x, Data, reps = 1){
  # Determine year
  y <- max(Data@Year) - Data@LHYear + 1 
  
  # If new model needs to be fit (year 1, 6, 11, etc.)
  if(y%%5==1)
  {
    # Fit model until it converges
    compiler::enableJIT(0)
    mySCA <- try(SCA_fit(x, Data), silent = TRUE)
    while(inherits(mySCA, "try-error") || mySCA$conv != 0){
      mySCA <- try(SCA_fit(x, Data), silent = TRUE)
    }
    
    # Calculate SSB_MSY, F_MSY using last 5 years
    yrstart <- length(Data@Year)-4
    yrend <- length(Data@Year)
    model_sel <- 1/(1+exp((mySCA$a50-c(0:Data@MaxAge))/mySCA$k))
    DATA <- data.frame(years = rep(yrstart:yrend, each = Data@MaxAge + 1),
                       ages = rep(0:Data@MaxAge, length(yrstart:yrend)),
                       waa = as.vector(Data@Misc$StockPars$Wt_age[x,,yrstart:yrend]), 
                       mat = as.vector(Data@Misc$StockPars$Mat_age[x,,yrstart:yrend]),
                       m = as.vector(Data@Misc$StockPars$M_ageArray[x,,yrstart:yrend]),
                       selectivity = rep(model_sel, length(yrstart:yrend)))
    RP_DF <- average_time_period_RPs(DATA = DATA, years = "years", 
                                     ages = "ages", m = "m", 
                                     selectivity = "selectivity", 
                                     waa = "waa", mat = "mat", 
                                     srr = "beverton-holt", 
                                     alpha = mySCA$alpha, 
                                     beta = mySCA$beta,
                                     SPR = 0.4, 
                                     average_yrs = list(yrstart:yrend))
      
    # Store SSB_MSY, F_MSY, and model-estimated q until the next model fit
    CP_list[["F_MSY_Last5"]][y:(y+4),x] <<- RP_DF$F_MSY
    CP_list[["SSB_MSY_Last5"]][y:(y+4),x] <<- RP_DF$SSB_MSY
    CP_list[["q_Last5"]][y:(y+4),x] <<- mySCA$q
  }

  Rec <- new("Rec")
  # Load in stored values needed for HCR
  fcp <- CP_list[["F_MSY_Last5"]][y,x]
  q <- CP_list[["q_Last5"]][y,x]  
  lcp <- 0.4*CP_list[["SSB_MSY_Last5"]][y,x]
  ucp <- 0.8*CP_list[["SSB_MSY_Last5"]][y,x]    
    
  # Grab last index value (SSB in kt + obs error) and apply q
  indexEMunits <- q*sum(Data@Misc$StockPars$SSB_P[x,,y,1]*2)*exp(rnorm(1,0,Data@Obs$Isd[x]))
    
  # Implement HCR using above values
  slope <- fcp/(ucp-lcp) # Slope of the line on the HCR between (LCP,F=0) and 
                         # (UCP,F=Fmsy)
  if(indexEMunits < lcp){fmsy_adj <- 1e-6}
  if(indexEMunits > ucp){fmsy_adj <- fcp}
  if(indexEMunits >= lcp & indexEMunits <= ucp){fmsy_adj <- 
                                                (slope*indexEMunits-slope*lcp)}
  Rec@Effort <- fmsy_adj/Data@OM$FinF[x]

  print(x)
  print(y)
  print(Rec)
  return(Rec)
}  
class(MP_last5) <- "MP"

# Uses MSY-based reference pts calculated with entire time series
MP_entirets <- function (x, Data, reps = 1){
  # Determine year
  y <- max(Data@Year) - Data@LHYear + 1 
  
  # If new model needs to be fit (year 1, 6, 11, etc.)
  if(y%%5==1)
  {
    # Fit model until it converges
    compiler::enableJIT(0)
    mySCA <- try(SCA_fit(x, Data), silent = TRUE)
    while(inherits(mySCA, "try-error") || mySCA$conv != 0){
      mySCA <- try(SCA_fit(x, Data), silent = TRUE)
    }
    
    # Calculate SSB_MSY, F_MSY using entire time series
    yrstart <- 1
    yrend <- length(Data@Year)
    model_sel <- 1/(1+exp((mySCA$a50-c(0:Data@MaxAge))/mySCA$k))
    DATA <- data.frame(years = rep(yrstart:yrend, each = Data@MaxAge + 1),
                       ages = rep(0:Data@MaxAge, length(yrstart:yrend)),
                       waa = as.vector(Data@Misc$StockPars$Wt_age[x,,yrstart:yrend]), 
                       mat = as.vector(Data@Misc$StockPars$Mat_age[x,,yrstart:yrend]),
                       m = as.vector(Data@Misc$StockPars$M_ageArray[x,,yrstart:yrend]),
                       selectivity = rep(model_sel, length(yrstart:yrend)))
    RP_DF <- average_time_period_RPs(DATA = DATA, years = "years", 
                                     ages = "ages", m = "m", 
                                     selectivity = "selectivity", 
                                     waa = "waa", mat = "mat", 
                                     srr = "beverton-holt", 
                                     alpha = mySCA$alpha, 
                                     beta = mySCA$beta,
                                     SPR = 0.4, 
                                     average_yrs = list(yrstart:yrend))
      
    # Store SSB_MSY, F_MSY, and model-estimated q until the next model fit
    CP_list[["F_MSY_EntireTS"]][y:(y+4),x] <<- RP_DF$F_MSY
    CP_list[["SSB_MSY_EntireTS"]][y:(y+4),x] <<- RP_DF$SSB_MSY
    CP_list[["q_EntireTS"]][y:(y+4),x] <<- mySCA$q
  }

  Rec <- new("Rec")
  # Load in stored values needed for HCR
  fcp <- CP_list[["F_MSY_EntireTS"]][y,x]
  q <- CP_list[["q_EntireTS"]][y,x]  
  lcp <- 0.4*CP_list[["SSB_MSY_EntireTS"]][y,x]
  ucp <- 0.8*CP_list[["SSB_MSY_EntireTS"]][y,x]  
    
  # Grab last index value (SSB in kt + obs error) and apply q
  indexEMunits <- q*sum(Data@Misc$StockPars$SSB_P[x,,y,1]*2)*exp(rnorm(1,0,Data@Obs$Isd[x]))
    
  # Implement HCR using above values
  slope <- fcp/(ucp-lcp) # Slope of the line on the HCR between (LCP,F=0) and 
                         # (UCP,F=Fmsy)
  if(indexEMunits < lcp){fmsy_adj <- 1e-6}
  if(indexEMunits > ucp){fmsy_adj <- fcp}
  if(indexEMunits >= lcp & indexEMunits <= ucp){fmsy_adj <- 
    (slope*indexEMunits-slope*lcp)}
  Rec@Effort <- fmsy_adj/Data@OM$FinF[x]

  print(x)
  print(y)
  print(Rec)
  return(Rec)
}  
class(MP_entirets) <- "MP"

# Dynamic HCR based on MSY reference pts
MP_mixed <- function (x, Data, reps = 1){
  # Determine year
  y <- max(Data@Year) - Data@LHYear + 1 
  
  # If new model needs to be fit (year 1, 6, 11, etc.)
  if(y%%5==1)
  {
    # Fit model until it converges
    compiler::enableJIT(0)
    mySCA <- try(SCA_fit(x, Data), silent = TRUE)
    while(inherits(mySCA, "try-error") || mySCA$conv != 0){
      mySCA <- try(SCA_fit(x, Data), silent = TRUE)
    }
    
    # Calculate SSB_MSY, F_MSY using last 5 years
    yrstart <- length(Data@Year)-4
    yrend <- length(Data@Year)
    model_sel <- 1/(1+exp((mySCA$a50-c(0:Data@MaxAge))/mySCA$k))
    DATA <- data.frame(years = rep(yrstart:yrend, each = Data@MaxAge + 1),
                       ages = rep(0:Data@MaxAge, length(yrstart:yrend)),
                       waa = as.vector(Data@Misc$StockPars$Wt_age[x,,yrstart:yrend]), 
                       mat = as.vector(Data@Misc$StockPars$Mat_age[x,,yrstart:yrend]),
                       m = as.vector(Data@Misc$StockPars$M_ageArray[x,,yrstart:yrend]),
                       selectivity = rep(model_sel, length(yrstart:yrend)))
    RP_DF <- average_time_period_RPs(DATA = DATA, years = "years", 
                                     ages = "ages", m = "m", 
                                     selectivity = "selectivity", 
                                     waa = "waa", mat = "mat", 
                                     srr = "beverton-holt", 
                                     alpha = mySCA$alpha, 
                                     beta = mySCA$beta,
                                     SPR = 0.4, 
                                     average_yrs = list(yrstart:yrend))
      
    # Calculate SSB_MSY using first 5 years
    yrstart <- 1
    yrend <- 5
    DATA <- data.frame(years = rep(yrstart:yrend, each = Data@MaxAge + 1),
                       ages = rep(0:Data@MaxAge, length(yrstart:yrend)),
                       waa = as.vector(Data@Misc$StockPars$Wt_age[x,,yrstart:yrend]), 
                       mat = as.vector(Data@Misc$StockPars$Mat_age[x,,yrstart:yrend]),
                       m = as.vector(Data@Misc$StockPars$M_ageArray[x,,yrstart:yrend]),
                       selectivity = rep(model_sel, length(yrstart:yrend)))
    RP_DF2 <- average_time_period_RPs(DATA = DATA, years = "years", 
                                      ages = "ages", m = "m", 
                                      selectivity = "selectivity", 
                                      waa = "waa", mat = "mat", 
                                      srr = "beverton-holt", 
                                      alpha = mySCA$alpha, 
                                      beta = mySCA$beta,
                                      SPR = 0.4, 
                                      average_yrs = list(yrstart:yrend))
      
    # Store, until the next model fit: last 5 years SSB_MSY, 
    # last 5 years F_MSY, first 5 years SSB_MSY, and model-estimated q
    CP_list[["F_MSY_Last5Mixed"]][y:(y+4),x] <<- RP_DF$F_MSY
    CP_list[["SSB_MSY_Last5Mixed"]][y:(y+4),x] <<- RP_DF$SSB_MSY
    CP_list[["SSB_MSY_First5Mixed"]][y:(y+4),x] <<- RP_DF2$SSB_MSY
    CP_list[["q_Mixed"]][y:(y+4),x] <<- mySCA$q
  }

  Rec <- new("Rec")
  # Load in stored values needed for HCR
  fcp <- CP_list[["F_MSY_Last5Mixed"]][y,x]
  q <- CP_list[["q_Mixed"]][y,x]  
  lcp <- 0.4*CP_list[["SSB_MSY_First5Mixed"]][y,x]
  ucp <- 0.8*CP_list[["SSB_MSY_Last5Mixed"]][y,x]
    
  # Grab last index value (SSB in kt + obs error) and apply q
  indexEMunits <- q*sum(Data@Misc$StockPars$SSB_P[x,,y,1]*2)*exp(rnorm(1,0,Data@Obs$Isd[x]))
    
  # Implement HCR using above values
  slope <- fcp/(ucp-lcp) # Slope of the line on the HCR between (LCP,F=0) and 
                         # (UCP,F=Fmsy)
  if(indexEMunits < lcp){fmsy_adj <- 1e-6}
  if(indexEMunits > ucp){fmsy_adj <- fcp}
  if(indexEMunits >= lcp & indexEMunits <= ucp){fmsy_adj <- 
    (slope*indexEMunits-slope*lcp)}
  Rec@Effort <- fmsy_adj/Data@OM$FinF[x]

  print(x)
  print(y)
  print(Rec)
  return(Rec)
}  
class(MP_mixed) <- "MP"

# Dynamic HCR based on MSY reference pts, with BPA as ucp
MP_mixed_BPA <- function (x, Data, reps = 1){
  # Determine year
  y <- max(Data@Year) - Data@LHYear + 1 
  
  # If new model needs to be fit (year 1, 6, 11, etc.)
  if(y%%5==1)
  {
    # Fit model until it converges
    compiler::enableJIT(0)
    mySCA <- try(SCA_fit(x, Data), silent = TRUE)
    while(inherits(mySCA, "try-error") || mySCA$conv != 0){
      mySCA <- try(SCA_fit(x, Data), silent = TRUE)
    }
    
    # Calculate F_MSY using last 5 years
    yrstart <- length(Data@Year)-4
    yrend <- length(Data@Year)
    model_sel <- 1/(1+exp((mySCA$a50-c(0:Data@MaxAge))/mySCA$k))
    DATA <- data.frame(years = rep(yrstart:yrend, each = Data@MaxAge + 1),
                       ages = rep(0:Data@MaxAge, length(yrstart:yrend)),
                       waa = as.vector(Data@Misc$StockPars$Wt_age[x,,yrstart:yrend]), 
                       mat = as.vector(Data@Misc$StockPars$Mat_age[x,,yrstart:yrend]),
                       m = as.vector(Data@Misc$StockPars$M_ageArray[x,,yrstart:yrend]),
                       selectivity = rep(model_sel, length(yrstart:yrend)))
    RP_DF <- average_time_period_RPs(DATA = DATA, years = "years", 
                                     ages = "ages", m = "m", 
                                     selectivity = "selectivity", 
                                     waa = "waa", mat = "mat", 
                                     srr = "beverton-holt", 
                                     alpha = mySCA$alpha, 
                                     beta = mySCA$beta,
                                     SPR = 0.4, 
                                     average_yrs = list(yrstart:yrend))
    
    # Calculate SSB_MSY using first 5 years
    yrstart <- 1
    yrend <- 5
    DATA <- data.frame(years = rep(yrstart:yrend, each = Data@MaxAge + 1),
                       ages = rep(0:Data@MaxAge, length(yrstart:yrend)),
                       waa = as.vector(Data@Misc$StockPars$Wt_age[x,,yrstart:yrend]), 
                       mat = as.vector(Data@Misc$StockPars$Mat_age[x,,yrstart:yrend]),
                       m = as.vector(Data@Misc$StockPars$M_ageArray[x,,yrstart:yrend]),
                       selectivity = rep(model_sel, length(yrstart:yrend)))
    RP_DF2 <- average_time_period_RPs(DATA = DATA, years = "years", 
                                      ages = "ages", m = "m", 
                                      selectivity = "selectivity", 
                                      waa = "waa", mat = "mat", 
                                      srr = "beverton-holt", 
                                      alpha = mySCA$alpha, 
                                      beta = mySCA$beta,
                                      SPR = 0.4, 
                                      average_yrs = list(yrstart:yrend))
    
    # Store, until the next model fit: last 5 years F_MSY, 
    # first 5 years SSB_MSY, terminal year logSSB sigma, and model-estimated q
    CP_list[["F_MSY_Last5_MixedBPA"]][y:(y+4),x] <<- RP_DF$F_MSY
    CP_list[["SSB_MSY_First5_MixedBPA"]][y:(y+4),x] <<- RP_DF2$SSB_MSY
    CP_list[["q_MixedBPA"]][y:(y+4),x] <<- mySCA$q
    CP_list[["SSB_Sigma_MixedBPA"]][y:(y+4),x] <<- mySCA$logSSB_terminal_se
  }
  
  Rec <- new("Rec")
  # Load in stored values needed for HCR
  fcp <- CP_list[["F_MSY_Last5_MixedBPA"]][y,x]
  q <- CP_list[["q_MixedBPA"]][y,x]  
  lcp <- 0.4*CP_list[["SSB_MSY_First5_MixedBPA"]][y,x]
  ucp <- lcp*exp(1.645*sqrt(CP_list[["SSB_Sigma_MixedBPA"]][y,x]^2+Data@Obs$Isd[x]^2))
  
  # Grab last index value (SSB in kt + obs error) and apply q
  indexEMunits <- q*sum(Data@Misc$StockPars$SSB_P[x,,y,1]*2)*exp(rnorm(1,0,Data@Obs$Isd[x]))
  
  # Implement HCR using above values
  slope <- fcp/(ucp-lcp) # Slope of the line on the HCR between (LCP,F=0) and 
  # (UCP,F=Fmsy)
  if(indexEMunits < lcp){fmsy_adj <- 1e-6}
  if(indexEMunits > ucp){fmsy_adj <- fcp}
  if(indexEMunits >= lcp & indexEMunits <= ucp){fmsy_adj <- 
    (slope*indexEMunits-slope*lcp)}
  Rec@Effort <- fmsy_adj/Data@OM$FinF[x]
  
  print(x)
  print(y)
  print(Rec)
  return(Rec)
}  
class(MP_mixed_BPA) <- "MP"