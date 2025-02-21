#
# Some miscellaneous helper functions
#
################################################################################

# Load necessary libraries
library(dplyr); library(MARS)

# Function returns the length of an organism at a specific time t using the 
# von Bertalanffy growth equation, which takes parameters for the organism's 
# asymptotic length (l_inf), growth coefficient (k), and starting time (t0).
#
# t = the vector of timepoints with which to output lengths
# l_inf = the asymptotic length parameter
# k = the growth coefficient
# t0 = the starting time parameter
von_bertalanffy <- function(t, l_inf, k, t0) 
{
  l_t <- l_inf * (1 - exp(-k * (t - t0)))
}

# Function gives a vector that describes survivorship-at-age.
#
# n_ages = the number of age classes considered
# m = the vector or single numeric value that describes mortality-at-age
# selectivity = a vector describing selectivity-at-age
# f = a single numeric fishing mortality value
survivorship <- function(n_ages, m, selectivity, f)
{
  # Basic parameter validation
  if (!is.numeric(n_ages) || length(n_ages) != 1 || n_ages <= 0) {
    stop("n_ages must be a positive numeric value.")
  }
  if (!is.numeric(f) || length(f) != 1) {
    stop("f must be a numeric value.")
  }
  if (!is.numeric(m) || (length(m) != 1 && length(m) != n_ages)) {
    stop("m must be either a single numeric value or a numeric vector of 
         length n_ages.")
  }
  if (!is.numeric(selectivity) || length(selectivity) != n_ages) {
    stop("selectivity must be a numeric vector with length n_ages.")
  }
  
  # If m is given as a constant, turn this into a vector with a value for each age
  if (length(m) == 1) {
    m <- rep(m, 1, n_ages)
  } 
  
  # Age = recruit
  survivorship_at_age <- numeric(n_ages)
  survivorship_at_age[1] <- 1
  # Ages 1 to the plus class (but not including the plus class)
  for (i in 2:(n_ages-1)) {
    survivorship_at_age[i] <- survivorship_at_age[i-1] * exp(-(m[i-1] + f * selectivity[i-1]))
  }
  # Last age class
  final_survivorship <- survivorship_at_age[n_ages-1] * 
    exp(-(m[n_ages-1] + f * selectivity[n_ages-1])) / (1 - exp(-(m[n_ages] + 
                                                                   f * selectivity[n_ages])))
  survivorship_at_age[n_ages] <- final_survivorship
  
  return(survivorship_at_age)
}

# Function returns commonly used equilibrium reference points (including SSB_0, SSB_MSY, 
# F_MSY, and F_X%SPR).
# 
# n_ages = the number of age classes considered
# m = the vector or single numeric value that describes mortality-at-age
# selectivity = a vector describing selectivity-at-age
# waa = a vector describing weight-at-age
# mat = a vector describing maturity-at-age
# srr = either "ricker" to assume a Ricker SRR, "beverton-holt"
# (default) to assume a Beverton-Holt SRR, or NULL to assume no SRR and use a single mean 
# value (see recruitment_mean below)
# recruitment_mean = required if srr is NULL (and ignored if not), this is the single mean 
# value used in place of a SRR
# alpha and beta = numeric Beverton-Holt or Ricker SRR parameter values; ignored 
# if srr is NULL (and thus no SRR is assumed), or if R0 and h are provided
# R0 and h = numeric values Beverton-Holt or Ricker SRR parameter values, ignored
# if srr is NULL (and thus no SRR is assumed)
# SPR = the decimal used to calculate the SPR-based reference point (default = 0.4) 
equilibrium_RPs <- function(n_ages, m, selectivity, waa, 
                                         mat, srr = "beverton-holt", 
                                         recruitment_mean = NULL, 
                                         alpha = NULL, beta = NULL, R0 = NULL, h = NULL, 
                                         SPR = 0.4)
{
  # Basic parameter validation
  if (!is.numeric(n_ages) || n_ages < 1 || n_ages != as.integer(n_ages)) {
    stop("n_ages must be a positive integer.")
  }
  if (!is.numeric(selectivity) || length(selectivity) != n_ages) {
    stop("selectivity must be a numeric vector of length n_ages.")
  }
  if (!is.numeric(selectivity) || length(selectivity) != n_ages) {
    stop("waa must be a numeric vector of length n_ages.")
  }
  if (!is.numeric(selectivity) || length(selectivity) != n_ages) {
    stop("mat must be a numeric vector of length n_ages.")
  }
  if (!is.numeric(m) || (length(m) != 1 && length(m) != n_ages)) {
    stop("m must be either a single numeric value or a numeric vector of 
         length n_ages.")
  }
  if (srr != "ricker" && srr!= "beverton-holt" && !is.null(srr)) {
    stop("srr must either be 'ricker', 'beverton-holt', or NULL.")
  }
  if (!is.null(srr) && ((length(alpha) != 1 || length(beta) != 1) && 
                        (length(R0) != 1 || length(h) != 1))) {
    stop("Either alpha and beta or R0 and h must be provided as single numeric values.")
  }
  if (length(SPR) != 1 || !is.numeric(SPR)) {
    stop("SPR must be a single numeric value.")
  }
  if (is.null(srr) && (length(recruitment_mean) != 1 || 
      !is.numeric(recruitment_mean))) {
    stop("recruitment_mean must be provided as a single numeric value.")
  }
  
  # List to store reference points in
  ref_pts <- list()
  
  # Vector of fishing mortality rates
  f <- seq(0, 5, 0.001)
  
  # Calculate SSB-per-recruit for all values of f
  SSB_per_recruit <- sapply(f, function(fj) sum(survivorship(n_ages = n_ages, m = m, 
                                                             selectivity = selectivity, 
                                                             f = fj) * waa * mat))
  
  # Calculate equilibrium recruitment based on the selected SRR
  if (is.null(srr)) {
    equilibrium_recruitment <- rep(recruitment_mean, length(f))
  } 
  else
  {
    if(!is.null(h)&&!is.null(R0)) # Use h and R0 if given, and if SRR isn't NULL
    {
      alpha <- 4*h/(SSB_per_recruit[1]*(1-h))
      beta <- (alpha*SSB_per_recruit[1]-1)/(R0*SSB_per_recruit[1])
    }
    if(srr == "ricker")
    {
      equilibrium_recruitment <- log(alpha * SSB_per_recruit) / (beta * SSB_per_recruit)
    }
    else # beverton-holt
    {
      equilibrium_recruitment <- (alpha * SSB_per_recruit - 1) / (beta * SSB_per_recruit)
    }
  }
  
  # Calculate yield-per-recruit for all values of f
  yield_per_recruit <- sapply(f, function(fj) 
    sum(survivorship(n_ages = n_ages, m = m, selectivity = selectivity, f = fj) * 
          waa * (1 - exp(-(m + fj * selectivity))) * (fj * selectivity) / (m + fj * selectivity)))
  
  # Calculate equilibrium yield and SPR
  equilibrium_yield <- equilibrium_recruitment * yield_per_recruit
  SPR_val <- SSB_per_recruit / SSB_per_recruit[1]
  
  # Store F_MSY... can handle rare function use case where there is multiple sol'n
  max_index <- which.max(equilibrium_yield)
  if (length(max_index) > 0) {
    ref_pts["F_MSY"] <- f[max_index]
  } else { ref_pts["F_MSY"] <- f[1]}

  # Store F_X%SPR... can handle rare function use case where there is multiple sol'n
  min_index <- which.min(abs(SPR_val - SPR))
  if (length(min_index) > 0) {
    ref_pts["F_X%SPR"] <- f[min_index]
  } else {
    ref_pts["F_X%SPR"] <- f[1]
  }
  
  # Store other calculations in list
  ref_pts["SSB_0"] <- SSB_per_recruit[1]*equilibrium_recruitment[1]
  ref_pts["SSB_MSY"] <- equilibrium_recruitment[which(f == ref_pts["F_MSY"])]*
    SSB_per_recruit[which(f == ref_pts["F_MSY"])]
  ref_pts["MSY"] <- equilibrium_yield[which(f == ref_pts["F_MSY"])]
  ref_pts["F_MSY_SPR"] <- SPR_val[which(f == ref_pts["F_MSY"])]*100 # X that gives F_X%SPR = F_MSY
  ref_pts["SSB_X%SPR"] <- equilibrium_recruitment[which(f == ref_pts["F_X%SPR"])]*
    SSB_per_recruit[which(f == ref_pts["F_X%SPR"])]
  ref_pts["phi_F_MSY"]<- SSB_per_recruit[which(f == ref_pts["F_MSY"])]
  ref_pts["phi_F_X%SPR"]<- SSB_per_recruit[which(f == ref_pts["F_X%SPR"])]
  ref_pts["phi_0"]<- SSB_per_recruit[1]
  
  # If have negative SSB ref pts, make them 0... in SSB_0 case, it means 1/phi0
  # intersects the SRR at a negative SSB
  if(ref_pts["SSB_0"] < 0) {ref_pts["SSB_0"] <- 0}
  if(ref_pts["SSB_MSY"] < 0) {ref_pts["SSB_MSY"] <- 0}
  if(ref_pts["SSB_X%SPR"] < 0) {ref_pts["SSB_X%SPR"] <- 0}
  
  return(ref_pts)
}

# Function calculates commonly used equilibrium reference points (SSB_0, SSB_MSY, F_MSY, 
# F_X%SPR, etc.) by averaging empirical values of m, selectivity, weight-at-age, and 
# maturity-at-age over one or more desired timeframes. If phi_0_at_age is set to TRUE, the 
# function instead calculates phi_0-at-age (again, by averaging over the specified 
# timeframes).
# 
# DATA = a dataframe with a row for each year-age combination, and corresponding 
# estimates of the productivity parameters
# years = the name of the dataframe's year column
# average_yrs_list = a list of vectors, where each vector contains a years to average across
# ages = the name of the dataframe's ages column
# m = a single numeric value for constant mortality-at-age, or the name of the dataframe's 
# mortality-at-age column
# selectivity = the name of the dataframe's selectivity column
# waa = the name of the dataframe's weight-at-age column
# mat = the name of the dataframe's maturity-at-age column
# srr = either "ricker" to assume a Ricker SRR, "beverton-holt" (default) to assume a 
# Beverton-Holt SRR, or NULL to assume no SRR and use a single mean value (see 
# recruitment_mean below); ignored if phi_0_at_age == FALSE
# recruitment_mean = required if srr is NULL and phi_0_at_age == FALSE (and ignored 
# if not), this is the single mean value used in place of a SRR
# alpha and beta = numeric Beverton-Holt or Ricker SRR parameter values; ignored 
# if srr is NULL (and thus no SRR is assumed), or if R0 and h are provided
# R0 and h = numeric values Beverton-Holt or Ricker SRR parameter values, ignored
# if srr is NULL (and thus no SRR is assumed)
# VB_params = In place of waa, a list containing 'l_inf', 'k', 't0', 'a', and 'b' - i.e., the 
# three von Bertalanffy parameters plus a and b from a log10 length-weight relationship
# SPR = the decimal used to calculate the SPR-based reference point (default = 0.4); 
# ignored if phi_0_at_age = TRUE
# phi_0_at_age = If changed to TRUE from its default of FALSE, the function instead 
# returns only phi_0-at-age (again, by averaging over the specified timeframes)
average_time_period_RPs <- function(DATA, years, average_yrs_list, ages, 
                                              m, selectivity, waa=NULL, mat, 
                                              srr = "beverton-holt", recruitment_mean = NULL,
                                              alpha=NULL, beta=NULL, R0=NULL, h=NULL, VB_params=NULL,
                                              SPR = 0.4, phi_0_at_age = FALSE)
{
  # Basic parameter validation
  if (!is.data.frame(DATA)) {
    stop("DATA must be a data frame.")
  }
  if (!is.character(years) || !years %in% names(DATA)) {
    stop("years must be a character string matching a column name in the data frame.")
  }
  if (!is.character(ages) || !ages %in% names(DATA)) {
    stop("ages must be a character string matching a column name in the data frame.")
  }
  if (!is.numeric(m) && (!is.character(m) || !m %in% names(DATA))) {
    stop("m must be either a numeric constant or a character string matching a 
           column name in the data frame.")
  }
  if (!is.character(selectivity) || !selectivity %in% names(DATA)) {
    stop("selectivity must be a character string matching a column name in the data frame.")
  }
  if (!is.character(waa)&&!is.null(waa)) {
    stop("waa must be a character string matching a column name in the data frame.")
  }
  if (!is.character(mat) || !mat %in% names(DATA)) {
    stop("mat must be a character string matching a column name in the data frame.")
  }
  if (phi_0_at_age == FALSE && srr != "ricker" 
      && srr!= "beverton-holt" && !is.null(srr)) {
    stop("srr must be provided as either 'ricker', 'beverton-holt', or NULL.")
  }
  if (!is.null(srr) && phi_0_at_age == FALSE && 
      ((length(alpha) != 1 || length(beta) != 1) && 
      (length(R0) != 1 || length(h) != 1))) {
    stop("Either alpha and beta or R0 and h must be provided as single numeric values.")
  }
  if (length(SPR) != 1 || !is.numeric(SPR)) {
    stop("SPR must be a single numeric value.")
  }
  if (is.null(srr) && phi_0_at_age == FALSE && (length(recruitment_mean) != 1 || 
      !is.numeric(recruitment_mean))) {
    stop("recruitment_mean must be provided as a single numeric value.")
  }
  if (is.null(VB_params) && is.null(waa)) {
    stop("Either 'VB_params' or 'waa' must be provided.")
  }
  if (length(phi_0_at_age)!= 1 || is.logical(phi_0_at_age) == FALSE){
    stop("'phi_0_at_age' must be a single logical value (i.e., TRUE or FALSE)")
  }
  required_params <- c("l_inf", "k", "t0", "a", "b")
  if (!is.null(VB_params)) {
    if (!all(required_params %in% names(VB_params))) {
      stop("The list 'VB_params' must contain 'l_inf', 'k', 't0', 'a', and 'b' parameters.")
    }
    if (!all(sapply(VB_params[required_params], function(x) is.numeric(x) && length(x) == 1))) {
      stop("Parameters in 'VB_params' must be numeric values of length 1.")
    }
  }
  
  # If m is given as a constant, add m column to DATA with all rows as that constant
  if (is.numeric(m))
  {
    DATA$"M_value" <- m
    m <- "M_value"
  }
  
  # If given von Bertalanffy parameters plus a and b from length-weight curve,
  # use that to fit a weight for each age class
  if(is.null(VB_params)==FALSE)
  {
     waa <- "waa"
     DATA[waa] <- rep(10^(log10(VB_params[["a"]])+VB_params[["b"]]*
                            log10(von_bertalanffy(DATA[1:nrow(unique(DATA[ages])), ages], 
                                        VB_params[["l_inf"]], 
                                        VB_params[["k"]], 
                                        VB_params[["t0"]]))), 
                     length(unique(DATA[years]))) 
  }
  
  # Create a list to store the results for each set of average_yrs
  results <- list()
  
  # Loop through average_yrs_list
  for (i in seq_along(average_yrs_list)) {
    average_yrs <- average_yrs_list[[i]]  # Get the current set of years
    
    # Create filter condition for the years in the current set
    filter_condition <- DATA[[years]] %in% average_yrs
    
    # Average m, selectivity, weight-at-age, and maturity-at-age over desired timeframe
    AVERAGE_DATA <- DATA %>% filter(filter_condition) %>% group_by_at({{ ages }}) %>%
      summarize(M_average = mean(.data[[m]]), 
                selectivity_average = mean(.data[[selectivity]]),
                waa_average = mean(.data[[waa]]), 
                mat_average = mean(.data[[mat]]))
    AVERAGE_DATA$selectivity_average <- AVERAGE_DATA$selectivity_average/
      max(AVERAGE_DATA$selectivity_average) # Normalize selectivity to a max of 1
    
    # Act based on whether phi_0_at_age is requested or SSB_0, SSB_MSY, etc. are
    if(phi_0_at_age == TRUE)
    {
      # Do necessary calculations using the averaged values
      survivorship_at_age <- survivorship(n_ages = nrow(unique(DATA[ages])), 
                                          m = AVERAGE_DATA$M_average, 
                                          selectivity = AVERAGE_DATA$selectivity_average, 
                                          f = 0)
      calc <- survivorship_at_age*AVERAGE_DATA$waa_average*AVERAGE_DATA$mat_average
    }
    else
    {
      # Do necessary calculations using the averaged values
      calc <- equilibrium_RPs(n_ages = nrow(unique(DATA[ages])), 
                                           m = AVERAGE_DATA$M_average, 
                                           selectivity = AVERAGE_DATA$selectivity_average, 
                                           waa = AVERAGE_DATA$waa_average, 
                                           mat = AVERAGE_DATA$mat_average,
                                           alpha = alpha, beta = beta, h = h,
                                           R0 = R0, SPR = SPR, 
                                           srr = srr, 
                                           recruitment_mean = recruitment_mean)
    }
    # Store the results in the list
    results[[paste("Time Period", i)]] <- calc
  }
  
  # Convert the list to a dataframe if phi_0_at_age == FALSE
  if(phi_0_at_age == FALSE) {results <- bind_rows(results, .id = "Time Period")}
  
  return(results)
}

# ab_adjust depends on BHnll and is used to "transform" a and b (alpha and beta) to 
# those from the SRR for the recruitment shift
BHnll <- function(theta,S,R){ # a b parameterization. 
  bha <- exp(theta[1])
  bhb <- exp(theta[2])
  bhsig <- exp(theta[3])
  log_R <- log(R)
  P_log_R <- log(bha*S/(1+bhb*S))
  negloglike <- sum(log(bhsig)+(log_R-P_log_R)^2/(2*bhsig^2))
  return (negloglike)
}
ab_adjust <- function(BHa,BHb,meanPerry){
  SRR <- data.frame(SSB=1:1000)
  SRR$R <- mean(BHa)*SRR$SSB/(1+mean(BHb)*SRR$SSB)
  SRR$Rb <- meanPerry*SRR$R#exp(log(meanPerry)+log(SRR$R))  
  mle_BH <- optim(fn=BHnll,par=c(log(mean(BHa)),log(mean(BHb)),log(0.5)),S=SRR$SSB,
                  R=SRR$Rb,hessian = T)
  return(exp(mle_BH$par[1:2]))
}

# Sets parameters for observation model portion of OM
obs_change <- function(OM){
  OM@Cobs <- c(0,0.05); OM@Cbiascv <- 0.01; OM@CAA_nsamp <- c(1000,2000) 
  OM@CAA_ESS <- c(800,1200); OM@CAL_nsamp <- c(1000,2000)
  OM@Iobs <- c(0.1,0.25); OM@Btobs <- c(0,0.05); OM@Btbiascv <- 0.01 
  OM@beta <- c(1,1); LenMbiascv <- 0.0; OM@Mbiascv <- 0.001; OM@Kbiascv <- 0.01 
  OM@t0biascv <- 0.01; OM@Linfbiascv <- 0.01; OM@LFCbiascv <- 0.01; OM@LFSbiascv <- 0.01
  OM@FMSY_Mbiascv <- 0.01; OM@BMSY_B0biascv <- 0.01; OM@Irefbiascv <- 0.01 
  OM@Brefbiascv <- 0.01; OM@Crefbiascv <- 0.01; OM@Dbiascv <- 0.01
  OM@Dobs <- c(0,0.05); OM@hbiascv <- 0.01; OM@Recbiascv <- c(0,0.05) 
  OM@sigmaRbiascv<-0; OM@Eobs <- c(0,0.05); OM@Ebiascv <- 0.001
  return(OM)
}

# Sets projections from the OM to use the average of the last five years for productivity 
# components
cpars_pro5 <- function(OM){
  for(i in 1:100){
    OM@cpars$Wt_age[i,,(OM@nyears+1):dim(OM@cpars$Wt_age)[3]] <- 
      apply(OM@cpars$Wt_age[1,,(OM@nyears - 4):OM@nyears],1,mean)
    OM@cpars$V[i,,(OM@nyears+1):dim(OM@cpars$V)[3]] <- 
      apply(OM@cpars$V[i,,(OM@nyears - 4):OM@nyears],1,mean) / max(apply(OM@cpars$V[i,,(OM@nyears - 4):OM@nyears],1,mean))
    OM@cpars$Mat_age[i,,(OM@nyears+1):dim(OM@cpars$Mat_age)[3]] <- 
      apply(OM@cpars$Mat_age[1,,(OM@nyears - 4):OM@nyears],1,mean) 
    OM@cpars$M_ageArray[i,,(OM@nyears+1):dim(OM@cpars$M_ageArray)[3]] <- 
      apply(OM@cpars$M_ageArray[1,,(OM@nyears - 4):OM@nyears],1,mean)
  }
  return(OM)
}

# Sets projections from the OM to use the flipped historical time series for productivity 
# components
cpars_proflip <- function(OM){
  
  # Function to reverse a row
  reverse_rows <- function(row, constant) {
    rev_row <- rev(row)
    rev_row <- c(rev_row, rep(row[1], 100))
    return(rev_row)
  }
  
  for(i in 1:100){
    OM@cpars$Wt_age[i,,(OM@nyears+1):dim(OM@cpars$Wt_age)[3]] <- 
      t(apply(OM@cpars$Wt_age[i,,1:OM@nyears], 1, reverse_rows))[,1:OM@proyears]
    OM@cpars$V[i,,(OM@nyears+1):dim(OM@cpars$V)[3]] <- 
      t(apply(OM@cpars$V[i,,1:OM@nyears], 1, reverse_rows))[,1:OM@proyears]
    OM@cpars$Mat_age[i,,(OM@nyears+1):dim(OM@cpars$Mat_age)[3]] <- 
      t(apply(OM@cpars$Mat_age[i,,1:OM@nyears], 1, reverse_rows))[,1:OM@proyears] 
    OM@cpars$M_ageArray[i,,(OM@nyears+1):dim(OM@cpars$M_ageArray)[3]] <- 
      t(apply(OM@cpars$M_ageArray[i,,1:OM@nyears], 1, reverse_rows))[,1:OM@proyears]
  }
  return(OM)
}

# Sets OM mat in the projections to gradually shift by a specified annual_shift
# (where a - annual_shift indicates a left shift)
mat_shift_proj <- function(OM, annual_shift){
  
  # Calculate shifted mat ogive for each projection year based on baseline (last historical year)
  mat_data <- as.data.frame(cbind(mat = OM@cpars$Mat_age[1,,OM@nyears], 
                                  age = 0:OM@maxage))
  mat_fit <- glm(mat ~ age, data = mat_data, family = "binomial")
  mat_matrix_proj <- OM@cpars$Mat_age[1,,1:OM@proyears]
  for (i in 1:OM@proyears) {
    mat_matrix_proj[, i] <- plogis(predict(mat_fit, 
                                   newdata = data.frame(age = 0:OM@maxage-
                                                        annual_shift*i)))
  }

  # Populate sims with resulting mat matrix
  for(i in 1:100)
  {
    OM@cpars$Mat_age[i,,(OM@nyears+1):(OM@proyears+OM@nyears)] <- mat_matrix_proj
  }
  
  return(OM)
}

# Sets OM waa in the projections to change exponentially based on a specified annual_rate
waa_trend_proj <- function(OM, annual_rate)
{
  # Calculate changes in waa  
  waa_matrix_proj <- OM@cpars$Wt_age[1,,1:OM@proyears]
  for (i in 1:nrow(waa_matrix_proj)){
    waa_matrix_proj[i, ] <- OM@cpars$Wt_age[1,i,OM@nyears]*(annual_rate)^(1:OM@proyears)
  }
  
  # Populate sims with resulting waa matrix
  for(i in 1:OM@nsim){
    OM@cpars$Wt_age[i,,(OM@nyears+1):(OM@proyears+OM@nyears)] <- waa_matrix_proj
  }
  
  return(OM)
}

# Takes a dataframe of MSE results and appends, to each row, the correct sim's 
# LRPs, USRs, and first 5 years SSB0
sim_LRP_SSB0 <- function(OM, MSE, DF, srr = mysrr)
{
  # Calculate mean phi0 across first 5 yrs
  phi0_first5 <- numeric(5)
  for (i in 1:5) 
  {
    survivorship_unfished <- survivorship(n_ages = length(0:OM@maxage), 
                                          m = OM@cpars$M_ageArray[1,,i], f = 0, 
                                          selectivity = rep(1, length(0:OM@maxage)))
    phi0_first5[i] <- sum(as.vector(OM@cpars$Wt_age[1,,i])*survivorship_unfished
                          *as.vector(OM@cpars$Mat_age[1,,i]))
  }
  mean_phi0_first5 <- mean(phi0_first5)
  
  # Calculate the LRPs and USRs for each sim-projection year combo, and calculate
  # the first 5 years SSB0 for each sim
  RPs <- data.frame()
  for(i in 1:100)
  {
    # Calculate alpha and beta using the mean phi0 from first 5 yrs
    R0 <- OM@cpars$R0[i]
    h <- mean(MSE@RefPoint$ByYear$h[1,1:5])
    if(srr == 'beverton-holt'){
      alpha <- 4*h/(mean_phi0_first5*(1-h))
      beta <- (alpha*mean_phi0_first5-1)/(R0*mean_phi0_first5)
    }
    if(srr == 'ricker'){
      alpha <- ((5*h)^1.25)/(mean_phi0_first5)
      beta <- log(alpha*mean_phi0_first5)/(R0*mean_phi0_first5)
    }
    
    # Calculate the LRPs and USRs for each projection year
    for(y in 1:OM@proyears)
    {
      # Grab necessary data from the OM
      DATA <- data.frame(years = rep(1:(OM@nyears+y), each = OM@maxage+1),
                         ages = rep(0:OM@maxage, OM@nyears+y),
                         waa = as.vector(OM@cpars$Wt_age[1,,1:(OM@nyears+y)]), 
                         mat = as.vector(OM@cpars$Mat_age[1,,1:(OM@nyears+y)]),
                         m = as.vector(OM@cpars$M_ageArray[1,,1:(OM@nyears+y)]),
                         selectivity = as.vector(OM@cpars$V[i,,1:(OM@nyears+y)]))
      
      # Calculate and store the LRPs and USRs plus the first 5 SSB0
      RP_DF <- average_time_period_RPs(DATA = DATA, years = "years", 
                                       ages = "ages", m = "m", 
                                       selectivity = "selectivity", 
                                       waa = "waa", mat = "mat", 
                                       srr = srr, alpha = alpha, beta = beta,
                                       SPR = 0.4, 
                                       average_yrs = list(1:5,
                                                          (OM@nyears+y-4):(OM@nyears+y),
                                                          1:(OM@nyears+y)))
      RP_DF$Sim <- i
      RP_DF$Yr <- y
      RP_DF$`Time Period` <- c("First5","Last5","EntireTS")
      RPs <- rbind(RPs, RP_DF)
    }
  }
  
  # Add correct sim's LRPs, USRs, and SSB0 to DF
  set.seed(OM@seed)
  sim_used <- sample(1:100, OM@nsim)
  for (i in 1:length(sim_used)){
    for(y in 1:OM@proyears)
    {
      DF[DF$Yr == y&DF$Sim == i, c("LRP_First5","USR_First5")] <- 
        RPs$SSB_MSY[RPs$Sim == sim_used[i]&RPs$Yr == y&RPs$`Time Period`== "First5"]/1000
      DF[DF$Yr == y&DF$Sim == i, c("LRP_Last5","USR_Last5")] <- 
        RPs$SSB_MSY[RPs$Sim == sim_used[i]&RPs$Yr == y&RPs$`Time Period`== "Last5"]/1000
      DF[DF$Yr == y&DF$Sim == i, c("LRP_EntireTS","USR_EntireTS")] <- 
        RPs$SSB_MSY[RPs$Sim == sim_used[i]&RPs$Yr == y&RPs$`Time Period`== "EntireTS"]/1000
    }
    DF$ssb0_first5_sim[DF$Sim == i] <- 
      RPs$SSB_0[RPs$Sim == sim_used[i]&RPs$Yr == 1&RPs$`Time Period`== "First5"]/1000
  }
  DF[ , grep("LRP", names(DF))] <- DF[ , grep("LRP", names(DF))]*0.4
  DF[ , grep("USR", names(DF))] <- DF[ , grep("USR", names(DF))]*0.8
  
  return(DF)
}

# A CppAD-friendly alternative to max()
find_max <- function(values)
{
  max_value <- values[1]
  for (i in 2:length(values)) 
  {
    max_value <- CondExpGt(values[i], max_value, values[i], max_value)
  }
  return(max_value)
}

# Logistic or dome-shaped selectivity for a given age... dome if both optional 
# parameters a50_2 and k_2 provided, otherwise logistic with parameters a50 and k
calculate_selectivity <- function(age, a50, k, a50_2 = NULL, k_2 = NULL) {
  if(is.null(a50_2)&is.null(k_2)) {
    return((1/(1+exp((a50-age)/k))))
  } else {
    return((1/(1+exp((a50-age)/k))) * (1/(1+exp(-(a50_2-age)/k_2))))
  }
}