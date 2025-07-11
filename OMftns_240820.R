library(abind); library(MSEtool); library(SAMtool)


# Takes a model fit using the SCA_fit function and draws from the covariance matrix
# to set up an OM
SCA2OM <- function (obj, dat, nsim = OM@nsim, proyears = 50, interval = 1, Name = NULL, 
                      WLa = 1, WLb = 3, Obs = MSEtool::Imprecise_Unbiased, myseed,
                      Imp = MSEtool::Perfect_Imp, nyr_par_mu = 5, LowerTri = 2,
                      plusgroup = T, altinit = 0, fixq1 = T, report = FALSE, silent = FALSE,
                      ...) 
{
  # Simulate from the covariance matrix of the model
  set.seed(myseed)
  SD <- TMB::sdreport(obj, getJointPrecision = TRUE)
  mu <- SD$par.fixed
  covm <- SD$cov.fixed
  numu <- rep(NA, length(mu))
  munams <- names(mu)
  nams <- unique(munams)
  ni <- length(nams)
  covnams <- rownames(covm)
  ind <- 1:length(mu)
  for (i in 1:ni) {
    indto <- ind[covnams == nams[i]]
    indfrom <- ind[munams == nams[i]]
    numu[indto] <- mu[indfrom]
  }
  names(numu) <- covnams
  samps <- mvtnorm::rmvnorm(nsim, numu, covm, checkSymmetry = F)
  report_internal_fn <- function(x, samps, obj) {
    obj$report(samps[x, ])
  }
  output <- lapply(1:nsim, report_internal_fn, samps = samps, obj = obj)
  dims <- dim(output[[1]]$logN)
  ny <- dims[1]
  na <- dims[2]
  
  # Set up naa, recruitment, etc. for each sim 
  naa <- aperm(array(unlist(lapply(output, FUN = function(x) exp(x$logN))), 
                     c(ny, na, nsim)), c(3, 2, 1))
  R0 <- exp(samps[,colnames(samps)=="logR0"])
  h <- 0.8*plogis(samps[,colnames(samps)=="trans_h"])+0.2
  Maa <- aperm(array(unlist(lapply(output, FUN = function(x) dat$M)), 
                     c(ny, na, nsim)), c(3, 2, 1))
  Mataa <- aperm(array(unlist(lapply(output, FUN = function(x) dat$MO)), 
                       c(ny, na, nsim)), c(3, 2, 1))
  waa <- aperm(array(unlist(lapply(output, FUN = function(x) dat$SW)), 
                     c(ny, na, nsim)), c(3, 2, 1))
  laa <- (waa/WLa)^(1/WLb)
  
  # Calculate annual F for each sim
  AFM <- array(NA, dim = c(nsim, nyears))
  for(s in 1:nsim)
  {
    for(y in 1:(nyears))
    {
      AFM[s,y] <- exp(sum(samps[s,colnames(samps) == "logFM"][1:y]))
    }
  }
  
  # Calculate, for each sim, annual FAA using selectivity and 
  # annual F
  a50 <- plogis(samps[,colnames(samps)=="trans_a50"])*max_age
  k <- exp(samps[,colnames(samps)=="trans_k"])
  FAA <- Maa  
    for(s in 1:nrow(AFM)){
      for(y in 1:ncol(AFM)){
        FAA[s,,y] <- 1/(1+exp((a50[s]-c(1:max_age))/k[s]))*AFM[s,y]
      }
    }
  
  # Use sims to set up and return OM
  CurrentYr <- dat$year[length(dat$year)]
  OM <- VPA2OM(Name = stock, proyears, interval, CurrentYr, h = h, myseed = myseed,
                 Obs, Imp, naa, FAA, waa, Mataa, Maa, laa, nyr_par_mu, 
                 LowerTri, recind = 1, plusgroup, altinit = 2, fixq1 = fixq1, 
                 report = report, silent = silent, R0 = R0, spawn_time_frac = 0)
  OM@seed <- myseed
  return(OM)
}

################################################################################
################################################################################
################################################################################
sample_recruitment <- function (Perr_hist, proyears, procsd, AC, seed) 
{
  if (!missing(seed)) 
    set.seed(seed)
  nsim <- nrow(Perr_hist)
  if (length(procsd) == 1) 
    procsd <- rep(procsd, nsim)
  if (length(AC) == 1) 
    AC <- rep(AC, nsim)
  procmu <- -0.5 * procsd^2 * (1 - AC)/sqrt(1 - AC^2)
  Perr_delta <- rnorm(nsim * proyears, procmu, procsd) %>% 
    matrix(nrow = nsim, ncol = proyears)
  Perr_proj <- matrix(NA_real_, nsim, proyears)
  Perr_proj[, 1] <- AC * Perr_hist[, ncol(Perr_hist)] + Perr_delta[,1] * sqrt(1 - AC^2)
  for (y in 2:ncol(Perr_proj)) Perr_proj[, y] <- AC * Perr_proj[,y - 1] + Perr_delta[, y] * sqrt(1 - AC^2)
  return(Perr_proj)
}

LinInterp <- function (x, y, xlev, ascending = FALSE, zeroint = FALSE) {
  if (zeroint) {
    x <- c(0, x)
    y <- c(0, y)
  }
  if (ascending) {
    x_out <- x[1:which.max(x)]
    y_out <- y[1:which.max(x)]
  }
  else {
    x_out <- x
    y_out <- y
  }
  if (any(xlev < min(x_out))) 
    warning("There are xlev values less than min(x).")
  if (any(xlev > max(x_out))) 
    warning("There are xlev values greater than max(x).")
  approx(x_out, y_out, xlev, rule = 2, ties = "ordered")$y
}


split_along_dim <- function (a, n){
  stats::setNames(lapply(split(a, arrayInd(seq_along(a), dim(a))[, 
                                                                 n]), array, dim = dim(a)[-n], dimnames(a)[-n]), dimnames(a)[[n]])
}

calc_survival <- function (x, StockPars, plusgroup = TRUE, inc_spawn_time = FALSE) {
  dd <- dim(StockPars$M_ageArray)
  all_years <- dd[3]
  sapply(1:all_years, calc_survival_yr, x = x, StockPars = StockPars, 
         plusgroup = plusgroup, inc_spawn_time = inc_spawn_time)
}

calc_survival_yr <- function (yr, x, StockPars, plusgroup = TRUE, inc_spawn_time = FALSE) { 
  n_age <- StockPars$n_age
  lst.age <- dim(StockPars$M_ageArray)[2]
  if (inc_spawn_time) {
    spawn_time_frac <- StockPars$spawn_time_frac[x]
  } else {
    spawn_time_frac <- 0
  }
  surv <- (rep(NA, n_age))
  surv[1] <- 1 * exp(-StockPars$M_ageArray[x, 1, yr] * spawn_time_frac)
  for (a in 2:n_age) {
    surv[a] <- surv[a - 1] * exp(-(StockPars$M_ageArray[x, a - 1, yr] * (1 - spawn_time_frac) + StockPars$M_ageArray[x,a, yr] * spawn_time_frac))
  }
  if (plusgroup) 
     surv[lst.age] <- surv[n_age]/(1 - exp(-StockPars$M_ageArray[x, a, yr]))
    if (lst.age < n_age) 
    surv[(lst.age + 1):n_age] <- 0
  surv
}

CalcUnfishedRefs <- function (x, ageM, N0_a, SSN0_a, SSB0_a, B0_a, VB0_a, SSBpRa, SSB0a_a){
  avg.ind <- 1:(ceiling(ageM[x, 1]) + 1)
  nyears <- dim(N0_a)[2]
  if (length(avg.ind) > nyears) 
    avg.ind <- 1:nyears
  N0 <- mean(N0_a[x, avg.ind])
  SSN0 <- mean(SSN0_a[x, avg.ind])
  SSB0 <- mean(SSB0_a[x, avg.ind])
  B0 <- mean(B0_a[x, avg.ind])
  VB0 <- mean(VB0_a[x, avg.ind])
  SSBpR <- mean(SSBpRa[x, avg.ind])
  if (length(avg.ind) > 1) {
    SSB0a <- apply(SSB0a_a[x, avg.ind, ], 2, mean)
  }
  else {
    SSB0a <- SSB0a_a[x, avg.ind, ]
  }
  list(N0 = N0, SSN0 = SSN0, SSB0 = SSB0, B0 = B0, VB0 = VB0, 
       SSBpR = SSBpR, SSB0a = SSB0a)
}


CalcMSYRefs <- function (x, MSY_y, FMSY_y, SSBMSY_y, BMSY_y, VBMSY_y, ageM, nyears){
  n.yrs <- ceiling(ageM[x, nyears])
  nyears1 <- dim(ageM)[2]
  minY <- floor(n.yrs/2)
  maxY <- n.yrs - minY - 1
  avg.ind <- (nyears - minY):(nyears + maxY)
  avg.ind <- avg.ind[avg.ind > 0]
  if (max(avg.ind) > nyears1) 
    avg.ind <- avg.ind[avg.ind < nyears1]
  MSY <- mean(MSY_y[x, avg.ind])
  FMSY <- mean(FMSY_y[x, avg.ind])
  SSBMSY <- mean(SSBMSY_y[x, avg.ind])
  BMSY <- mean(BMSY_y[x, avg.ind])
  VBMSY <- mean(VBMSY_y[x, avg.ind])
  data.frame(MSY = MSY, FMSY = FMSY, SSBMSY = SSBMSY, BMSY = BMSY, 
             VBMSY = VBMSY)
}

popdynCPP<- function (nareas, maxage, Ncurr, pyears, M_age, Asize_c, MatAge, 
          WtAge, FecAge, Vuln, Retc, Prec, movc, SRrelc, Effind, Spat_targc, 
          hc, R0c, SSBpRc, aRc, bRc, Qc, Fapic, maxF, MPA, control, 
          SSB0c, SRRfun, SRRpars, plusgroup = 0L, spawn_time_frac = 0) {
  .Call("_MSEtool_popdynCPP", PACKAGE = "MSEtool", nareas, 
        maxage, Ncurr, pyears, M_age, Asize_c, MatAge, WtAge, 
        FecAge, Vuln, Retc, Prec, movc, SRrelc, Effind, Spat_targc, 
        hc, R0c, SSBpRc, aRc, bRc, Qc, Fapic, maxF, MPA, control, 
        SSB0c, SRRfun, SRRpars, plusgroup, spawn_time_frac)
}

per_recruit_F_calc <- function (x, yr.ind = 1, StockPars, V, SPR_target = seq(0.2,0.6, 0.05)) {
  maxage <- StockPars$maxage
  plusgroup <- StockPars$plusgroup
  if (length(yr.ind) == 1) {
    M_at_Age <- StockPars$M_ageArray[x, , yr.ind]
    Wt_at_Age <- StockPars$Wt_age[x, , yr.ind]
    Mat_at_Age <- StockPars$Mat_age[x, , yr.ind]
    Fec_at_Age <- StockPars$Fec_Age[x, , yr.ind]
    V_at_Age <- V[x, , yr.ind]
  }
  else {
    M_at_Age <- apply(StockPars$M_ageArray[x, , yr.ind], 
                      1, mean)
    Wt_at_Age <- apply(StockPars$Wt_age[x, , yr.ind], 1, 
                       mean)
    Mat_at_Age <- apply(StockPars$Mat_age[x, , yr.ind], 1, 
                        mean)
    Fec_at_Age <- apply(StockPars$Fec_Age[x, , yr.ind], 1, 
                        mean)
    V_at_Age <- apply(V[x, , yr.ind], 1, mean)
  }
  if (max(which(M_at_Age != 0)) != (maxage + 1)) {
    ind <- which(M_at_Age > 0)
    M_at_Age <- M_at_Age[ind]
    Wt_at_Age <- Wt_at_Age[ind]
    Mat_at_Age <- Mat_at_Age[ind]
    Fec_at_Age <- Fec_at_Age[ind]
    V_at_Age <- V_at_Age[ind]
    maxage <- length(ind) - 1
  }
  boundsF <- c(0.001, 3)
  F_search <- exp(seq(log(min(boundsF)), log(max(boundsF)), 
                      length.out = 50))
  Ref_search <- Ref_int_cpp(F_search, M_at_Age = M_at_Age, 
                            Wt_at_Age = Wt_at_Age, Mat_at_Age = Mat_at_Age, Fec_at_Age = Fec_at_Age, 
                            V_at_Age = V_at_Age, StockPars$relRfun, StockPars$SRRpars[[x]], 
                            maxage = maxage, plusgroup = plusgroup, spawn_time_frac = StockPars$spawn_time_frac[x])
  YPR_search <- Ref_search[1, ]
  SPR_search <- Ref_search[2, ]
  RPS <- Ref_search[3, ]
  FSPR <- vapply(SPR_target, function(xx) LinInterp_cpp(SPR_search,F_search, xlev = xx), numeric(1))
  dYPR_dF <- (YPR_search[-1] - YPR_search[-length(YPR_search)])/(F_search[-1] - F_search[-length(F_search)])
  F01 <- LinInterp_cpp(dYPR_dF, F_search[-length(YPR_search)], xlev = 0.1 * dYPR_dF[1])
  Fmax <- F_search[which.max(YPR_search)]
  SSB <- apply(StockPars$SSB[x, , , ], 2, sum)
  R <- apply(StockPars$N[x, 1, , ], 1, sum)
  Fmed <- LinInterp_cpp(RPS, F_search, xlev = median(R/SSB))
  if (StockPars$SRrel[x] == 3) {
    if (!is.null(StockPars$SPRcrashfun)) {
      if (!is.null(formals(StockPars$SPRcrashfun))) {
        SPRcrash <- StockPars$SPRcrashfun(SSBpR0 = StockPars$SSBpR[x, 1], StockPars$SRRpars[[x]])
        Fcrash <- LinInterp_cpp(SPR_search, F_search, xlev = SPRcrash)
      }
      else {
        SPRcrash <- Fcrash <- NA
      }
    }
    else {
      SPRcrash <- Fcrash <- NA
    }
    return(list(FSPR, FYPR = c(YPR_F01 = F01, YPR_Fmax = Fmax, 
                               SPRcrash = SPRcrash, Fcrash = Fcrash, Fmed = Fmed)))
  }
  if (StockPars$SRrel[x] == 1) {
    CR <- 4 * StockPars$hs[x]/(1 - StockPars$hs[x])
  }
  else if (StockPars$SRrel[x] == 2) {
    CR <- (5 * StockPars$hs[x])^1.25
  }
  alpha <- CR/StockPars$SSBpR[x, 1]
  if (min(RPS) >= alpha) {
    SPRcrash <- min(1, RPS[1]/alpha)
    Fcrash <- 0
  }
  else if (max(RPS) <= alpha) {
    SPRcrash <- local({
      slope <- (SPR_search[length(SPR_search)] - SPR_search[length(SPR_search) - 
                                                              1])/(RPS[length(SPR_search)] - RPS[length(SPR_search) - 
                                                                                                   1])
      out <- SPR_search[length(SPR_search)] - slope * (RPS[length(SPR_search)] - 
                                                         alpha)
      max(out, 0.01)
    })
    Fcrash <- max(F_search)
  }
  else {
    SPRcrash <- LinInterp_cpp(RPS, SPR_search, xlev = alpha)
    Fcrash <- LinInterp_cpp(RPS, F_search, xlev = alpha)
  }
  return(list(FSPR, FYPR = c(YPR_F01 = F01, YPR_Fmax = Fmax, 
                             SPRcrash = SPRcrash, Fcrash = Fcrash, Fmed = Fmed)))
}

Ref_int_cpp <- function (F_search, M_at_Age, Wt_at_Age, Mat_at_Age, Fec_at_Age, 
          V_at_Age, relRfun, SRRpars, maxage, plusgroup = 1L, spawn_time_frac = 0) {
  .Call("_MSEtool_Ref_int_cpp", PACKAGE = "MSEtool", F_search, 
        M_at_Age, Wt_at_Age, Mat_at_Age, Fec_at_Age, V_at_Age, 
        relRfun, SRRpars, maxage, plusgroup, spawn_time_frac)
}

LinInterp_cpp <- function (x, y, xlev) {
  .Call("_MSEtool_LinInterp_cpp", PACKAGE = "MSEtool", x, y, 
        xlev)
}

CalcSPReq <- function (FM, StockPars, n_age, nareas, nyears, proyears, nsim, Hist = FALSE) {
  if (Hist) {
    yind <- 1:nyears
  }
  else {
    yind <- 1:proyears + nyears
  }
  n_y <- length(yind)
  M <- replicate(nareas, StockPars$M_ageArray[, , yind])
  Wt_age <- replicate(nareas, StockPars$Wt_age[, , yind])
  Mat_age <- replicate(nareas, StockPars$Mat_age[, , yind])
  Fec_age <- replicate(nareas, StockPars$Fec_Age[, , yind])
  spawn_time_frac <- StockPars$spawn_time_frac
  spawn_time_frac <- replicate(StockPars$maxage + 1, spawn_time_frac)
  spawn_time_frac <- replicate(n_y, spawn_time_frac)
  spawn_time_frac <- replicate(nareas, spawn_time_frac)
  ind <- which(M[1, , 1, 1] == 0)
  if (length(ind) > 0) {
    ind2 <- which(M[1, , 1, 1] != 0)
    M <- M[, ind2, , ]
    Wt_age <- Wt_age[, ind2, , ]
    Mat_age <- Mat_age[, ind2, , ]
    Fec_age <- Fec_age[, ind2, , ]
    FM <- FM[, ind2, , ]
    n_age <- dim(M)[2]
  }
  Z <- FM + M
  initdist <- replicate(n_y, StockPars$initdist[, 1, ]) %>% aperm(c(1, 3, 2))
  NPR_M <- NPR_F <- array(NA_real_, c(nsim, n_age, n_y, nareas))
  NPR_M[, 1, , ] <- NPR_F[, 1, , ] <- initdist
  NPR_M[, 1, , ] <- NPR_M[, 1, , ] * exp(-M[, 1, , ] * spawn_time_frac[, 1, , ])
  NPR_F[, 1, , ] <- NPR_F[, 1, , ] * exp(-M[, 1, , ] * spawn_time_frac[, 1, , ])
  for (a in 2:n_age) {
    NPR_M[, a, , ] <- NPR_M[, a - 1, , ] * exp(-((M[, a - 1, , ] * (1 - spawn_time_frac[, a, , ])) + M[, a,  , ] * spawn_time_frac[, a, , ]))
    NPR_F[, a, , ] <- NPR_F[, a - 1, , ] * exp(-((Z[, a - 1, , ] * (1 - spawn_time_frac[, a, , ])) + Z[, a,  , ] * spawn_time_frac[, a, , ]))
  }
  if (StockPars$plusgroup) {
    NPR_M[, n_age, , ] <- NPR_M[, n_age, , ]/(1 - exp(-M[, n_age, , ]))
    NPR_F[, n_age, , ] <- NPR_F[, n_age, , ]/(1 - exp(-Z[, n_age, , ]))
  }
  SPReq <- apply(NPR_F * Fec_age, c(1, 3), sum)/apply(NPR_M * Fec_age, c(1, 3), sum)
  return(SPReq)
}


getBlow <- function (x, N, Asize, SSBMSY, SSBpR, MPA, SSB0, nareas, retA, 
          MGThorizon, Find, Perr, M_ageArray, hs, Mat_age, Wt_age, 
          Fec_age, R0a, V, nyears, maxage, mov, Spat_targ, SRrel, aR, 
          bR, Bfrac = 0.5, maxF, ploty = F, plusgroup = 0, SRRfun, 
          SRRpars, spawn_time_frac) {
  Blow <- 0
  return(Blow)
}

Blow_opt <- function (lnq, N, Asize_c, SSBMSYc, SSBpRc, MPA, SSB0c, nareas, 
          retAc, MGThorizonc, Fc, Perrc, Mc, hc, Mac, Wac, Fecac, R0c, 
          Vc, nyears, maxage, movc, Spat_targc, SRrelc, aRc, bRc, Bfrac, 
          maxF, mode = 1, plusgroup = 0, SRRfun, SRRpars, spawn_time_frac = 0) {
 
    return(0)
 
}
optYield <- function (logFa, Asize_c, nareas, maxage, Ncurr, pyears, M_age, 
          MatAge, WtAge, FecAge, WtAgeC, Vuln, Retc, Prec, movc, SRrelc, 
          Effind, Spat_targc, hc, R0c, SSBpRc, aRc, bRc, Qc, MPA, maxF, 
          SSB0c, plusgroup = 0, SRRfun, SRRpars, spawn_time_frac = 0) {
  FMSYc <- exp(logFa)
  simpop <- popdynCPP(nareas, maxage, Ncurr, pyears, M_age, 
                      Asize_c, MatAge, WtAge, FecAge, Vuln, Retc, Prec, movc, 
                      SRrelc, Effind, Spat_targc, hc, R0c, SSBpRc, aRc, bRc, 
                      Qc = 0, Fapic = FMSYc, MPA = MPA, maxF = maxF, control = 2, 
                      SSB0c = SSB0c, SRRfun = SRRfun, SRRpars = SRRpars, plusgroup = plusgroup, 
                      spawn_time_frac = spawn_time_frac)
  Cn <- simpop[[6]]/simpop[[8]] * simpop[[1]] * (1 - exp(-simpop[[8]]))
  Cb <- Cn[, pyears, ] * array(WtAgeC[, pyears], dim = dim(Cn[, pyears, ]))
  -mean(apply(Cb, 2, sum, na.rm = TRUE))
}

calcRefYield <- function (x, StockPars, FleetPars, pyears, Ncurr, nyears, proyears) {
  opt <- optimize(optYield, log(c(0.001, 10)), Asize_c = StockPars$Asize[x, 
  ], StockPars$nareas, StockPars$maxage, Ncurr = Ncurr[x, 
                                                       , ], pyears = pyears, M_age = StockPars$M_ageArray[x, 
                                                                                                          , (nyears):(nyears + proyears)], MatAge = StockPars$Mat_age[x, 
                                                                                                                                                                      , (nyears):(nyears + proyears)], WtAge = StockPars$Wt_age[x, 
                                                                                                                                                                                                                                , (nyears):(nyears + proyears)], FecAge = StockPars$Fec_Age[x, 
                                                                                                                                                                                                                                                                                            , (nyears):(nyears + proyears)], WtAgeC = FleetPars$Wt_age_C[x, 
                                                                                                                                                                                                                                                                                                                                                         , (nyears):(nyears + proyears)], Vuln = FleetPars$V_real[x, 
                                                                                                                                                                                                                                                                                                                                                                                                                  , (nyears):(nyears + proyears)], Retc = FleetPars$retA_real[x, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                              , (nyears):(nyears + proyears)], Prec = StockPars$Perr_y[x, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       (nyears):(nyears + proyears + StockPars$maxage)], movc = split_along_dim(StockPars$mov[x, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              , , , (nyears):(nyears + proyears)], 4), SRrelc = StockPars$SRrel[x], 
  Effind = FleetPars$Find[x, ], Spat_targc = FleetPars$Spat_targ[x], 
  hc = StockPars$hs[x], R0c = StockPars$R0a[x, ], SSBpRc = StockPars$SSBpR[x, 
  ], aRc = StockPars$aR[x, ], bRc = StockPars$bR[x, 
  ], MPA = FleetPars$MPA, maxF = StockPars$maxF, SSB0c = StockPars$SSB0[x], 
  plusgroup = StockPars$plusgroup, SRRfun = StockPars$SRRfun, 
  SRRpars = StockPars$SRRpars[[x]], spawn_time_frac = StockPars$spawn_time_frac[x])
  -opt$objective
}

CalculateQ <- function (x, StockPars, FleetPars, pyears, bounds = c(1e-05, 
                                                      15), control) {
  opt <- optimize(optQ, log(bounds), depc = StockPars$D[x], 
                  SSB0c = StockPars$SSB0[x], StockPars$nareas, StockPars$maxage, 
                  Ncurr = StockPars$N[x, , 1, ], pyears, M_age = StockPars$M_ageArray[x, 
                          , ], MatAge = StockPars$Mat_age[x, , ], Asize_c = StockPars$Asize[x, 
                          ], WtAge = StockPars$Wt_age[x, , ], FecAge = StockPars$Fec_Age[x, 
                          , ], Vuln = FleetPars$V_real[x, , ], Retc = FleetPars$retA_real[x, 
                          , ], Prec = StockPars$Perr_y[x, ], movc = split_along_dim(StockPars$mov[x, 
                          , , , ], 4), SRrelc = StockPars$SRrel[x], Effind = FleetPars$Find[x, 
                          ], Spat_targc = FleetPars$Spat_targ[x], hc = StockPars$hs[x], 
                  R0c = StockPars$R0a[x, ], SSBpRc = StockPars$SSBpR[x, 
                  ], aRc = StockPars$aR[x, ], bRc = StockPars$bR[x, 
                  ], maxF = StockPars$maxF, MPA = FleetPars$MPA, plusgroup = StockPars$plusgroup, 
                  StockPars$VB0[x], SBMSYc = StockPars$SSBMSY[x], SRRfun = StockPars$SRRfun, 
                  SRRpars = StockPars$SRRpars[[x]], control, spawn_time_frac = StockPars$spawn_time_frac[x])
  return(exp(opt$minimum))
}

message_warn<-function (...) {
  if (requireNamespace("usethis", quietly = TRUE)) {
    x <- paste(...)
    usethis::ui_warn(x)
  }
  else if (requireNamespace("crayon", quietly = TRUE)) {
    x <- paste(...)
    return(base::message(crayon::red(paste(base::strwrap(x), 
                                           collapes = "\n"))))
  }
  else {
    return(base::message(...))
  }
}

cp<-function(theOM,oldOM){
  theOM@cpars$Wt_age <- oldOM@cpars$Wt_age[,,1:(theOM@nyears+theOM@proyears)]
  theOM@cpars$M_ageArray <- array(dim=dim(theOM@cpars$Wt_age))
  theOM@cpars$M_ageArray[,,1:theOM@nyears] <- oldOM@cpars$M_ageArray[,,1:theOM@nyears]
  theOM@cpars$Len_age <- oldOM@cpars$Len_age[,,1:(theOM@nyears+theOM@proyears)]
  theOM@cpars$Mat_age <- oldOM@cpars$Mat_age[,,1:(theOM@nyears+theOM@proyears)]
  theOM@cpars$V <- oldOM@cpars$V[,,1:(theOM@nyears+theOM@proyears)]
  theOM@cpars$Perr_y <- oldOM@cpars$Perr_y[,1:(theOM@maxage+theOM@nyears+theOM@proyears)]
  theOM@cpars$Wt_age_C <- oldOM@cpars$Wt_age_C[,,1:(theOM@nyears+theOM@proyears)]
  return(theOM)
}


################################################################################
################################################################################
################################################################################

MAKEmeOM <- function (Name = "Me stock", proyears = 50, interval = 1, myseed=666,
                      CurrentYr = as.numeric(format(Sys.Date(), "%Y")), h = 0.999, 
                      Obs = MSEtool::Imprecise_Unbiased, Imp = MSEtool::Perfect_Imp, 
                      naa, faa, waa, Mataa, Maa, laa, nyr_par_mu = 3, LowerTri = 1, 
                      recind = 1, plusgroup = TRUE, altinit = 0, fixq1 = TRUE, R0=R0,
                      report = T, silent = FALSE, ...) 
{
  set.seed(myseed)
  simup <- function(param, OM, do_slot = TRUE) {
    paramnam <- deparse(substitute(param))
    if (length(param) == 1) {
      OM@cpars[[paramnam]] <- rep(param, nsim)
    }
    else if (length(param) < nsim) {
      stop(paste("parameter vector", paramnam, "was length", 
                 length(param), "which is neither length 1 or nsim long"))
    }
    else {
      OM@cpars[[paramnam]] <- param
    }
    if (do_slot) 
      slot(OM, paramnam) <- rep(param[1], 2)
    OM
  }
  dots <- list(...)
  cond <- vapply(list(faa, waa, Mataa, Maa, laa), function(x) all(dim(naa) == dim(x)), logical(1))
  if (!length(cond) || any(is.na(cond)) || any(!cond)) {
    stop("One or more of the following arrays do not have the same shape: naa, faa, waa, Mataa, Maa, Laa")
  }
  if (!is.null(dots$fecaa)) {
    fecaa <- dots$fecaa
    if (!all(dim(naa) == dim(fecaa))) 
      stop("Dimension of fecaa not equal to that for naa")
  }
  else {
    fecaa <- waa * Mataa
  }
  nyears <- dim(naa)[3]
  nsim <- dim(naa)[1]
  if (recind == 1) {
    ageind <- 1:dim(naa)[2]
    dims <- c(nsim, 1, nyears)
    zeros <- array(0, dims)
    rec <- naa[, 1, 2:nyears]
    muR0 <- apply(rec, 1, mean)
    N0 <- array(cbind(rec, muR0), dims)
    naa <- abind(N0, naa, along = 2)
    faa <- abind(zeros, faa, along = 2)
    waa <- abind(zeros, waa, along = 2)
    Mataa <- abind(zeros, Mataa, along = 2)
    laa <- abind(zeros, laa, along = 2)
    fecaa <- abind(zeros, fecaa, along = 2)
    Maa <- abind(array(0.000001 + .Machine$double.eps, dim(Maa[,1, ])), Maa, along = 2)
    message("Age zero positions for arrays were created with the following assumptions: N(0) = N(1) * exp(M(1)), N0 in most recent year is mean(R0), F(0) = weight(0) = maturity(0) = length(0) = 0, M(0) = M(1)")
  }
  n_age <- dim(naa)[2]
  maxage <- dim(naa)[2] - 1
  if (!is.null(dots$SRrel)) {
    SRrel <- dots$SRrel
  }
  else {
    SRrel <- 1
  }
  if (!is.null(dots$R0)) {
    if (length(dots$R0) == 1) {
      R0 <- rep(dots$R0, nsim)
    }
    else {
      R0 <- dots$R0
    }
  }
  else {
    R0 <- naa[, 1, 1]
  }
  if (length(h) == 1) 
    h <- rep(h, nsim)
  OM <- new("OM")
  if (!is.null(dots$spawn_time_frac)) {
    if (lengths(dots$spawn_time_frac) == 1) {
      spawn_time_frac <- rep(dots$spawn_time_frac, nsim)
    }
    else {
      spawn_time_frac <- dots$spawn_time_frac
    }
  }
  else {
    spawn_time_frac <- rep(0, nsim)
  }
  SBsurv <- lapply(1:nsim, calc_survival, StockPars = list(M_ageArray = Maa,n_age = n_age, spawn_time_frac = spawn_time_frac), plusgroup = plusgroup, 
                   inc_spawn_time = TRUE) %>% abind(along = 3) %>% aperm(c(3, 1, 2))
  SSBpR <- apply(SBsurv * fecaa, c(1, 3), sum)
  SSBpR_out <- vapply(1:nsim, function(x) {
    ageM <- min(MSEtool:::LinInterp(Mataa[x, , 1], y = 1:n_age, 0.5), 
                maxage - 1)
    SSBpR[x, 1:(ceiling(ageM) + 1)] %>% mean()
  }, numeric(1))
  if (!is.null(dots$phi0)) {
    new_SR <- local({
      if (length(dots$phi0) == 1) {
        phi0 <- rep(dots$phi0, nsim)
      }
      else {
        phi0 <- dots$phi0
      }
      if (SRrel == 1) {
        Arec <- 4 * h/(1 - h)/phi0
        Brec <- (5 * h - 1)/(1 - h)/R0/phi0
        K <- Arec * SSBpR_out
        h_out <- K/(4 + K)
        R0_out <- (5 * h_out - 1)/(1 - h_out)/Brec/SSBpR_out
      }
      else {
        Arec <- (5 * h)^1.25/phi0
        Brec <- 1.25 * log(5 * h)/R0/phi0
        h_out <- 0.2 * (Arec * SSBpR_out)^0.8
        R0_out <- 1.25 * log(5 * h_out)/Brec/SSBpR_out
      }
      list(h = h_out, R0 = R0_out)
    })
    R0 <- new_SR$R0
    h <- new_SR$h
  }
  SSB0 <- R0 * SSBpR_out
  SSB <- apply(naa * exp(-spawn_time_frac * (Maa + faa)) * fecaa, c(1, 3), sum, na.rm = TRUE)
  D <- SSB[, nyears]/SSB0
  OM@Name <- "4x5YHaddock"
  OM@M <- OM@L50 <- OM@L50_95 <- OM@L5 <- OM@LFS <- OM@Vmaxlen <- c(1,1)
  OM@a <- OM@b <- 1
  OM@nsim <- nsim
  OM@nyears <- nyears
  OM@maxage <- maxage
  OM@proyears <- proyears
  OM@interval <- interval
  OM <- simup(D, OM)
  hs <- h
  OM <- simup(hs, OM, do_slot = FALSE)
  OM@h <- rep(hs[1], 2)
  OM <- simup(R0, OM)
  OM@Size_area_1 <- OM@Frac_area_1 <- OM@Prob_staying <- rep(0.5,2)
  LenCV = 0.1
  OM <- simup(LenCV, OM)
  OM@SRrel <- SRrel
  OM@isRel <- FALSE
  OM@CurrentYr <- CurrentYr
  OM@Msd <- OM@Ksd <- OM@Linfsd <- OM@qcv <- OM@Esd <- OM@AC <- rep(0,2)
  OM@pstar <- 0.5
  OM@reps <- 1
  OM@EffYears <- c(1, 5, OM@nyears)
  OM@EffLower <- rep(1, 3)
  OM@EffUpper <- rep(1.1, 3)
  OM@Spat_targ <- c(1, 1)
  OM@qinc <- rep(0, 2)
  OM <- suppressMessages(Replace(OM, Obs, Sub = "Obs"))
  OM <- suppressMessages(Replace(OM, Imp, Sub = "Imp"))
  Wt_age <- M_ageArray <- Len_age <- Mat_age <- V <- Fec_age <- array(NA,c(nsim, n_age, nyears + proyears))
  Find <- apply(faa, c(1, 3), max)
  Wt_age[, , 1:nyears] <- waa
  M_ageArray[, , 1:nyears] <- Maa
  Len_age[, , 1:nyears] <- laa
  Mat_age[, , 1:nyears] <- Mataa
  Fmax <- aperm(array(rep(Find, n_age), c(nsim, nyears, n_age)), c(1, 3, 2))
  V[, , 1:nyears] <- faa/Fmax
  adjustV <- function(Vi) {
    nan.ind <- apply(Vi, 2, max) %>% is.nan() %>% which()
    if (length(nan.ind) > 0) {
      non.nan.ind <- max(nan.ind) + 1
      Vi[, nan.ind] <- Vi[, non.nan.ind]
    }
    Vi
  }
  V <- sapply(1:nsim, function(i) adjustV(V[i, , ]), simplify = "array") %>% 
    aperm(c(3, 1, 2))
  parmu <- function(arr, nyears, proyears, nyr_par_mu) {
    arr[, , nyears + (1:proyears)] <- array(rep(apply(arr[,, nyears - (0:(nyr_par_mu - 1))], 1:2, mean), proyears), c(dim(arr)[1:2], proyears))
    arr
  }
  Wt_age <- parmu(Wt_age, nyears, proyears, nyr_par_mu)
  M_ageArray <- parmu(M_ageArray, nyears, proyears, nyr_par_mu)
  Len_age <- parmu(Len_age, nyears, proyears, nyr_par_mu)
  Mat_age <- parmu(Mat_age, nyears, proyears, nyr_par_mu)
  V <- parmu(V, nyears, proyears, nyr_par_mu)
  OM@cpars$Wt_age <- Wt_age
  OM@cpars$M_ageArray <- M_ageArray
  OM@cpars$Len_age <- Len_age
  OM@cpars$Mat_age <- Mat_age
  OM@cpars$V <- V
  OM@cpars$Find <- Find
  if (!is.null(dots$fecaa)) {
    Fec_age[, , 1:nyears] <- fecaa
    OM@cpars$Fec_age <- parmu(Fec_age, nyears, proyears, nyr_par_mu)
  }
  if (SRrel == 1) {
    recd <- vapply(1:nsim, function(x) (0.8 * R0[x] * h[x] * SSB[x, ])/(0.2 * SSB0[x] * (1 - h[x]) + (h[x] - 0.2) * SSB[x, ]), numeric(nyears)) %>% t()
  }
  else {
    recd <- local({
      a <- (5 * h)^1.25/SSBpR_out
      b <- 1.25 * log(5 * h)/SSBpR_out/R0
      vapply(1:nsim, function(x) a[x] * SSB[x, ] * exp(-b[x] * SSB[x, ]), numeric(nyears)) %>% t()
    })
  }
  recdevs <- log(naa[, 1, ]/recd)
  recdevs[is.na(recdevs)] <- 0
  if (!is.null(dots$Perr)) {
    if (length(dots$Perr) == 1) {
      procsd <- rep(dots$Perr, nsim)
    }
    else {
      procsd <- dots$Perr
    }
  }
  else {
    procsd <- apply(recdevs, 1, sd)
  }
  if (!is.null(dots$AC)) {
    if (length(dots$AC) == 1) {
      AC <- rep(dots$AC, nsim)
    }
    else {
      AC <- dots$AC
    }
  }
  else {
    AC <- apply(recdevs, 1, function(x) stats::acf(x, plot = FALSE)$acf[2,1, 1])
  }
  #surv <- lapply(1:nsim, calc_survival, StockPars = list(M_ageArray = Maa, n_age = n_age), plusgroup = plusgroup, inc_spawn_time = FALSE) %>% abind(along = 3) %>% aperm(c(3, 1, 2))
  surv <- lapply(1:nsim, calc_survival, StockPars = list(M_ageArray = Maa, n_age = n_age, spawn_time_frac=rep(0.4,nsim)), plusgroup = plusgroup, inc_spawn_time = T) %>% abind(along = 3) %>% aperm(c(3, 1, 2))
  surv3 <- lapply(1:nsim, calc_survival, StockPars = list(M_ageArray = array(1e-6,dim=dim(Maa)), n_age = n_age, spawn_time_frac=rep(0.4,nsim)), plusgroup = plusgroup, inc_spawn_time = F) %>% abind(along = 3) %>% aperm(c(3, 1, 2))
  Perr <- array(NA_real_, c(nsim, maxage + nyears - LowerTri))
  if (altinit < 2) {
    Perr[, n_age:1] <- log(naa[, , 1]/(R0 * surv[, , 1]))
  }
  if (altinit == 3) {
    Perr[, n_age:1] <- log(naa[, , 1]/(R0 * surv3[, , 1]))
  }
  else if (altinit == 2) {
    Perr[, n_age:2] <- log(naa[, 1:(n_age - 1), 1]/(R0 * SBsurv[, 1:(n_age - 1), 1]))
    survDLMtool <- aperm(exp(-apply(Maa[, , 1], 1, cumsum)), c(2, 1))
    fac <- surv[, n_age, 1]/survDLMtool[, n_age]
    Perr[, 1] <- log(naa[, n_age, 1]/(R0 * surv[, n_age, 1] * fac))
  }
  Perr[, maxage + 2:(nyears - LowerTri)] <- recdevs[, 2:(nyears - LowerTri)]
  Perr_pro <- sample_recruitment(Perr_hist = Perr, proyears = proyears + LowerTri, procsd = procsd, AC = AC)
  OM@cpars$Perr_y <- exp(cbind(Perr, Perr_pro))
  OM@Perr <- rep(mean(procsd), 2)
  OM@AC <- rep(mean(AC), 2)
  OM@cpars$AC <- AC
  OM@cpars$Perr <- procsd
  if (fixq1) 
    OM@cpars$qs <- rep(1, nsim)
  OM@maxF <- ceiling(max(Find))
  OM@Linf <- rep(1, 2)
  OM@K <- rep(0.3, 2)
  OM@t0 <- rep(0, 2)
  OM@DR <- rep(0, 2)
  OM@MPA <- FALSE
  if (any(spawn_time_frac > 0)) 
    OM@cpars$spawn_time_frac <- spawn_time_frac
  if (!plusgroup) 
    OM@cpars$plusgroup <- 0L
  if (report) {
    if (!silent) 
      message("\nRunning historical simulations to compare VPA output and OM conditioning...\n")
    Hist <- runMSE(OM, Hist = TRUE, silent = TRUE)
    nc <- ceiling((maxage + 1)/3)
    nr <- ceiling((maxage + 1)/nc)
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    par(mfrow = c(nr, nc), mai = c(0.4, 0.4, 0.3, 0.05), 
        omi = c(0.25, 0.25, 0.01, 0.01))
    yrs <- CurrentYr - ((nyears - 1):0)
    cols <- rep("black", nyears)
    pch <- rep(1, nyears)
    cols[nyears - (0:(LowerTri - 1))] <- "blue"
    pch[nyears - (0:(LowerTri - 1))] <- 4
    for (a in 0:maxage) {
      N_OM <- apply(Hist@AtAge$Number[1, a + 1, , ], 1, sum)
      ylim = c(0, max(naa[1, a + 1, ], N_OM, na.rm = TRUE) * 1.05)
      if (a == 0) 
        plot(yrs, naa[1, a + 1, ], xlab = "", ylab = "", col = cols, pch = pch, ylim = ylim, yaxs = "i")
      if (a > 0) 
        plot(yrs, naa[1, a + 1, ], xlab = "", ylab = "", ylim = ylim, yaxs = "i")
      lines(yrs, N_OM, col = "green")
      mtext(paste("Age ", a), 3, line = 0.5, cex = 0.9)
      if (a == 0) 
        legend("top", legend = c("Assessment", "OM", "Recr. ignored (LowerTri)"), text.col = c("black", "green", "blue"), bty = "n")
      res <- N_OM - naa[1, a + 1, ]
      plotres <- abs(res) > mean(naa[1, a + 1, ] * 0.025, 
                                 na.rm = TRUE)
      if (any(plotres, na.rm = TRUE)) {
        for (y in 1:nyears) {
          if (!is.na(plotres[y]) && plotres[y]) 
            lines(rep(yrs[y], 2), c(naa[1, a + 1, y], 
                                    N_OM[y]), col = "red")
        }
      }
    }
    mtext("Year", 1, line = 0.3, outer = T)
    mtext("Numbers", 2, line = 0.4, outer = T)
  }
  return(OM)
  
}
