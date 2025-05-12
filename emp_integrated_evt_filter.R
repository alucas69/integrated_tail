library("pracma")
library("progress")



## dynamic_tau_filter() filters the dynamic (1-kappa) quantile:
## [note: so for kappa = 0.1 it is the UPPER tail quantile]
##
## INPUTS:
## =======
## ydata: vector of length T holding the data
## omega: scalar intercept of the transition equation
## beta: scalar loading of lagged tau
## alpha1: scalar loading of exceedance indicator
## alpha2: scalar loading of exceedance size
## f1: initial value of the threshold
## kappa: the right tail quantile percentage (0.1 for a 10% right tail prob)
## steep: if > 0, then the steepness of the sigmoid logistic curve as an approximation of the step function
## sharp: if TRUE, then use the sharp indicator function rather than the smooth sigmoid approximation
##
## OUTPUTS:
## =======
## llik: the tick-loss objective function
## f_out: vector of lenght T holding the dynamic thresholds
##
dynamic_tau_filter <- function(
    ydata, omega, beta, alpha1, alpha2, f1,
    kappa = 0.1, steep = 3, sharp = FALSE) {
  
  ## initialize
  f_out = ydata
  llik = 0
  dimn = length(ydata)
  
  ## check arguments
  if ((!sharp) & (steep <= 0)) error(paste0("*** ERROR IN dynamic_tau_filter(): sigmoid approximation to the indicator is flat or decreasing (steep = ", steep))
  
  ## loop over the filter
  for (i1 in 1:dimn) {
    
    ## store the filtered value
    f_out[i1] = f1
    
    ## update the objective
    PoT = ydata[i1] - f1
    sigma = ifelse(sharp, as.integer(PoT > 0), 1 / (1 + exp(-steep * PoT)))
    llik = llik - kappa * (1 - sigma) * PoT + (1 - kappa) * sigma * PoT
    
    ## update the dynamic threshold
    f1 = omega + beta * f1 -
      alpha1 * ( 
        - (1 - kappa) * sigma 
        + kappa * (1 - sigma)  
        - sigma * (1 - sigma) * PoT * steep
      ) +
      alpha2 * PoT
  }
  
  ## check on integrity of objective value
  if (is.na(llik) | is.infinite(llik)) llik = 1e10
  
  return(list(
    llik = llik / length(ydata),
    f_out = f_out
  ))
}




## par.trafo.tau() transforms the raw parameters of the threshold filter
##
## INPUT:
## =======
## x: vector of length 4 holding the transformed parameters
##
## OUTPUT:
## =======
## vector holding the required parameters of the filter by name (omega, beta, alpha1, alpha2)
par.trafo.tau <- function(x) {
  return(unlist(list(
    omega = x[1],
    beta = 1 / (1 + exp(-x[2])),
    alpha1 = -5 + 10 / (1 + exp(-x[3])),
    alpha2 = -1 + 2 / (1 + exp(-x[4]))
    # alpha2 = 1 / (1 + exp(-x[4]))
  )))
}





## local.tau.filter() a helper function for estimating the (1-kappa) dynamic 
## threshold model
## [note: so for kappa = 0.1 it is the UPPER tail quantile]
##
## INPUTS:
## =======
## x: vector with parameters
## givemore: see under OUTPUTS 
## steep: if > 0, then the steepness of the sigmoid logistic curve as an approximation of the step function
## sharp: if TRUE, then use the sharp indicator function rather than the smooth sigmoid approximation
##
## OUTPUTS:
## =======
## if givemore is FALSE: 
##    only a double, holding the tickloss objective
## if givemore is TRUE: 
##    a list(), holding:
##      llik: double, tick-loss function
##      path: the filtered dynamic thresholds
local.tau.filter = function(ydata, x, kappa, givemore = FALSE, steep = 3, sharp = FALSE) {
  ## initialize the first quantile to the unconditional one
  f1 = quantile(ydata, 1 - kappa) 
  ## transform parameters and compute tickloss
  para = par.trafo.tau(x)

  f_out = dynamic_tau_filter_cpp(
    ydata, omega = para['omega'], beta = para['beta'], 
    alpha1 = para['alpha1'], alpha2 = para['alpha2'], 
    f1, kappa, steep, sharp)
  ## decide what to return
  if (givemore) {
    return(list(llik = f_out$llik, path = f_out$f_out))
  } else return(f_out$llik)
}




## par.trafo.tau() transforms the raw parameters of the threshold filter
##
## INPUT:
## =======
## x: vector of length 2 + (EVT_OPTIONS$FIXED_OMEGA) holding the transformed parameters
##
## OUTPUT:
## =======
## vector holding the required parameters of the filter by name (omega, beta, alpha1, alpha2)
par.trafo.tailindex <- function(x, EVT_OPTIONS) {
  return(unlist(list(
    omega = ifelse(EVT_OPTIONS$FIXED_OMEGA, EVT_OPTIONS$FIXED_OMEGA_VALUE, exp(x[2])),
    beta = 1,
    alpha = 1 / (1 + exp(-x[1]))
  )))
}



## tail.filter() filter the dynamic tail parameter using the integrated filter
## [NB: the UPPER tail shape dynamics are estimated]
##
## INPUTS:
## =======
## POTdata: vector with Peaks Over Thresholds (POT) data
## omega: the intercept of the transition equation
## beta: the loading of the previous tail parameter in the transition equation
## alpha: the loading of the scaled score in the transition equation
## f1: the initial tail parameter value
##
## OUTPUTS:
## =======
## list():
##    llik: the log-likelihood value of the GPD approximation
##    llik.t: the vector of log-likelihood contributions
##    f_out: vector of filtered tail parameters
tail.filter <- function(POTdata, omega, beta, alpha, f1) {
  
  ## initialize
  dimn = length(POTdata)
  f_out = POTdata
  llik = 0
  
  ## run the filter over all observations
  for (i1 in 1:dimn) {
    ## store value
    f_out[i1] = f1
    ## update the parameter
    f1 = omega + beta * f1 + alpha * (log(1+POTdata[i1]) - f1)
  }
  ## compute the likelihood
  llik.t = log(f_out) + (1 + 1/f_out) * log(1+POTdata)
  llik = sum(llik.t)
  
  ## check log-lik for integrity
  if (is.na(llik) | is.infinite(llik)) llik = 1e10
  
  return(list(
    llik = llik,
    llik.t = llik.t,
    f_out = f_out
  ))
}




## compute_std_err() computes standard error of score-driven model
##
## INTPUT:
## ========
## x: vector of parameters at optimum
## POTdata: vector with Peaks Over Thresholds (POT) data
## f1: the initial tail parameter value
##
## OUTPUT:
## ========
## list():
##    Hinv: inverse Hessian based covariance matrix
##    sandwich: sandwich estimator of covariance matrix
##    Hinv.raw: inverse Hessian based covariance matrix of theta (parameter log lik)
##    sandwich.raw: sandwich estimator of covariance matrix of theta (parameter log lik)
compute_std_err <- function(x, POTdata, f1, EVT_OPTIONS) {
  
  ## parameter transformation part for Delta method
  Deltapart = jacobian(par.trafo.tailindex, x, EVT_OPTIONS = EVT_OPTIONS)[c(1,3), ]
  if (EVT_OPTIONS$FIXED_OMEGA) Deltapart = Deltapart[2]

  ## compute jacobian and OPG
  matJ = jacobian(function(x) {
    par_in = par.trafo.tailindex(x, EVT_OPTIONS)
    value.out = tail.filter(POTdata, par_in['omega'], par_in['beta'], par_in['alpha'], f1)
    return(value.out$llik.t)
  }, x)
  matJ = (t(matJ) %*% matJ) / length(POTdata)
  matJ = (matJ + t(matJ)) / 2
  
  ## compute Hessian 
  Hinv = hessian(function(x) {
    par_in = par.trafo.tailindex(x, EVT_OPTIONS)
    value.out = tail.filter(POTdata, par_in['omega'], par_in['beta'], par_in['alpha'], f1)
    return(mean(value.out$llik.t))
  }, x)
  Hinv = pinv(Hinv)
  Hinv.raw = (Hinv + t(Hinv)) / 2

  ## construct sandwich  
  matJ.raw = (Hinv.raw %*% matJ %*% t(Hinv.raw))
  matJ.raw = (matJ.raw + t(matJ.raw)) / 2
  matJ = Deltapart %*% matJ.raw %*% t(Deltapart)
  matJ = (matJ + t(matJ)) / 2
  Hinv = Deltapart %*% Hinv.raw %*% t(Deltapart)
  Hinv = (Hinv + t(Hinv)) / 2
  nms = "alpha"; if (!EVT_OPTIONS$FIXED_OMEGA) nms = c("omega", nms)
  colnames(Hinv) = rownames(Hinv) = colnames(matJ) = rownames(matJ) = nms
  colnames(Hinv.raw) = rownames(Hinv.raw) = colnames(matJ.raw) = rownames(matJ.raw) = paste0("theta", 1:NROW(Hinv))
  return(list(
    Hinv = Hinv / length(POTdata),
    sandwich = matJ / length(POTdata),
    Hinv.raw = Hinv.raw / length(POTdata),
    sandwich.raw = matJ.raw / length(POTdata)
  ))
}




## simulate_filter_bands() computes bands around filter reflecting estimation uncertainty
## 
simulate_filter_bands <- function(x, POTdata, covmatrix, f1, nrsims = 100, pct = 0.5, raw = TRUE, returnsims = FALSE, EVT_OPTIONS) {
  if (nrsims <= 0) {
    bands = paths = NULL
  } else {
    ## initialize
    if (raw) Lt = eigen(covmatrix$sandwich.raw, TRUE) else {
      Lt = eigen(covmatrix$sandwich, TRUE)
      x = par.trafo.tailindex(x, EVT_OPTIONS)[colnames(covmatrix$sandwich)]
    }
    if (NROW(Lt$vectors) == 1) Lt = matrix(sqrt(Lt$values), ncol = 1, nrow = 1) else Lt = Lt$vectors %*% diag(sqrt(Lt$values)) %*% t(Lt$vectors)
    dimk = nrow(Lt)
    dimT = length(POTdata)
    paths = matrix(0, nrow = dimT, ncol = nrsims)
    bands = matrix(0, nrow = dimT, ncol = 2)
    
    ## run simulations
    lifesign = progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                                total = nrsims,
                                complete = "=",   # Completion bar character
                                incomplete = "-", # Incomplete bar character
                                current = ">",    # Current bar character
                                clear = FALSE,    # If TRUE, clears the bar when finish
                                width = 100) 
    i1 = 1; skipped = 0; cat('\nSimulating band: '); while (i1 <= nrsims) {
      if (raw) par_in = par.trafo.tailindex(x + rnorm(dimk) %*% Lt, EVT_OPTIONS) else {
        par_in = x + rnorm(dimk) %*% Lt
        if (EVT_OPTIONS$FIXED_OMEGA) {
          par_in = c(EVT_OPTIONS$FIXED_OMEGA_VALUE, par_in, 1)
        } else {
          par_in = c(par_in, 1)
        }
        names(par_in) = c('omega', 'alpha', 'beta')
      }
      
      ## check correctness of parameters
      if ((par_in['alpha'] > 0) & (par_in['omega'] > 0)) {
        paths[ , i1] = tail.filter(POTdata, par_in['omega'], par_in['beta'], par_in['alpha'], f1)$f_out
        i1 = i1 + 1
        lifesign$tick()
      } else {
        skipped = skipped + c(par_in['omega'] <= 0, par_in['alpha'] <= 0)
      }
    }
    cat('\nDiscarded (omega,alpha) = (', paste0(round(100 * skipped / nrsims), collapse = ' % , '), '%) band simulations\n')
    
    ## compute bands
    for (i1 in sequence(dimT)) bands[i1, ] = quantile(paths[i1, ], probs = c((1 - pct)/2, 1 - (1 - pct) / 2), na.rm = TRUE)
    colnames(bands) = paste0(round(100 * pct), c('-lower', '-upper'))
  }
  
  if (returnsims) return(list(
    bands = bands,
    paths = paths
  )) else return(list(
    bands = bands
  ))
}






## estimate_full_model() produces output for the tau and f.t estimation
## [note: so for kappa = 0.1 it is the UPPER tail quantile]
##
## INPUT:
## =======
## my_data: dataframe()
##    $y: ydata: vector of raw data
##    $dates: vector of corresponding dates
## kappa: upper tail probability for the thresholds for the POTs
## pct: confidence level for the simulated bands for the tail index filter
## nrsims: number of simulations for the confidence bands
## raw.sim: if TRUE, simulate from the parameters used to optimize the log-lik, otherwise simulate omega and alpha directly using the delta method
## smoothings: vector containing the smoothing approximations to the step function (0 is sharp step, -1 is sharp plus Nelder-Mead optimizer)
## EVT_OPTIONS: list()
##      $FIXED_OMEGA: if TRUE, do not estimate omega
##      $FIXED_OMEGA_VALUE: value of omega if $FIXED_OMEGA == TRUE
## verbosity: integer 0..5, higher for more intermediate output;
##
## OUTPUT:
## =======
## list():
##    $my_data: original dataframe plus additional columns
##    $tau.optimizer.output: output from the optim() for the thresholds
##    $tail.optimizer.output: output of optim() for the integrated score-driven tail model
##    $bands: matrix, holding the POT times f.t and the upper and lower confidence bands
##    $covmatrices: list of covariance matrices for the tail index model (.raw for the untransformed parameters, otherwise for the omega and alpha)
estimate_full_model <- function(my_data, smoothings = c(5,20,-1),
                                EVT_OPTIONS = list(),
                                verbosity = 0,
                                lower_tail = FALSE) {
  ## initialize
  tail_multiplier = ifelse(lower_tail, 1, -1)
  tmp_data = tail_multiplier * my_data$y
  pct = ifelse(is.null(EVT_OPTIONS$BANDS_PCT), 0.8, EVT_OPTIONS$BANDS_PCT)
  nrsims = ifelse(is.null(EVT_OPTIONS$BANDS_NRSIMS), 100, EVT_OPTIONS$BANDS_NRSIMS)
  raw.sim = ifelse(is.null(EVT_OPTIONS$BANDS_RAW), FALSE, EVT_OPTIONS$BANDS_RAW)
  if (is.null(EVT_OPTIONS$TAU_TAIL_PCT)) {kappa = 0.05; EVT_OPTIONS$TAU_TAIL_PCT = kappa} else {kappa = EVT_OPTIONS$TAU_TAIL_PCT}
  if (is.null(EVT_OPTIONS$FIXED_OMEGA)) EVT_OPTIONS$FIXED_OMEGA = FALSE
  if (is.null(EVT_OPTIONS$FIXED_OMEGA_VALUE)) EVT_OPTIONS$FIXED_OMEGA_VALUE = 0
  externaltau = EVT_OPTIONS$EXTERNAL_TAU
  
  dimT = nrow(my_data)
  cat('\n#########################\n## EVT model optimization\n#########################\n')
  
  ## loop over different smoothness sigmoids of PZC LOWER (!) tail threshold estimations
  if (is.null(externaltau)) {
    cat("\n\n... NOW ESTIMATING THE DYNAMIC THRESHOLDS\n")
    par0 = c( 0.0003117466,  5.7797848912, 0.0079226142, 0.0057528009) ## warm for pos/neg alpha parameterization
    if (verbosity >= 1) cat('\nInitial parameter:\n', unlist(par.trafo.tau(par0)))
    for (steep in smoothings) {
      if (steep < 0) {
        steep = max(smoothings)
        sharp = TRUE
        my_report = list(maxit = 2000)
        if (verbosity >= 4) {my_report$trace = 1; my_report$REPORT = 1}
        o.out = optim(
          par0, function(x) {return(local.tau.filter(tmp_data, x, kappa = kappa, sharp = TRUE))},
          method = "Nelder-Mead", control = my_report)
        steep = "indicator"
      } else {
        my_report = list(maxit = 500)
        if (verbosity >= 4) {my_report$trace = 1; my_report$REPORT = 1} else {
          if (verbosity >= 3) {my_report$trace = 1; my_report$REPORT = 10}}
        o.out = optim(
          par0, function(x) {return(local.tau.filter(tmp_data, x, kappa = kappa, sharp = FALSE, steep = steep))},
          method = "BFGS", control = my_report)
      }
      if (verbosity >= 2) cat('\nEVT tau parameters output for smoothing = ', steep, '\n', o.out$par,
                              '\n', unlist(par.trafo.tau(o.out$par)))
      z.out = local.tau.filter(tmp_data, o.out$par, kappa = kappa, givemore = TRUE)
      if (steep == smoothings[1]) tau.path = z.out$path else tau.path = cbind(tau.path, z.out$path)

      par0 = o.out$par
    }
    tau.path = tail_multiplier * tau.path[ , length(smoothings)]
  } else tau.path = externaltau
  my_data$EVT_tau = tau.path
  
  ## now ensure we are looking at the upper EVT tail (opposite of PZC!)
  tail_multiplier = ifelse(lower_tail, -1, 1)
  tmp_data = tail_multiplier * my_data$y
  tau.path = tail_multiplier * tau.path
  # ## ensure positive lower bound on tau.t
  tau_eps = 1e-2
  if (min(tau.path) < tau_eps) {
    cat('\n\nTrimming ', length(which(tau.path < tau_eps)), 'observations due to too low tau\n')
    tau.path[which(tau.path < tau_eps)] = tau_eps
  }
  
  ## now estimate the time varying tail index model
  cat("\n... NOW ESTIMATING THE TAIL INDEX")
  POTdates_idx = cumsum(as.integer(tmp_data > tau.path)); POTdates_idx[ which(POTdates_idx < 1) ] = 1
  POTidx = which(tmp_data > tau.path)
  POTdata = (tmp_data/tau.path - 1)[POTidx]
  POTdates = tmp_data[POTidx]
  par0 = c(-4) ## transformed alpha and omega
  if (!EVT_OPTIONS$FIXED_OMEGA) par0 = c(par0, -4)
  f1.quick = mean(log(1+POTdata[1:50]))
  
  my_report = list(maxit = 1000, reltol = 1e-10)
  if (verbosity >= 2) {my_report$trace = 1; my_report$REPORT = 1} else {
    if (verbosity >= 1) {my_report$trace = 1; my_report$REPORT = 10}}
  tail.out = optim(
    par0, 
    function(x) {
      par_in = par.trafo.tailindex(x, EVT_OPTIONS)
      f1 = f1.quick
      tmp = tail.filter(POTdata, par_in['omega'], par_in['beta'], par_in['alpha'], f1)$llik
      return(tmp / length(POTdata))
    },
    method = "BFGS", control = my_report
  )
  if (verbosity >= 1) cat('\nFinal parameters:\n', par.trafo.tailindex(tail.out$par, EVT_OPTIONS))
  
  ## estimate the tail index path and simulated confidence band
  par.t = par.trafo.tailindex(tail.out$par, EVT_OPTIONS)
  omega = par.t['omega']; alpha = par.t['alpha']
  tail.index = tail.filter(POTdata, omega, 1, alpha, f1.quick)
  covmatrices = compute_std_err(tail.out$par, POTdata, f1.quick, EVT_OPTIONS)
  simulations = simulate_filter_bands(tail.out$par, POTdata, covmatrices, f1.quick, nrsims = nrsims, pct = pct, raw = raw.sim, EVT_OPTIONS = EVT_OPTIONS)

  my_data$EVT_ft = tail.index$f_out[POTdates_idx]
  my_data$EVT_bandL = simulations$bands[POTdates_idx, 1]
  my_data$EVT_bandU = simulations$bands[POTdates_idx, 2]
  ## compute our VaR and EL and the true VaR
  for (alpha in EVT_OPTIONS$ALPHAS_EXTREME) {
    my_data$EVT_VaR_alpha = my_data$EVT_tau * ( (alpha / EVT_OPTIONS$TAU_TAIL_PCT) ^ -my_data$EVT_ft)
    my_data$EVT_EL_alpha = my_data$EVT_VaR_alpha / (1 - my_data$EVT_ft)
    my_data$EVT_EL_alpha[ which( my_data$EVT_ft >= 1) ] = NA  ## remove non-existent EL
    names(my_data)[(-1:0) + length(my_data)] = paste0(names(my_data)[(-1:0) + length(my_data)], alpha)
  }
  
  if (!exists('o.out')) o.out = NULL
  return(list(
    my_data = my_data,
    tau.optimizer.output = o.out,
    tail.optimizer.output = tail.out,
    covmatrices = covmatrices
  ))
}




EVT_PZC_plot <- function(my_data, alpha_extreme, sub_idx = NULL) {
  ## construct new dataframe for plotting
  dim_n = nrow(my_data)
  my_data = subset(my_data, y <= 0)
  vmy_data = rbind( 
    data.frame(idx = my_data$dates, y = my_data$y, VaRy = my_data[ , paste0("VaRtrue_alpha" , alpha_extreme)], ELy = my_data[ , paste0("ELtrue_alpha" , alpha_extreme)], VaR = 'true', ES = 'true'),
    data.frame(idx = my_data$dates, y = NA,        VaRy = my_data[ , paste0("EVT_VaR_alpha" , alpha_extreme)], ELy = my_data[ , paste0("EVT_EL_alpha" , alpha_extreme)], VaR = 'EVT',  ES = 'EVT'),
    data.frame(idx = my_data$dates, y = NA,        VaRy = my_data[ , paste0("PZC_VaR0Neld_alpha", alpha_extreme)], ELy = my_data[ , paste0("PZC_EL0Neld_alpha", alpha_extreme)], VaR = 'PZC',  ES = 'PZC')
  )
  vmy_data$VaR = factor(vmy_data$VaR, levels = c('PZC', 'EVT', 'true'), ordered = TRUE)
  vmy_data$ES = factor(vmy_data$ES, levels = c('PZC', 'EVT', 'true'), ordered = TRUE)
  if (is.null(sub_idx)) sub_idx = round(dim_n * 0.5):round(dim_n * 0.6)
  
  ## plot results
  g1 = ggplot(data = vmy_data, aes(x = idx, colour = VaR)) +
    geom_point(aes(y = y), color = "gray") +
    geom_line(aes(y = VaRy)) +
    scale_color_manual(values=c('blue', 'red', 'black')) +
    xlab('') + ylab('')
  # plot(g1)
  g1s = ggplot(data = vmy_data[ vmy_data$idx %in% sub_idx, ], aes(x = idx, colour = VaR)) +
    geom_point(aes(y = y), color = "gray") +
    geom_line(aes(y = VaRy)) +
    scale_color_manual(values=c('blue', 'red', 'black')) +
    xlab('') + ylab('')
  # plot(g1s)

  g2 = ggplot(data = vmy_data, aes(x = idx, colour = ES)) +
    geom_point(aes(y = y), color = "gray") +
    geom_line(aes(y = ELy)) +
    scale_color_discrete(name = '') +
    scale_color_manual(values=c('blue', 'red', 'black')) +
    xlab('') + ylab('')
  # plot(g2)
  g2s = ggplot(data = vmy_data[ vmy_data$idx %in% sub_idx, ], aes(x = idx, colour = ES)) +
    geom_point(aes(y = y), color = "gray") +
    geom_line(aes(y = ELy)) +
    scale_color_discrete(name = '') +
    scale_color_manual(values=c('blue', 'red', 'black')) +
    xlab('') + ylab('')
  # plot(g2s)
  
  gg = ggarrange(g1, g1s, g2, g2s, nrow = 2, ncol = 2, common.legend = TRUE)
  plot(gg)
}

