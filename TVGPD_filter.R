library("pracma")
library("progress")



## TVGPD.filter() filter the dynamic tail parameter using the integrated filter
## [NB: the UPPER tail shape dynamics are estimated]
##
## INPUTS:
## =======
## POTdata: vector with Peaks Over Thresholds (POT) data
## omega: the 2x1 intercept of the transition equation
## beta: the 2x1 loading of the previous tail parameter in the transition equation
## alpha: the 2x1 loading of the scaled score in the transition equation
## f1: the initial tail parameter value; order (logxi, logdelta)
##
## OUTPUTS:
## =======
## list():
##    llik: the log-likelihood value of the GPD approximation
##    llik.t: the vector of log-likelihood contributions
##    f_out: vector of filtered tail parameters
##    TVGPD_out: vector with columns xi and delta
TVGPD.filter <- function(POTdata, omega, beta, alpha, f1) {
  
  ## initialize
  dimn = length(POTdata)
  f_out = matrix(0, nrow = dimn, ncol = 2)
  colnames(f_out) = c('logxi', 'logdelta')
  llik = 0
  
  ## run the filter over all observations
  for (i1 in 1:dimn) {
    ## store value
    f_out[i1, ] = f1
    ## update the parameter using Eq.(7) from D'Innocenzo et al.
    xi = exp(f1[1]); delta = exp(f1[2]); xt = POTdata[i1]
    st = c(
      (1 + xi) / xi^2 * log(1 + xi * xt / delta) + 
        (delta - (xi + 3 + 1/xi) * xt) / (delta + xi * xt),
      sqrt(1 + 2 * xi) * (xt - delta) / (delta + xi * xt)
    )
    f1 = omega + beta * f1 + alpha * st
  }
  ## compute the likelihood
  TVGPD_out = exp(f_out)
  colnames(TVGPD_out) = c('xi', 'delta')
  llik.t = -f_out[ , 'logdelta'] - (1 + 1/TVGPD_out[ , 'xi']) * 
    log(1 + TVGPD_out[ , 'xi'] * POTdata / TVGPD_out[ , 'delta'])
  llik = sum(llik.t)
  
  ## check log-lik for integrity
  if (is.na(llik) | is.infinite(llik)) llik = -1e10
  
  return(list(
    llik = llik,
    llik.t = llik.t,
    f_out = f_out,
    TVGPD_out = TVGPD_out
  ))
}




TVGPD.par.trafo = function(x) {
  omega = x[1:2]
  alpha = 1 / (1 + exp(-x[3:4]))
  beta = 1 / (1 + exp(-x[5:6]))
  f1 = omega / (1 - beta)
  return(list(
    omega = omega,
    alpha = alpha,
    beta = beta,
    f1 = f1
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
## TVGPD_OPTIONS$TVGPD_tau: vector of tau (POT thresholds); must be there!
## verbosity: integer 0..5, higher for more intermediate output;
## in_sample_idx: a vector of rows of my_data$y that is used for estimation
## out_of_sample_idx: a vector of rows of my_data$y that is used for evaluation of the filter, VaR, and ES. The out_of_sample_idx may also fully contain the in_sample_idx (plus extra observations)
##
## OUTPUT:
## =======
## list():
##    $my_data: original dataframe plus additional columns
##    $tau.optimizer.output: output from the optim() for the thresholds
##    $tail.optimizer.output: output of optim() for the integrated score-driven tail model
##    $bands: matrix, holding the POT times f.t and the upper and lower confidence bands
##    $covmatrices: list of covariance matrices for the tail index model (.raw for the untransformed parameters, otherwise for the omega and alpha)
TVGPD_estimate_full_model <- function(
    my_data, TVGPD_OPTIONS, verbosity = 0,
    in_sample_idx = NULL, out_of_sample_idx = NULL
) {

  dimT = nrow(my_data)
  if (verbosity > 0) cat('\n#########################\n## TVGPD model optimization\n#########################\n')
  

  # ## ensure positive lower bound on tau.t
  tau_eps = 1e-2
  union_indices = union(in_sample_idx, out_of_sample_idx)
  if (min(my_data$TVGPD_tau[union_indices]) < tau_eps) {
    if (verbosity > 0) cat('\n\nTrimming ', length(which(my_data$TVGPD_tau[union_indices] < tau_eps)), 'observations due to too low tau\n')
    my_data$TVGPD_tau[union_indices] = pmax(my_data$TVGPD_tau, tau_eps)
  }
  
  ## now estimate the time varying tail index model
  if (verbosity > 0) cat("\n... NOW ESTIMATING THE TAIL INDEX")
  POTdates_idx = cumsum(as.integer(my_data$y > my_data$TVGPD_tau))
  POTdates_idx[ which(POTdates_idx < 1) ] = 1
  POTdates_idx_is = POTdates_idx[in_sample_idx]
  POTdates_idx_oos = POTdates_idx[out_of_sample_idx]
  POTidx = which(my_data$y > my_data$TVGPD_tau)
  POTidx_is = which(my_data$y[in_sample_idx] > my_data$TVGPD_tau[in_sample_idx])
  POTidx_oos = which(my_data$y[out_of_sample_idx] > my_data$TVGPD_tau[out_of_sample_idx])
  POTdata = (my_data$y - my_data$TVGPD_tau)[POTidx]
  POTdata_is = (my_data$y - my_data$TVGPD_tau)[POTidx_is]
  POTdata_oos = (my_data$y - my_data$TVGPD_tau)[POTidx_oos]
  ## transformed omega, alpha and beta
  par0 = c(0.02 * log(0.2), 0.02 * log(mean(POTdata_is)), -2, -2, 4, 4) 
  # ## set xi to same initial value as EVT
  # start_up_length = 50
  # f1.quick = mean(log(1 + POTdata_is[1:start_up_length] / TVGPD_tau[in_sample_idx[1:start_up_length]]))
  # ## then set delta to weighted scale
  # aid1 = 1; for (i1 in 1:5) aid = mean( 
  #   (1 + f1.quick) * POTdata_is[1:start_up_length] / (1 + f1.quick * POTdata_is[1:start_up_length] / aid1))
  # ## then take logs
  # f1.quick = log(c(f1.quick, aid1))
  
  my_report = list(maxit = 1000, reltol = 1e-10)
  if (verbosity >= 2) {my_report$trace = 1; my_report$REPORT = 1} else {
    if (verbosity >= 1) {my_report$trace = 1; my_report$REPORT = 10}}
  tail.out = optim(
    par0, 
    function(x) {
      par_in = TVGPD.par.trafo(x)
      tmp = -TVGPD.filter(POTdata_is, par_in$omega, par_in$beta, par_in$alpha, par_in$f1)$llik
      return(tmp / length(POTdata_is))
    },
    method = "BFGS", control = my_report
  )
  if (verbosity >= 1) cat('\nFinal parameters:\n', TVGPD.par.trafo(tail.out$par))
  
  ## estimate the tail index path and simulated confidence band
  par_in = TVGPD.par.trafo(tail.out$par)
  tail.index = TVGPD.filter(POTdata_oos, par_in$omega, par_in$beta, par_in$alpha, par_in$f1)
  my_data[out_of_sample_idx, "TVGPD_xi"] = tail.index$TVGPD_out[POTdates_idx_oos, 'xi']
  my_data[out_of_sample_idx, "TVGPD_delta"] = tail.index$TVGPD_out[POTdates_idx_oos, 'delta']
  ## compute our VaR and EL and the true VaR
  my_data[out_of_sample_idx, "TVGPD_VaR"] = my_data$TVGPD_tau[out_of_sample_idx] + my_data$TVGPD_delta[out_of_sample_idx] / my_data$TVGPD_xi[out_of_sample_idx] * ( (TVGPD_OPTIONS$TAU_TAIL_PCT / TVGPD_OPTIONS$ALPHA_EXTREME)^my_data$TVGPD_xi[out_of_sample_idx] - 1)
  my_data[out_of_sample_idx, "TVGPD_ES"] = (my_data$TVGPD_VaR[out_of_sample_idx] + my_data$TVGPD_delta[out_of_sample_idx] - my_data$TVGPD_xi[out_of_sample_idx] * my_data$TVGPD_tau[out_of_sample_idx]) / (1 - my_data$TVGPD_xi[out_of_sample_idx])
  my_data[which( my_data$TVGPD_xi[out_of_sample_idx] >= 1), "TVGPD_ES"] = NA  ## remove non-existent ES

  return(list(
    my_data = my_data,
    tail.optimizer.output = tail.out
  ))
}






TVGPD_optimize <- function(
    my_data, TVGPD_OPTIONS, 
    PZC_OPTIONS = NULL,
    verbosity = 0,
    in_sample_idx = NULL,
    out_of_sample_idx = NULL
) {
  
  ## check in and out of sample indices
  if (is.null(in_sample_idx)) in_sample_idx = 1:nrow(my_data)
  if (is.null(out_of_sample_idx)) out_of_sample_idx = in_sample_idx
  if (min(in_sample_idx) < 1) error("TVGPD_optimize(): in-sample starts before obs #1")
  if (min(out_of_sample_idx) < 1) error("TVGPD_optimize(): out-of-sample starts before obs #1")
  if (max(in_sample_idx) > nrow(my_data)) error("TVGPD_optimize(): in-sample ends after data ends")
  if (max(out_of_sample_idx) > nrow(my_data)) error("TVGPD_optimize(): out-of-sample ends after data ends")

  if (is.null(TVGPD_OPTIONS$EXTERNAL_TAU)) {
    ## estimate tau using PZC
    ## make sure that out_of_sample_idx encompasses or tightly follows in_sample_idx
    all_idx = sort(union(in_sample_idx, out_of_sample_idx))
    PZC_OPTIONS$ALPHA_EXTREME = TVGPD_OPTIONS$TAU_TAIL_PCT
    PZC_optimizer_outputs = 
      PZC_optimize(my_data, filter = PZC_filter_cpp, 
                   PZC_OPTIONS = PZC_OPTIONS, verbosity = 0,
                   in_sample_idx = in_sample_idx,
                   out_of_sample_idx = out_of_sample_idx)
    my_data$TVGPD_tau = PZC_optimizer_outputs$my_data$PZC_VaR; PZC_optimizer_outputs$my_data = NULL
  } else my_data$TVGPD_tau = TVGPD_OPTIONS$EXTERNAL_TAU

  aid1 = TVGPD_estimate_full_model(
    my_data, TVGPD_OPTIONS = TVGPD_OPTIONS, verbosity = verbosity,
    in_sample_idx = in_sample_idx, out_of_sample_idx = out_of_sample_idx)

  return(aid1)
}



