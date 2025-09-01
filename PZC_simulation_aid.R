library(ggplot2)
library(ggpubr)



get_index_from_date_vector = function(indate, from_dates, start = TRUE) {
  ## start ==
  ## TRUE: the first date index in "from_dates" on or after "indate"
  ## FALSE: the last date index in "from_dates" on or before "indate"
  indate = as.POSIXct(indate, tz = tz(from_dates[1]))
  if (start) return(min(which(from_dates >= indate)))
  return(max(which(from_dates <= indate)))
}



LFZ0_obj <- function(ydata, v_values, e_values, alpha , smoothing = 0) {
  ## smooth the step function
  if (smoothing <= 0) stepper = as.integer(ydata <= v_values) else
    stepper = 1 / (1 + exp( (ydata - v_values) / smoothing))
  
  ## compute the objective
  obj = - stepper / (alpha * e_values) * (v_values - ydata) +
    v_values / e_values + log(-e_values) - 1
  return(obj)
}



## ------------------------------------------------------
## ------------------------------------------------------




PZC_filter <- function(ydata, omega, A_mat, B_mat, f0, alpha, smoothing = 0, 
                       short_return = FALSE) {
  
  ## initialize
  if (is.null(A_mat) | is.null(B_mat) | is.null(omega)) return(1e20)
  dim_n = length(ydata)
  f_values = matrix(0, nrow = dim_n, ncol = 2)
  f_current = f0
  
  ## run the filter
  for (i1 in 1:dim_n) {
    ## store current value of time varying parameters
    f_values[i1, ] = f_current
    
    ## smooth the step function
    if (smoothing <= 0) stepper = as.integer(ydata[i1] <= f_current[1]) else
      stepper = 1 / (1 + exp( (ydata[i1] - f_current[1]) / smoothing))
    
    ## update using (11)+(12) followed by (9)+(16)
    score = c(
      -f_current[1] * (stepper - alpha),
      stepper * ydata[i1] / alpha - f_current[2]
    )
    f_current = omega + B_mat %*% f_current + A_mat %*% score
    # if (any(is.nan(f_current))) return(1e20)
    # if (max(f_current) >= 0) return(1e20)
    
  }
  
  ## return results
  obj.t = LFZ0_obj(ydata, f_values[ , 1], f_values[ , 2], alpha , smoothing)
  obj = mean(obj.t)
  if (is.nan(obj) | is.infinite(obj)) return(1e20)
  if (short_return) return(obj) else
    return(list(
      obj = obj,
      obj.t = obj.t,
      ydata = ydata,
      v = f_values[ , 1], 
      e = f_values[ , 2]
    ))
}



## ------------------------------------------------------
## ------------------------------------------------------




get_parameters <- function(theta, 
                           diagonal_A = TRUE, 
                           diagonal_B = TRUE, 
                           rw_specification = FALSE,
                           reverse = FALSE) {
  trafo_A <- function(x) {ee = exp(x); return(0.1 * (ee - 1)/(ee + 1))}
  itrafo_A <- function(x) {x = x/0.1; return(log((1 + x)/(1 - x)))}
  trafo_B <- function(x) {ee = exp(x); return((ee - 1)/(ee + 1))}
  itrafo_B <- function(x) return(log((1 + x)/(1 - x)))
  trafo_O <- function(x) {ee = exp(x); return(10 * (ee - 1)/(ee + 1))}
  itrafo_O <- function(x) {x = x/10; return(log((1 + x)/(1 - x)))}
  if (!reverse) {
    idx = 0
    omega = trafo_O(theta[ idx + (1:2)]); idx = idx + 2
    if (diagonal_A) {A_mat = diag(trafo_A(theta[ idx + (1:2)])); idx = idx + 2} else
    {A_mat = matrix(trafo_A(theta[ idx + (1:4)]), ncol = 2); idx = idx + 4}
    if (rw_specification) B_mat = diag(2) else {
      if (diagonal_B) {B_mat = diag(trafo_B(theta[ idx + (1:2)])); idx = idx + 2} else {
        B_mat = matrix(trafo_B(theta[ idx + (1:4)]), ncol = 2); idx = idx + 4}
    }
    return(list(omega = omega, A_mat = A_mat, B_mat = B_mat))
  } else {
    theta_out = itrafo_O(c(theta$omega))
    if (diagonal_A) theta_out = c(theta_out, diag(itrafo_A(theta$A_mat))) else
      theta_out = c(theta_out, itrafo_A(theta$A_mat))
    if (!rw_specification) {
      if (diagonal_B) theta_out = c(theta_out, diag(itrafo_B(theta$B_mat))) else
        theta_out = c(theta_out, itrafo_B(theta$B_mat))
    }
    return(theta_out)
  }
}



## ------------------------------------------------------
## ------------------------------------------------------




## PZC_initialize_parameters(): initialize the static parameters
##     for the Patton/Ziegel/Chen estimation
PZC_initialize_parameters <- function(ydata, alpha, first_pct = 0.1,
                                      diagonal_A, diagonal_B, rw_specification) {
  
  dim_n = length(ydata)
  y_sub = ydata[1:median(dim_n, 250, 0.1 * dim_n)]; 
  f0 = quantile(y_sub, alpha); f0 = c(f0, mean( y_sub[ y_sub <= f0] ))
  B0 = 0.995 * diag(2)
  A0 = 1e-3 * matrix(1, ncol = 2, nrow = 2)
  omega0 = quantile(ydata, alpha); omega0 = c(omega0, mean( ydata[ ydata <= omega0] )); omega0 = (diag(2) - B0) %*% f0
  theta0 = get_parameters(list(omega = omega0, A_mat = A0, B_mat = B0),
                          diagonal_A = diagonal_A,
                          diagonal_B = diagonal_B,
                          rw_specification = rw_specification,
                          reverse = TRUE)
  return(list(theta0 = theta0, f0 = f0))
}




## ------------------------------------------------------
## ------------------------------------------------------




PZC_optimize <- function(
    my_dataframe,
    filter = PZC_filter,
    smoothings = c(1/10, 1/25, 0, -1),
    init_par = NULL,
    PZC_OPTIONS = list(),
    VaR_output_name = "PZC_VaR",
    ES_output_name = "PZC_ES",
    verbosity = 0,
    in_sample_idx = NULL,
    out_of_sample_idx = NULL
) {
  
  ## IMPORTANT:
  ## the code follows the paper, inclusive of its appendix on the repeated 
  ## smoothed estimation and final non Newton optimization, all for the LEFT
  ## tail. As the code is integrated in an approach that concentrates on the
  ## RIGHT tail, returns are flipped at the beginning, negative VaR and ES
  ## filtered, and then the data, VaR, and ES flipped back to positive 
  ## territory
  my_dataframe$y = -my_dataframe$y # flip data
  
  ## check in and out of sample indices
  if (is.null(in_sample_idx)) in_sample_idx = 1:nrow(my_dataframe)
  if (is.null(out_of_sample_idx)) out_of_sample_idx = in_sample_idx
  if (min(in_sample_idx) < 1) error("PZC_optimize(): in-sample starts before obs #1")
  if (min(out_of_sample_idx) < 1) error("PZC_optimize(): out-of-sample starts before obs #1")
  if (max(in_sample_idx) > nrow(my_dataframe)) error("PZC_optimize(): in-sample ends after data ends")
  if (max(out_of_sample_idx) > nrow(my_dataframe)) error("PZC_optimize(): out-of-sample ends after data ends")
  
  ## check parameters for optimization and set default values
  diagonal_A = ifelse(is.null(PZC_OPTIONS$diagonal_A), FALSE, PZC_OPTIONS$diagonal_A)
  diagonal_B = ifelse(is.null(PZC_OPTIONS$diagonal_B), TRUE, PZC_OPTIONS$diagonal_B)
  rw_specification = ifelse(is.null(PZC_OPTIONS$rw_specification), FALSE, PZC_OPTIONS$rw_specification)
  if (is.null(PZC_OPTIONS$ALPHA_EXTREME)) alpha = 0.05 else alpha = PZC_OPTIONS$ALPHA_EXTREME
  all.opt.out = list(); all.opt.count = 1
  
  ###############################################
  ## initialize PZC parameters for estimation
  ###############################################
  if (is.null(init_par)) theta0 = 
    PZC_initialize_parameters(
      my_dataframe$y[in_sample_idx], alpha, 
      diagonal_A = diagonal_A, 
      diagonal_B = diagonal_B, 
      rw_specification = PZC_OPTIONS$rw_specification) else theta0 = init_par
  f0 = theta0$f0
  theta0 = theta0$theta0
  
  
  ## estimate for smoothed indicator function (cascade)
  for (smoothing in smoothings) {
    
    if (verbosity > 0) {
      cat("\n#######################")
      cat("\n## alpha     = ", alpha)
      cat("\n## smoothing = ", smoothing)
      cat("\n#######################\n\n")
    }
    
    if (smoothing < 0) {
      my_report = list(maxit = 2000)
      # my_report = list(maxit = 2000, REPORT = 100, trace = 1)
      if (verbosity >= 2) {my_report$trace = 1; my_report$REPORT = 1}
      a.out = optim(theta0,
                    function(x) {
                      x = get_parameters(x, diagonal_A = diagonal_A, diagonal_B = diagonal_B, rw_specification = rw_specification)
                      return(filter(my_dataframe$y[in_sample_idx], x$omega, x$A_mat, x$B_mat, f0, alpha, smoothing)$obj)
                    },
                    method = "Nelder-Mead",
                    control = my_report)
    } else { 
      # my_report = list(maxit = 50)
      my_report = list(maxit = 2000)
      if (verbosity >= 2) {my_report$trace = 1; my_report$REPORT = 1} else {
        if (verbosity >= 1) {my_report$trace = 1; my_report$REPORT = 10}
      }
      a.out = optim(theta0,
                    function(x) {
                      x = get_parameters(x, diagonal_A = diagonal_A, diagonal_B = diagonal_B, rw_specification = rw_specification)
                      return(filter(my_dataframe$y[in_sample_idx], x$omega, x$A_mat, x$B_mat, f0, alpha, smoothing)$obj)
                    },
                    method = "BFGS",
                    control = my_report)
    }
    a.out$alpha = alpha; a.out$smoothing = smoothing; all.opt.out[[all.opt.count]] = a.out; all.opt.count = all.opt.count + 1;
    theta0 = a.out$par
    p.out = get_parameters(theta0, diagonal_A, diagonal_B, rw_specification = rw_specification)
    if (verbosity >= 1) print(p.out)
    
    ## store FLIPPED to right versions of the just computed left-tail VaR and ES 
    ## (only for the last smoothing)
    b.out = filter(my_dataframe$y[out_of_sample_idx], p.out$omega, p.out$A_mat, p.out$B_mat, f0, alpha, smoothing)
    if (is.numeric(b.out)) browser()
    my_dataframe[out_of_sample_idx, VaR_output_name] = -b.out$VaR
    my_dataframe[out_of_sample_idx, ES_output_name] = -b.out$EL
  }
  
  my_dataframe$y = -my_dataframe$y # flip data back
  return(list(
    my_data = my_dataframe,
    optimizer_output = all.opt.out
  ))
}




## ------------------------------------------------------
## ------------------------------------------------------




plot_VaR <- function(my_dataframe, alphas, sub_frame_idx = NULL,
                     gg_plt_extras = NULL, make_plot = TRUE) {
  ## plotting
  if (is.null(sub_frame_idx)) sub_frame_idx = 1:nrow(my_dataframe)
  gg = lapply(
    unique(alphas),
    function(alpha) {
      out1 = ggplot(data = my_dataframe[sub_frame_idx, ], aes(x = dates)) +
        geom_point(aes(y = y), color = 'gray', fill = 'gray') +
        geom_line(aes(y = .data[[paste0('PZC_VaR0Neld_alpha', alpha)]]), color = 'red') +
        geom_line(aes(y = .data[[paste0('VaRtrue_alpha', alpha)]]), color = 'black') +
        geom_line(aes(y = .data[[paste0('ELtrue_alpha', alpha)]]), color = 'green') +
        geom_line(aes(y = .data[[paste0('PZC_EL0Neld_alpha', alpha)]]), color = 'forestgreen')
      if (!is.null(gg_plt_extras)) out1 = out1 + gg_plt_extras
      return(out1)
    })
  gg = ggarrange(plotlist = gg, nrow = 1)
  if (make_plot) plot(gg)
  return(gg)
}


## ------------------------------------------------------
## ------------------------------------------------------




save_PZC_EVT_frame <- function(
    my_dataframe, store_filename = NULL, 
    sub_frame_idx = NULL, select_names = NULL, new_names = NULL) {
  ## saving file
  if (!is.null(store_filename)) {
    if (is.null(sub_frame_idx)) sub_frame_idx = 1:nrow(my_dataframe)
    if (is.null(select_names)) select_names = names(my_dataframe)
    if (is.null(new_names)) new_names = select_names
    
    # select the subdataframe and rename variables
    my_dataframe = my_dataframe[sub_frame_idx, select_names]
    for (i1 in 1:length(select_names)) names(my_dataframe)[i1] = new_names[i1]
    
    # store output and make graph
    # do not remove braces
    while (file.exists(store_filename)) {
      answer1<-readline(prompt=paste("File", store_filename, "exists: overwrite (y/n)? "))
      if ((answer1 != "y") & (answer1 != "Y")) {
        answer1<-readline(prompt=paste("Provide new filename (without .csv): "))
        store_filename = paste0(answer1, ".csv")
      } else break
    } 
    write.csv(my_dataframe, file = store_filename, row.names = FALSE)
  } else {
    cat("\n\n*** NO FILENAME PROVIDED IN save_PZC_EVT_frame(), SO NOT STORING ANYHTING ***")
  }
}






NZ_table = function(my_data, data_name, evaluation_idx = NULL,
                    VaR_names, ES_names, alphas, method_names = NULL,
                    risk_measures_are_for_right_tail = TRUE,
                    VaR_true_name = NULL, 
                    ES_true_name = NULL,
                    scale = 1) {
  # alphas are close to one, and all the risk measures should be for the right tail of the data!!!
  
  if (is.null(evaluation_idx)) evaluation_idx = 1:nrow(my_data)
  my_data = my_data[ evaluation_idx, c(data_name, VaR_names, ES_names, VaR_true_name, ES_true_name)]
  my_data = my_data[ complete.cases(my_data), ]
  out_table = data.frame()
  
  ## check arguments on consistency of formats
  if (length(VaR_names) != length(ES_names)) error("NZ_table(): VaR_names and ES_names do not have the same length")
  if (length(alphas) == 1) alphas = rep(alphas, length(VaR_names))
  if (length(VaR_names) != length(alphas)) error("NZ_table(): VaR_names and alphas do not have the same length")
  tail_multiplies = ifelse(risk_measures_are_for_right_tail, 1, -1)
  if (is.null(method_names)) method_names = VaR_names
  if (length(VaR_names) != length(method_names)) error("NZ_table(): VaR_names and method_names do not have the same length")
  
  
  evaluation_idx = 1:nrow(my_data)
  for (i1 in 1:length(VaR_names)) {
    x = my_data[ , data_name]
    r1 = my_data[ , VaR_names[i1]]
    r2 = my_data[ , ES_names[i1]]

    tmp1 = data.frame()
    tmp1[evaluation_idx, "NZ_2.19"] = scale * ( (1 - alphas[i1]) * r1 + ifelse(x > r1, x - r1, 0) )
    tmp1[evaluation_idx, "NZ_2.20"] = scale * ( (1 - alphas[i1]) * log(r1) + ifelse(x > r1, log(pmax(1e-8, x/r1)), 0) )
    tmp1[evaluation_idx, "NZ_2.23"] = scale * ( (x > r1) * 0.5 * (x - r1) / sqrt(r2) + (1 - alphas[i1]) * 0.5 * (r1 + r2) / sqrt(r2) )
    tmp1[evaluation_idx, "NZ_2.24"] = scale * ( (x > r1) * (x - r1) / r2 + (1 - alphas[i1]) * (r1 / r2 - 1 + log(r2)) )
    if (i1 == 1) benchmark_method_performance = tmp1
    
    for (nms in c("NZ_2.19", "NZ_2.20", "NZ_2.23", "NZ_2.24")) {
      out_table[ method_names[i1], nms] = mean(tmp1[ , nms])
      if (i1 > 1) {
        aid1 = dm.test(benchmark_method_performance[ , nms], tmp1[ , nms], alternative = "two.sided", h = 1, power = 1)
        out_table[ paste0(method_names[i1], "_DM"), nms] = aid1$statistic
        out_table[ paste0(method_names[i1], "_DMp"), nms] = aid1$p.value
      }
    }

    if (!is.null(VaR_true_name)) {
      r1true = my_data[evaluation_idx, VaR_true_name]
      out_table[ method_names[i1], "VaR_MSE"] = mean( (r1 - r1true)^2 )
      out_table[ method_names[i1], "VaR_MAE"] = mean( abs(r1 - r1true) )
    }
    if (!is.null(ES_true_name)) {
      r2true = my_data[evaluation_idx, ES_true_name]
      out_table[ method_names[i1], "ES_MSE"] = mean( (r2 - r2true)^2 )
      out_table[ method_names[i1], "ES_MAE"] = mean( abs(r2 - r2true) )
    }
    
    if (i1 < length(VaR_names)) out_table = rbind(out_table, NA)
    
  }
  
  return(out_table)
}






read_empirical_data = function(filename, price_column_name, date_column_name,
                               prices_are_already_returns = FALSE, reverse = FALSE,
                               date_format = NULL, date_POSX_format = NULL, tz = "",
                               flip_return = FALSE, begin_date = NULL, end_date = NULL) {
  ##
  ## IMPORTANT: THE SERIES ARE RETURNED SUCH THAT THE ANALYSIS IS DONE 
  ##            ON THE RIGHT(!!) TAIL
  ##            TAU, VaR, ES are reported also for the RIGHT tail
  ##
  ## For USDRUB, BTCUSD, ETHUSD we look into the RIGHT tail of the returns
  ## For EURUSD, we want to look into the left tail, so for that series 
  ##             we need to set flip_return = TRUE, such that we look into 
  ##             the RIGHT tail of the negative returns
  ##

  ## read the data
  extension = stringi::stri_sub(filename, from = -4)
  if (extension == ".csv") {
    my_data = as.data.frame(read.csv(filename)[ , c(date_column_name, price_column_name)])
  } else if (extension == "xlsx") {
    my_data = as.data.frame(readxl::read_xlsx(filename)[ , c(date_column_name, price_column_name)])
  } else error(paste0("ERROR in read_empirical_data(): wrong file type of ", filename))
  ## label columns with fixed names and convert dates
  names(my_data) = c("dates", "y")
  if (reverse) my_data = my_data[nrow(my_data):1, ]
  if (!is.null(date_format)) my_data$dates = as.Date(my_data$dates, date_format)
  if (!is.null(date_POSX_format)) my_data$dates = as.POSIXct(my_data$dates, date_POSX_format, tz = tz)
  ## sample subselection
  if (!is.null(begin_date)) my_data = subset(my_data, dates >= begin_date)
  if (!is.null(end_date)) my_data = subset(my_data, dates <= end_date)
  ## construct 100 * log returns, flip them if needed, and remove missings observations
  if (!prices_are_already_returns) my_data$y = 
    c(NA, ifelse(flip_return, -100, 100) * diff(log(my_data$y)))
  return(my_data[complete.cases(my_data), ])
}






GARCH_optimizer = function(
    my_data,
    verbosity = 0,
    in_sample_idx = NULL,
    out_of_sample_idx = NULL,
    GARCH_optimizer_outputs = NULL
) {
  
  # Specify GARCH(1,1) with Normal innovations
  if (verbosity > 0) cat('\nEstimating GARCH-N')
  if (is.null(GARCH_optimizer_outputs)) {
    GARCH_optimizer_outputs_new = list()
    spec_garch <- ugarchspec(
      variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
      mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
      distribution.model = "norm"
    )
    fit_garch <- ugarchfit(spec = spec_garch, data = my_data$y[in_sample_idx], 
                           solver = 'gosolnp', solver.control = list(n.restarts = 10, parallel = TRUE, pkg = 'multicore', cores = 14))
    fit_garch = as.list(coef(fit_garch))
    GARCH_optimizer_outputs_new$fit_GARCH = fit_garch
  } else {
    fit_garch = GARCH_optimizer_outputs$fit_GARCH
    GARCH_optimizer_outputs_new = GARCH_optimizer_outputs
  }
  spec_garch <- ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
    mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
    distribution.model = "norm",
    fixed.pars = fit_garch
  )
  fit_garch <- ugarchfilter(spec=spec_garch, data=my_data$y[out_of_sample_idx])
  my_data[out_of_sample_idx, "GARCH_VaR"] = coef(fit_garch)["mu"] - qnorm(alpha_extreme) * sigma(fit_garch)
  my_data[out_of_sample_idx, "GARCH_ES"] = coef(fit_garch)["mu"] + (dnorm(qnorm(alpha_extreme)) / alpha_extreme) * sigma(fit_garch)
  
  # Specify GARCH(1,1) with Student t innovations
  if (verbosity > 0) cat('\nEstimating GARCH-t')
  if (is.null(GARCH_optimizer_outputs)) {
    spec_garch <- ugarchspec(
      variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
      mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
      distribution.model = "std"
    )
    fit_garch <- ugarchfit(spec = spec_garch, data = my_data$y[in_sample_idx], 
                           solver = 'gosolnp', solver.control = list(n.restarts = 10, parallel = TRUE, pkg = 'snowfall', cores = 14))
    fit_garch = as.list(coef(fit_garch))
    GARCH_optimizer_outputs_new$fit_GARCHt = fit_garch
  } else fit_garch = GARCH_optimizer_outputs$fit_GARCHt
  spec_garch <- ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
    mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
    distribution.model = "std",
    fixed.pars = fit_garch
  )
  fit_garch <- ugarchfilter(spec=spec_garch, data=my_data$y[out_of_sample_idx])
  nu = coef(fit_garch)["shape"]
  my_data[out_of_sample_idx, "GARCHt_VaR"] = coef(fit_garch)["mu"] - qt(alpha_extreme, df = nu) * sigma(fit_garch)
  my_data[out_of_sample_idx, "GARCHt_ES"] = coef(fit_garch)["mu"] + (
    dt(qt(alpha_extreme, df = nu), df = nu) / alpha_extreme * (nu + (qt(alpha_extreme, df = nu))^2) / (nu - 1)
  ) * sigma(fit_garch)
  
  return(list(
    my_data = my_data,
    GARCH_optimizer_outputs = GARCH_optimizer_outputs_new
  ))
}






print_NZ_table = function(in_sample_data, out_of_sample_data,
                          alpha_extreme, asset_name, with_print = TRUE,
                          NZ_table_in = NULL, NZ_table_sep = "  ---  ") {
  out_print = c(
    "=========================================",
    paste0("IN-SAMPLE RESULTS (", alpha_extreme, ",", asset_name, ")"),
    "=========================================",
    capture.output(round(NZ_table(
      my_data, data_name = "y", evaluation_idx = out_of_sample_idx,
      VaR_names = c('EVT_VaR', 'PZC_VaR', 'GARCH_VaR', 'GARCHt_VaR'), 
      ES_names = c('EVT_ES', 'PZC_ES', 'GARCH_ES', 'GARCHt_ES'),
      alphas = 1 - alpha_extreme, 
      method_names = c("EVT", "PZC", "GARCH", "GARCHt")),3)),
    "=========================================",
    paste0("OUT-OF-SAMPLE RESULTS (", alpha_extreme, ",", asset_name, ")"),
    "=========================================",
    paste0(asset_name, " NZ Results (tau alpha = ", alpha_tail, ")"),
    capture.output(round(NZ_table(
      my_data_oos, data_name = "y", 
      VaR_names = c('EVT_VaR', 'PZC_VaR', 'GARCH_VaR', 'GARCHt_VaR'), 
      ES_names = c('EVT_ES', 'PZC_ES', 'GARCH_ES', 'GARCHt_ES'),
      alphas = 1 - alpha_extreme, 
      method_names = c("EVT", "PZC", "GARCH", "GARCHt")),3))
  )
  out_print = paste(NZ_table_in, out_print, sep = NZ_table_sep)
  if(with_print) writeLines(out_print)
  return(out_print)
}



