library(ggplot2)
library(ggpubr)




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
    verbosity = 0
) {
  
  diagonal_A = ifelse(is.null(PZC_OPTIONS$diagonal_A), FALSE, PZC_OPTIONS$diagonal_A)
  diagonal_B = ifelse(is.null(PZC_OPTIONS$diagonal_B), TRUE, PZC_OPTIONS$diagonal_B)
  rw_specification = ifelse(is.null(PZC_OPTIONS$rw_specification), FALSE, PZC_OPTIONS$rw_specification)
  if (is.null(PZC_OPTIONS$ALPHAS_EXTREME)) alphas = 0.05 else alphas = PZC_OPTIONS$ALPHAS_EXTREME
  all.opt.out = list(); all.opt.count = 1
  
  for (alpha in unique(alphas)) {
    ###############################################
    ## initialize PZC parameters for estimation
    ###############################################
    if (is.null(init_par)) theta0 = 
        PZC_initialize_parameters(
          my_dataframe$y, alpha, 
          diagonal_A = diagonal_A, 
          diagonal_B = diagonal_B, 
          rw_specification = PZC_OPTIONS$rw_specification) else theta0 = init_par
    f0 = theta0$f0
    theta0 = theta0$theta0
    

    ## estimate for smoothed indicator function (cascade)
    for (smoothing in smoothings) {
      
      cat("\n#######################")
      cat("\n## alpha     = ", alpha)
      cat("\n## smoothing = ", smoothing)
      cat("\n#######################\n\n")
      
      if (smoothing < 0) {
        my_report = list(maxit = 2000)
        if (verbosity >= 2) {my_report$trace = 1; my_report$REPORT = 1}
        a.out = optim(theta0,
                      function(x) {
                        x = get_parameters(x, diagonal_A = diagonal_A, diagonal_B = diagonal_B, rw_specification = rw_specification)
                        return(filter(my_dataframe$y, x$omega, x$A_mat, x$B_mat, f0, alpha, smoothing)$obj)
                      },
                      method = "Nelder-Mead",
                      control = my_report)
      } else { 
        my_report = list(maxit = 50)
        if (verbosity >= 2) {my_report$trace = 1; my_report$REPORT = 1} else {
          if (verbosity >= 1) {my_report$trace = 1; my_report$REPORT = 10}
        }
        a.out = optim(theta0,
                      function(x) {
                        x = get_parameters(x, diagonal_A = diagonal_A, diagonal_B = diagonal_B, rw_specification = rw_specification)
                        return(filter(my_dataframe$y, x$omega, x$A_mat, x$B_mat, f0, alpha, smoothing)$obj)
                      },
                      method = "BFGS",
                      control = my_report)
      }
      a.out$alpha = alpha; a.out$smoothing = smoothing; all.opt.out[[all.opt.count]] = a.out; all.opt.count = all.opt.count + 1;
      theta0 = a.out$par
      p.out = get_parameters(theta0, diagonal_A, diagonal_B, rw_specification = rw_specification)
      if (verbosity >= 1) print(p.out)
      
      ## store VaR and EL for this smoothing
      b.out = filter(my_dataframe$y, p.out$omega, p.out$A_mat, p.out$B_mat, f0, alpha, smoothing)
      if (is.numeric(b.out)) browser()
      my_dataframe$VaR = b.out$VaR
      my_dataframe$EL = b.out$EL
      names(my_dataframe)[length(my_dataframe) + (-1:0)] = paste0(
        'PZC_', c('VaR', 'EL'), max(smoothing, 0), ifelse(smoothing < 0, 'Neld', ''),
        '_alpha', alpha)   
    }
  }
  
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


