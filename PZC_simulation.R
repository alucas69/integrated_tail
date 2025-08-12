###################
## Set up functions
###################
library(ggplot2)
library(ggpubr)
library(Rcpp)
library(RcppArmadillo)

source('PZC_simulation_aid.R')
source('Integrated_EVT_filter.R')
sourceCpp('PZC_simulation.cpp')
sourceCpp('Integrated_EVT_filter.cpp')


#####################################
## Set tail percentage for thresholds
## and extreme tail percentages for
## VaR and EL
#####################################
alpha_tail = 0.05
alphas_extreme = c(0.05, 1e-3)


###############################
## set model configurations PZC
## PZC = Patton/Ziegel/Chen
##       J.Econometrics
###############################
s_PZC_OPTIONS = list(
  diagonal_A = FALSE,        # boolean; TRUE: use diagonal A matrix 
  diagonal_B = TRUE,         # boolean; TRUE: use diagonal B matrix 
  ALPHAS_EXTREME = alphas_extreme,  # double, tail percentage for the thresholds tau
  rw_specification = FALSE,  # boolean; TRUE: use unit B matrix; overrides diagonal_B
  set_seed = 1234            # NULL or integer: if NULL, do not set the seed when generating data
)

###############################
## set model configurations EVT
###############################
s_EVT_OPTIONS = list(
  FIXED_OMEGA = FALSE,       # logical, if TRUE, do not estimate omega for tail dynamics, 
                             #          but fix to s_FIXED_OMEGA_VALUE
  FIXED_OMEGA_VALUE = 0,     # double, value of omega for tail dynamics if s_FIXED_OMEGA == TRUE
  BANDS_RAW = FALSE,         # logical, if FALSE base confidence bands on (omega, alpha),
                             #          otherwise base bands on free (theta1, theta2)
                             #          parameters and transform to (omega, alpha)
  BANDS_NRSIMS = 1000,       # integer, number of simulations for the confidence bands
  BANDS_PCT = 0.75,          # double, confidence level for the confidence bands
  TAU_TAIL_PCT = alpha_tail, # double, tail percentage for the thresholds tau
  ALPHAS_EXTREME = alphas_extreme,  # double, tail percentage for the thresholds tau
  EXTERNAL_TAU = NULL,       # can be overwritten if external thresholds available
  USE_PZC = TRUE,            # logical, if TRUE use Patton/Ziegel/Chen taus
  set_seed = s_PZC_OPTIONS$set_seed # NULL or integer: if NULL, do not set the seed when generating data
)



#######################################################
## Simulation routines: 
## examples where VaR/EL are in the *LEFT* tail
#######################################################
## STUDENT T WITH SCALE 1 AND TV NU
generate_data_student_t_v1 <- function(
    dim_n, alphas, df_min = 2.5, df_max = 10, nr_cycles = 3, 
    set_seed = NULL) {

  VaRtrue = function(sigma, df_inv, alpha) {
    qtn = qt(alpha, df = 1/df_inv)
    return(sigma * qtn)
  }
  ELtrue = function(sigma, df_inv, alpha) {
    qtn = qt(alpha, df = 1/df_inv)
    return(-sigma_sim * (1 + df_sim_inv * qtn^2)/(1 - df_sim_inv) * dt(qtn, df = 1/df_sim_inv) / alpha)
  }
  if (is.numeric(set_seed)) {cat(paste("\nSimulation seed set to", set_seed)); set.seed(set_seed)}
  ## actual generation
  df_sim_inv  = 1/df_max + 0.5 * (1/df_min - 1/df_max) * (1 - cos(nr_cycles * (1:dim_n) / dim_n * 2 * pi))
  sigma_sim = -(1 + 4 * abs((1:dim_n)/dim_n - 0.5)) / qt(0.05, df = 1 / df_sim_inv) # rep(1, dim_n)
  my_data = data.frame(dates = 1:dim_n, y = rt(dim_n, df = 1/df_sim_inv) * sigma_sim,
                       df_inv = df_sim_inv)
  for (alpha in unique(alphas)) { 
    my_data$VaRtrue = VaRtrue(sigma_sim, df_sim_inv, alpha)
    my_data$ELtrue = ELtrue(sigma_sim, df_sim_inv, alpha)
    names(my_data)[(-1:0) + length(my_data)] = paste0(c('VaRtrue', 'ELtrue'), '_alpha', alpha)
  }
  return(my_data)
}



################
## Generate data
################
dim_n = 200000
my_data = generate_data_student_t_v1(dim_n, alphas_extreme, 
                                     nr_cycles = 3,
                                     set_seed = s_PZC_OPTIONS$set_seed)


###############################################
## estimate PZC specification for extreme alpha
###############################################
PZC_optimizer_outputs = 
  PZC_optimize(my_data, filter = PZC_filter_cpp, 
               PZC_OPTIONS = s_PZC_OPTIONS, verbosity = 0)
my_data = PZC_optimizer_outputs$my_data; PZC_optimizer_outputs$my_data = NULL
###############################################
## plot PZC results
###############################################
plot_VaR(my_data, alphas = s_PZC_OPTIONS$ALPHAS_EXTREME, 
         sub_frame_idx = 75e3:80e3,
         gg_plt_extras = ylim(-20,0),
         make_plot = FALSE)





###############################################
## estimate EVT specification for extreme alpha
###############################################
EVT_optimizer_outputs = EVT_optimize(my_data, EVT_OPTIONS = s_EVT_OPTIONS, verbosity = 0)
my_data = EVT_optimizer_outputs$my_data; EVT_optimizer_outputs$my_data = NULL





###############################################
## plot result
###############################################
EVT_PZC_plot(my_data, min(alphas_extreme), sub_idx = 75e3:80e3)


###############################################
## print the result
###############################################
print("MSE Results")
NZ_MSE = function(x, r1, r2, r1true, r2true, idx = NULL, alpha, scale = 1e3) {
  # for alpha close to 1 and right tail VaR !!
  if (is.null(idx)) idx = 1:length(x)
  ee = data.frame(
    MAE_VaR = abs(r1 - r1true),
    MSE_VaR = (r1 - r1true)^2,
    NZ_2.19 = scale * ( (1 - alpha) * r1 + ifelse(x > r1, x - r1, 0) ),
    NZ_2.20 = scale * ( (1 - alpha) * log(r1) + ifelse(x > r1, log(x/r1), 0) ),
    MAE_ES = abs(r2 - r2true),
    MSE_ES = (r2 - r2true)^2,
    NZ_2.23 = scale * ( (x > r1) * 0.5 * (x - r1) / sqrt(r2) + (1 - alpha) * 0.5 * (r1 + r2) / sqrt(r2) ),
    NZ_2.24 = scale * ( (x > r1) * (x - r1) / r2 + (1 - alpha) * (r1 / r2 - 1 + log(r2)) )
  )
  return(ee)
}
out1 = NZ_MSE(-my_data$y, 
              -my_data$EVT_VaR_alpha0.001, -my_data$EVT_EL_alpha0.001,
              -my_data$VaRtrue_alpha0.001, -my_data$ELtrue_alpha0.001, 
              idx = NULL, 1 - min(alphas_extreme))
out2 = NZ_MSE(-my_data$y, 
              -my_data$PZC_VaR0Neld_alpha0.001, -my_data$PZC_EL0Neld_alpha0.001,
              -my_data$VaRtrue_alpha0.001, -my_data$ELtrue_alpha0.001, 
              idx = NULL, 1 - min(alphas_extreme))
out_tab = rbind( colMeans(out1), colMeans(out2) )
out_tab = rbind(out_tab, out_tab)
for (i1 in 1:ncol(out_tab)) {
  aid1 = t.test(out1[ , i1] - out2[ , i1])
  out_tab[3, i1] = aid1$statistic
  out_tab[4, i1] = aid1$p.value
}
rownames(out_tab) = c('EVT', 'PZC', 'DM', 'DM pval')
print(round(out_tab,3))



###############################################
## subselect figure data in separate dataframe
###############################################
store_filename = paste0("simulation_PZC_EVT_200k.csv")
save_PZC_EVT_frame(my_data, store_filename = store_filename,
                   select_names = c(
                     'dates', 'y', 'df_inv', 'VaRtrue_alpha0.05',
                     'VaRtrue_alpha0.001', 'ELtrue_alpha0.001',
                     'PZC_VaR0Neld_alpha0.001', 'PZC_EL0Neld_alpha0.001',
                     'EVT_tau',
                     # 'EVT_ft', 'EVT_bandL', 'EVT_bandU',
                     'EVT_VaR_alpha0.001', 'EVT_EL_alpha0.001',
                     NULL),
                   new_names = c(
                     'dates', 'y', 'df_inv', 'tau_true',
                     'VaRtrue_alpha0.001', 'ELtrue_alpha0.001',
                     'PZC_VaR_alpha0.001', 'PZC_EL_alpha0.001',
                     'EVT_tau',
                     # 'EVT_ft', 'EVT_bandL', 'EVT_bandU',
                     'EVT_VaR_alpha0.001', 'EVT_EL_alpha0.001',
                     NULL))

