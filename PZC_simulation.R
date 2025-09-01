###################
## Set up functions
###################
library(ggplot2)
library(ggpubr)
library(Rcpp)
library(RcppArmadillo)
library(forecast)

source('PZC_simulation_aid.R')
source('Integrated_EVT_filter.R')
sourceCpp('PZC_simulation.cpp')
sourceCpp('Integrated_EVT_filter.cpp')


#####################################
## Set tail percentage for thresholds
## and extreme tail percentages for
## VaR and EL
#####################################
# alpha_tail = 0.05
# alphas_extreme = c(0.05, 1e-3)
alpha_tail = 0.025
alphas_extreme = c(0.025, 2.5e-3)
asset_name = "simulation"

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
  set_seed = NULL            # NULL or integer: if NULL, do not set the seed when generating data
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
  BANDS_NRSIMS = 100,       # integer, number of simulations for the confidence bands
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
    dim_n, alphas, df_min = 2.5, df_max = 10, nr_nu_cycles = 1, nr_sigma_cycles = 1,
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
  df_sim_inv  = 1/df_max + 0.5 * (1/df_min - 1/df_max) * (1 - cos(nr_nu_cycles * (1:dim_n) / dim_n * 2 * pi))
  sigma_sim = -(1 + 4 * abs(((1:dim_n) %% (dim_n/nr_sigma_cycles))/(dim_n/nr_sigma_cycles) - 0.5)) / qt(0.05, df = 1 / df_sim_inv) # rep(1, dim_n)
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
dim_n = 40000
sub_frame = round(dim_n * c(0.75, 0.8))
my_data = generate_data_student_t_v1(dim_n, alphas_extreme, 
                                     nr_nu_cycles = 3,
                                     nr_sigma_cycles = 1,
                                     set_seed = s_PZC_OPTIONS$set_seed)
#############################
## Replace by empirical data
#############################
# ## EURUSD
# my_data2 = readxl::read_xlsx('../../Data/updated data/EURUSD.xlsx'); my_data2 = subset(my_data2, Date > "1999-01-01")
# # my_data2 = readxl::read_xlsx('../../Data/updated data/USDRUB.xlsx'); my_data2 = subset(my_data2, Date > "1999-01-01")
# names(my_data2)[2] = 'close'
# my_data2$Return = c(NA, 100 * diff(log(my_data2$close)))
# my_data2 = my_data2[complete.cases(my_data2), ]
# my_data = data.frame(
#   dates = 1:nrow(my_data2), y = my_data2$Return, df_inv = 0.1,
#   VaRtrue_alpha0.1 = 0, ELtrue_alpha0.1 = 0,
#   VaRtrue_alpha0.01 = 0, ELtrue_alpha0.01 = 0
# )
# dim_n = nrow(my_data); sub_frame = round(dim_n * c(0.75, 0.8))
# zz = 1; alpha_tail = 0.1 / zz; alphas_extreme = c(alpha_tail, alpha_tail / 10)
# s_PZC_OPTIONS$ALPHAS_EXTREME = alphas_extreme; s_EVT_OPTIONS$TAU_TAIL_PCT = alpha_tail; s_EVT_OPTIONS$ALPHAS_EXTREME = alphas_extreme
# names(my_data) = c('dates', 'y', 'df_inv', paste0('VaRtrue_alpha', alphas_extreme), paste0('ELtrue_alpha', alphas_extreme))
## BTCUSD
my_data2 = read.csv('../../Data/updated data/Bitfinex_BTCUSD_1h.csv')
# my_data2 = read.csv('../../Data/updated data/Bitfinex_ETHUSD_1h.csv')
my_data2 = my_data2[nrow(my_data2):1, ]
my_data2$Return = c(NA, -100 * diff(log(my_data2$close)))
my_data2 = my_data2[complete.cases(my_data2), ]
my_data = data.frame(
  dates = 1:nrow(my_data2), y = -my_data2$Return, df_inv = 0.1,
  VaRtrue_alpha0.025 = 0, ELtrue_alpha0.025 = 0, VaRtrue_alpha0.0025 = 0, ELtrue_alpha0.0025 = 0
)
dim_n = nrow(my_data); sub_frame = round(dim_n * c(0.75, 0.8))
zz = 4; alpha_tail = 0.1 / zz; alphas_extreme = c(alpha_tail, alpha_tail / 10)
s_PZC_OPTIONS$ALPHAS_EXTREME = alphas_extreme; s_EVT_OPTIONS$TAU_TAIL_PCT = alpha_tail; s_EVT_OPTIONS$ALPHAS_EXTREME = alphas_extreme
names(my_data) = c('dates', 'y', 'df_inv', paste0('VaRtrue_alpha', alphas_extreme), paste0('ELtrue_alpha', alphas_extreme))



###############################################
## estimate PZC specification for extreme alpha
###############################################
in_sample_idx = 1:(nrow(my_data) - 10000)
out_of_sample_idx = 1:nrow(my_data)
PZC_optimizer_outputs = 
  PZC_optimize(my_data, filter = PZC_filter_cpp, 
               PZC_OPTIONS = s_PZC_OPTIONS, verbosity = 0,
               in_sample_idx = in_sample_idx,
               out_of_sample_idx = out_of_sample_idx)
my_data = PZC_optimizer_outputs$my_data; PZC_optimizer_outputs$my_data = NULL


# # ## extra plots for empirical data
# ytmp = cbind(my_data_test[ , paste0('VaRtrue_alpha', alphas_extreme[1])],
#              my_data_test[ , paste0('PZC_VaR0Neld_alpha', alphas_extreme[1])])[out_of_sample_idx, ]
# plot(ytmp[ , 1], type = "l", ylim = range(ytmp)); lines(ytmp[ , 2], col = "red")
# ytmp = cbind(my_data_test[ , paste0('VaRtrue_alpha', alphas_extreme[2])],
#              my_data_test[ , paste0('PZC_VaR0Neld_alpha', alphas_extreme[2])])[out_of_sample_idx, ]
# plot(ytmp[ , 1], type = "l", ylim = range(ytmp)); lines(ytmp[ , 2], col = "red")
# stop(0)


###############################################
## plot PZC results
###############################################
# plot_VaR(my_data, alphas = s_PZC_OPTIONS$ALPHAS_EXTREME, 
#          sub_frame_idx = sub_frame,
#          gg_plt_extras = ylim(-20,0),
#          make_plot = FALSE)





###############################################
## estimate EVT specification for extreme alpha
###############################################
EVT_optimizer_outputs = EVT_optimize(my_data, EVT_OPTIONS = s_EVT_OPTIONS, verbosity = 0)
my_data = EVT_optimizer_outputs$my_data; EVT_optimizer_outputs$my_data = NULL



stop(0)

###############################################
## plot result
###############################################
sub_frame = sub_frame[1]:sub_frame[2]
# EVT_PZC_plot(my_data, min(alphas_extreme), sub_idx = sub_frame)


###############################################
## print the result
###############################################
print(paste0(asset_name, "MSE Results (tau alpha = ", alpha_tail, ")"))
out_tab = NZ_table(-my_data, data_name = "y", evaluation_idx = out_of_sample_idx,
                   VaR_names = c(paste0('EVT_VaR_alpha', min(alphas_extreme)), paste0('PZC_VaR0Neld_alpha', min(alphas_extreme))),
                   ES_names = c(paste0('EVT_EL_alpha', min(alphas_extreme)), paste0('PZC_EL0Neld_alpha', min(alphas_extreme))),
                   alphas = 1 - min(alphas_extreme), method_names = c("EVT", "PZC"))
print(round(out_tab,3))
stop(0)


###############################################
## subselect figure data in separate dataframe
###############################################
# store_filename = paste0("simulation_PZC_EVT_22k.csv")
# save_PZC_EVT_frame(my_data, store_filename = store_filename,
#                    select_names = c(
#                      'dates', 'y', 'df_inv', 'VaRtrue_alpha0.05',
#                      'VaRtrue_alpha0.001', 'ELtrue_alpha0.001',
#                      'PZC_VaR0Neld_alpha0.001', 'PZC_EL0Neld_alpha0.001',
#                      'EVT_tau',
#                      # 'EVT_ft', 'EVT_bandL', 'EVT_bandU',
#                      'EVT_VaR_alpha0.001', 'EVT_EL_alpha0.001',
#                      NULL),
#                    new_names = c(
#                      'dates', 'y', 'df_inv', 'tau_true',
#                      'VaRtrue_alpha0.001', 'ELtrue_alpha0.001',
#                      'PZC_VaR_alpha0.001', 'PZC_EL_alpha0.001',
#                      'EVT_tau',
#                      # 'EVT_ft', 'EVT_bandL', 'EVT_bandU',
#                      'EVT_VaR_alpha0.001', 'EVT_EL_alpha0.001',
#                      NULL))

