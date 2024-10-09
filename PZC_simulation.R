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
    dim_n, alphas, df_min = 2.5, df_max = 15, nr_cycles = 3, 
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
  sigma_sim = rep(1, dim_n)
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
dim_n = 100000
my_data = generate_data_student_t_v1(dim_n, alphas_extreme, 
                                     nr_cycles = 3,
                                     set_seed = s_PZC_OPTIONS$set_seed)


###############################################
## estimate PZC specification for extreme alpha
###############################################
my_data = PZC_optimize(my_data, filter = PZC_filter_cpp, 
                       PZC_OPTIONS = s_PZC_OPTIONS, verbosity = 0)
###############################################
## plot PZC results
###############################################
store_filename = paste0("simulation_LFZ0.csv")
plot_VaR(my_data, alphas = s_PZC_OPTIONS$ALPHAS_EXTREME, 
         # store_filename = store_filename,
         # sub_frame_idx = 40000:51000, 
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





