###################
## Set up functions
###################
library(ggplot2)
library(ggpubr)
library(Rcpp)
library(RcppArmadillo)
library(foreach)
library(doParallel)
library(doRNG)
library(optimParallel)
nr_free_cores = max(1, detectCores() - 1)
registerDoParallel(cores = nr_free_cores)

source('emp_support.R')
source('emp_integrated_evt_filter.R')
source('emp_pzc_support.R')
sourceCpp('PZC_simulation.cpp')
sourceCpp('Integrated_EVT_filter.cpp')


#####################################
## Select the asset and tail
#####################################
asset_nr = 1
asset_name = c("EUR", "RUB", "BTC", "ETH")[asset_nr]
lower_tail = c(F,F,T,T)[asset_nr]

# my_sub_frame = list(c("2004-10-01", "2012-12-31"), c("2010-01-01", "2024-12-31"), c("2021-10-01", "2023-01-31"), c("2021-10-01", "2023-01-31"))[[asset_nr]]
my_frame_limits = list( full_frame_min = c(-4.2,-14,-22,-25), full_frame_max = c(6,43,22,22),
                        sub_frame_min = c(-4.2, -18, -8, -15), sub_frame_max = c(6,30,15,20))
my_data_y_min = my_frame_limits[[ifelse(!is.null(my_frame_limits), "full_frame_min", "sub_frame_min")]][asset_nr]
my_data_y_max = my_frame_limits[[ifelse(!is.null(my_frame_limits), "full_frame_max", "sub_frame_max")]][asset_nr]


#####################################
## Set tail percentage for thresholds
## and extreme tail percentages for
## VaR and EL
#####################################
alpha_tail = 0.10
alphas_extreme = 0.01


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
  BANDS_NRSIMS = 0,       # integer, number of simulations for the confidence bands
  BANDS_PCT = 0.75,          # double, confidence level for the confidence bands
  TAU_TAIL_PCT = alpha_tail, # double, tail percentage for the thresholds tau
  ALPHAS_EXTREME = alphas_extreme,  # double, tail percentage for the thresholds tau
  EXTERNAL_TAU = NULL,       # can be overwritten if external thresholds available
  set_seed = s_PZC_OPTIONS$set_seed # NULL or integer: if NULL, do not set the seed when generating data
)



################
## Load data
################
my_data = load_data(asset_name)
x0 = my_data$x0
my_data = my_data$data
my_data$y = my_data$ret

# mydata = readxl::read_xlsx("input data/estimation_results.xlsx")
# my_data = subset(mydata, select = c("Return"))
# my_data$y = my_data$Return

###############################################
## estimate PZC specification for tail alpha
###############################################
aid1 = s_PZC_OPTIONS$ALPHAS_EXTREME; s_PZC_OPTIONS$ALPHAS_EXTREME = alpha_tail
my_data = PZC_optimize(
  my_data, filter = PZC_filter_cpp, PZC_OPTIONS = s_PZC_OPTIONS, verbosity = 0,
  only_store_last = TRUE, lower_tail = lower_tail)$my_data
my_data$tau = my_data$VaR
s_PZC_OPTIONS$ALPHAS_EXTREME = aid1


###############################################
## use the PZC VaR results as thresholds
###############################################
s_EVT_OPTIONS$EXTERNAL_TAU = my_data$tau
###############################################
## estimate EVT specification for extreme alpha
## using the PZC thesholds
###############################################
my_data = estimate_full_model(
  my_data, EVT_OPTIONS = s_EVT_OPTIONS, verbosity = 0, 
  lower_tail = lower_tail)$my_data




###############################################
## estimate PZC specification for extreme alpha
###############################################
my_data = PZC_optimize(
  my_data, filter = PZC_filter_cpp, 
  PZC_OPTIONS = s_PZC_OPTIONS, verbosity = 0,
  only_store_last = TRUE, lower_tail = lower_tail)$my_data




###############################################
## plot result
###############################################
my_idx = 1:nrow(my_data)
# if (!is.null(my_sub_frame)) my_idx = which((my_data$date >= my_sub_frame[1]) & (my_data$date <= my_sub_frame[2]))
gg_evt = emp_plot_VaR(
  my_data[my_idx, ], # sub_time = grid,
  series_names = c("y", "tau", paste0(c("EVT_VaR_alpha", "EVT_EL_alpha"), alphas_extreme)),
  series_plt_names = c("y", paste0("tau(a=", alpha_tail, ")"), 
                       paste0("VaR(a=", alphas_extreme, ")"),
                       paste0("EL(a=", alphas_extreme, ")")),
  gg_plt_extras = coord_cartesian(ylim = c(my_data_y_min,my_data_y_max)),
  make_plot = FALSE)
gg_pzc = emp_plot_VaR(
  my_data[my_idx, ], 
  series_names = c("y", "tau", paste0(c("VaR", "EL"))),
  series_plt_names = c("y", paste0("tau(a=", alpha_tail, ")"), 
                       paste0("VaR(a=", alphas_extreme, ")"),
                       paste0("EL(a=", alphas_extreme, ")")),
  gg_plt_extras = coord_cartesian(ylim = c(my_data_y_min,my_data_y_max)),
  make_plot = FALSE)
plot(ggarrange(gg_evt, gg_pzc, nrow = 1, ncol = 2))


###############################################
## backtest
###############################################
aid1 = ifelse(lower_tail, -1, 1)
print("Backtest")
print(sprintf("alpha = %2.4f", alphas_extreme))
print(sprintf("EVT   = %2.4f", mean(aid1 * my_data$y >= aid1 * my_data[ ,  paste0("EVT_VaR_alpha", alphas_extreme)])))
print(sprintf("PZC   = %2.4f", mean(aid1 * my_data$y >= aid1 * my_data$VaR)))
