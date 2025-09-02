##########################################
##########################################
## IMPORTANT: WE CONSIDER RISK MEASURES
##            IN THE RIGHT TAIL.
##            FOR THE LEFT TAIL, FLIP
##            THE RETURNS
##########################################
##########################################



###################
## Set up functions
###################
library(ggplot2)
library(ggpubr)
library(Rcpp)
library(RcppArmadillo)
library(forecast)
require(readxl)
require(writexl)
library(lubridate)
library(rugarch)

source('PZC_simulation_aid.R')
source('Integrated_EVT_filter.R')
source('TVGPD_filter.R')
sourceCpp('PZC_simulation.cpp')
sourceCpp('Integrated_EVT_filter.cpp')


#####################################
## Set tail percentage for thresholds
## and extreme tail percentages for
## VaR and EL
#####################################
verbosity_level = 0
NZ_table_out = total_table = NULL
## IMPORTANT: the alpha loop should be inside the asset loop, as the GARCH 
##            is NOT re-estimated for different alphas
alphas = 0.1/c(1,2,4); for (asset_name in c("XRP", "BTC", "ETH")) { NZ_table_out = NULL; GARCH_estimated = FALSE; for (alpha_tail in alphas) {alpha_extreme = alpha_tail / 10
# alphas = 0.1/c(1,2,4); for (asset_name in c("EUR", "RUB", "BTC", "ETH")) { NZ_table_out = NULL; GARCH_estimated = FALSE; for (alpha_tail in alphas) {alpha_extreme = alpha_tail / 10
# alpha_tail = 0.025
# alpha_extreme = alpha_extremes = 0.0025
# asset_name = c("EUR", "RUB", "BTC", "ETH")[1]
# NZ_table = NULL
# {{
  
###############################
## set model configurations PZC
## PZC = Patton/Ziegel/Chen
##       J.Econometrics
## IMPORTANT: PZC is orignally
##     for the left tail. We
##     recast it to the right
##     tail by implementing it
##     for the flipped 
##     empirical series
###############################
s_PZC_OPTIONS = list(
  diagonal_A = FALSE,        # boolean; TRUE: use diagonal A matrix 
  diagonal_B = TRUE,         # boolean; TRUE: use diagonal B matrix 
  ALPHA_EXTREME = alpha_extreme,  # double, tail percentage for the thresholds tau
  rw_specification = FALSE   # boolean; TRUE: use unit B matrix; overrides diagonal_B
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
  BANDS_NRSIMS = 1000,        # integer, number of simulations for the confidence bands
  BANDS_PCT = 0.75,          # double, confidence level for the confidence bands
  TAU_TAIL_PCT = alpha_tail, # double, tail percentage for the thresholds tau
  ALPHA_EXTREME = alpha_extreme, # double, tail percentage for the thresholds tau
  EXTERNAL_TAU = NULL        # can be overwritten if external thresholds available
)


##################################
## set model configurations TVGPD
##################################
s_TVGPD_OPTIONS = list(
  TAU_TAIL_PCT = alpha_tail, # double, tail percentage for the thresholds tau
  ALPHA_EXTREME = alpha_extreme, # double, tail percentage for the thresholds tau
  EXTERNAL_TAU = NULL        # can be overwritten if external thresholds available
)


#############################
## Read data
#############################
if (asset_name == "EUR") {
  ## EURUSD
  flip_left_to_right = FALSE
  my_data = read_empirical_data('../../Data/updated data/EURUSD.xlsx', 
                                price_column_name = "EURO", date_column_name = "Date",
                                prices_are_already_returns = FALSE,
                                flip_return = flip_left_to_right, 
                                begin_date = "1999-01-04")
  estimation_evaluation_dates = cbind(
    "1999-01-04", paste0(2011 + 1:14, "-12-31"),
    "1999-01-04", paste0(2012 + 1:14, "-12-31")
  )
} else if (asset_name == "RUB") {
  ## USDRUB
  flip_left_to_right = FALSE
  my_data = read_empirical_data('../../Data/updated data/USDRUB.xlsx', 
                                price_column_name = "RUSSIAN ROUBLE", date_column_name = "Date",
                                prices_are_already_returns = FALSE,
                                flip_return = flip_left_to_right, begin_date = "1999-01-01")
  estimation_evaluation_dates = cbind(
    "1999-01-04", paste0(2011 + 1:14, "-12-31"),
    "1999-01-04", paste0(2012 + 1:14, "-12-31")
  )
} else if (asset_name %in% c("BTC", "ETH", "XRP")) {
  ## CRYPTO-USD
  flip_left_to_right = TRUE
  my_data = read_empirical_data(paste0('../../Data/updated data/Bitfinex_', asset_name, 'USD_1h.csv'), 
                                price_column_name = "close", date_column_name = "date",
                                prices_are_already_returns = FALSE, reverse = TRUE, 
                                date_POSX_format = "%Y-%m-%d %H:%M:%S", tz = "UTC",
                                flip_return = flip_left_to_right)
  estimation_evaluation_dates = cbind(
    "2018-01-05", paste0(2020 + 1:5, "-12-31"),
    "2018-01-05", paste0(2021 + 1:5, "-12-31")
  )
} else {
  error(paste0("ERROR: WRONG ASSET NAME ", asset_name))
}
  


###############################################
## loop over the different estimation and
## evaluation periods
###############################################
initial_oos_run = TRUE
for (sample_count in 1:nrow(estimation_evaluation_dates)) {

  ###############################################
  ## set the sample
  ###############################################
  in_sample_idx = get_index_from_date_vector(estimation_evaluation_dates[sample_count, 1], my_data$dates, start = TRUE) :
    get_index_from_date_vector(estimation_evaluation_dates[sample_count, 2], my_data$dates, start = FALSE)
  out_of_sample_idx = get_index_from_date_vector(estimation_evaluation_dates[sample_count, 3], my_data$dates, start = TRUE) :
    get_index_from_date_vector(estimation_evaluation_dates[sample_count, 4], my_data$dates, start = FALSE)
  if (verbosity_level > -1) print(paste0(
    "(", asset_name, ", ", alpha_extreme, "): ",
    "In-sample: ", paste0(as.Date(my_data$dates[range(in_sample_idx)]), collapse = "/"),
    ";  Out-of-sample: ", paste0(as.Date(my_data$dates[range(out_of_sample_idx)]), collapse = "/")
  ))
  


  ###############################################
  ## estimate PZC specification for extreme alpha
  ## IMPORTANT: this gives the RIGHT tail VaR
  ##            and ES
  ###############################################
  PZC_optimizer_outputs = 
    PZC_optimize(my_data, filter = PZC_filter_cpp, 
                 PZC_OPTIONS = s_PZC_OPTIONS, verbosity = verbosity_level,
                 in_sample_idx = in_sample_idx,
                 out_of_sample_idx = out_of_sample_idx)
  my_data = PZC_optimizer_outputs$my_data; PZC_optimizer_outputs$my_data = NULL
  
  
  
  
  ###############################################
  ## estimate EVT specification for extreme alpha
  ## using PZC taus
  ###############################################
  EVT_optimizer_outputs = EVT_optimize(my_data, EVT_OPTIONS = s_EVT_OPTIONS, 
                                       PZC_OPTIONS = s_PZC_OPTIONS,
                                       verbosity = verbosity_level, 
                                       in_sample_idx = in_sample_idx,
                                       out_of_sample_idx = out_of_sample_idx)
  my_data = EVT_optimizer_outputs$my_data; EVT_optimizer_outputs$my_data = NULL
  
  
  
  
  ###############################################
  ## estimate TVGPD specification for extreme 
  ## alpha using PZC taus
  ###############################################
  s_TVGPD_OPTIONS$EXTERNAL_TAU = my_data$EVT_tau
  my_data = TVGPD_optimize(my_data, TVGPD_OPTIONS = s_TVGPD_OPTIONS, 
                           PZC_OPTIONS = s_PZC_OPTIONS, verbosity = 0,
                           in_sample_idx = in_sample_idx,
                           out_of_sample_idx = out_of_sample_idx)

  
  
  
  ###############################################
  ## estimate GARCH specifications for extreme
  ## alpha
  ###############################################
  ## economize on GARCH re-estimation for different alpha_extreme does not yet work
  ## in current nesting of loops
  if (!GARCH_estimated) {GARCH_optimizer_outputs = NULL; GARCH_estimated = FALSE}
  GARCH_optimizer_outputs = GARCH_optimizer(my_data, verbosity = verbosity_level,
                                            in_sample_idx = in_sample_idx,
                                            out_of_sample_idx = out_of_sample_idx,
                                            GARCH_optimizer_outputs = GARCH_optimizer_outputs)
  my_data = GARCH_optimizer_outputs$my_data; GARCH_optimizer_outputs = GARCH_optimizer_outputs$GARCH_optimizer_outputs
  
  
  
  ###############################################
  ## store the estimation results
  ###############################################
  if (max(out_of_sample_idx) > max(in_sample_idx)) {
    if (initial_oos_run) {my_data_oos = my_data + NA; initial_oos_run = FALSE}
    aid_idx = (max(in_sample_idx)+1):max(out_of_sample_idx)
    my_data_oos[ aid_idx, ] = my_data[aid_idx, ]
    for (nm in names(my_data)) my_data_oos[ aid_idx, nm] = my_data[aid_idx, nm]
  }
  
}

###############################################
## store file with IS and OOS filter values
###############################################
writexl::write_xlsx(my_data, path = paste0("results/", asset_name, " IS ", alpha_extreme, ".xlsx"))
writexl::write_xlsx(my_data_oos, path = paste0("results/", asset_name, " OOS ", alpha_extreme, ".xlsx"))


###############################################
## plot result (FULL SAMPLE only)
###############################################
if (max(in_sample_idx - 1:nrow(my_data)) < 1e-4) {
  flip_back_to_left = ifelse(flip_left_to_right, -1, 1)
  ## plot data and VaR levels EVT & PZC
  gg =
    ggplot(
      data = data.frame(
        x = rep(my_data$dates, 6),
        y = flip_back_to_left * c(my_data$PZC_VaR, my_data$EVT_VaR, my_data$EVT_tau, my_data$TVGPD_VaR, my_data$GARCH_VaR, my_data$GARCHt_VaR),
        VaR = rep(c("PZC", "EVT", "tau", "TVGPD", "GARCH", "GARCHt"), each = nrow(my_data))
      ),
      aes(x = x, y = y, color = VaR)
    ) +
    geom_line() +
    geom_line(data = my_data, aes(x = dates, y = flip_back_to_left * y), color = 'gray')
  plot(gg)
  ## plot data and VaR levels EVT & PZC
  gg =
    ggplot(
      data = data.frame(
        x = rep(my_data$dates, 5),
        y = flip_back_to_left * c(my_data$PZC_ES, my_data$EVT_ES, my_data$TVGPD_ES,  my_data$GARCH_ES, my_data$GARCHt_VaR),
        ES = rep(c("PZC", "EVT", "TVGPD", "GARCH", "GARCHt"), each = nrow(my_data))
      ),
      aes(x = x, y = y, color = ES)
    ) +
    geom_line() +
    geom_line(data = my_data, aes(x = dates, y = flip_back_to_left * y), color = 'gray')
  plot(gg)
  ## plot ft plus band
  gg = 
    ggplot(
      data = data.frame(
        x = rep(my_data$dates, 3), 
        y = c(my_data$EVT_ft, my_data$EVT_bandL, my_data$EVT_bandU),
        ES = rep(c("f(t)", "low", "high"), each = nrow(my_data))
      ),
      aes(x = x, y = y, color = ES)
    ) +
    geom_line()
  plot(gg)
}


###############################################
## print the NZ result
###############################################
NZ_table_out = print_NZ_table(in_sample_data, out_of_sample_data, alpha_extreme, 
                              asset_name, NZ_table_in = NZ_table_out, 
                              with_print = (alpha_tail == alphas[length(alphas)]))
if (alpha_tail == alphas[length(alphas)]) {
  outfile = file(paste0("results/", asset_name, alpha_extreme, ".txt"), "w"); writeLines(NZ_table_out, con = outfile); close(outfile)
  total_table = c(total_table, rep('', 3), NZ_table_out)
}

}}

###############################################
## print the total result and write to file
###############################################
writeLines(total_table)
outfile = file(paste0("results/total_result.txt"), "w"); writeLines(total_table, con = outfile); close(outfile)
