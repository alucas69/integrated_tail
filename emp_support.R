#######################################################
## load the data
#######################################################
load_data <- function(
    asset_name, folder_name = "input data") {
  
  org_call = deparse(sys.calls()[[sys.nframe()]])
  throw_error = function(extra = "") stop(paste0("Error in ", org_call, ": file does not exist; ", extra))
  
  ## hard coded file names for which the file works, 
  ## inclusive of starting values
  if (asset_name == "EUR") {
    ## read the data
    file_name = file.path(folder_name, "FX/EURUSD_042025.csv")
    if (!file.exists(file_name)) throw_error()
    my_data = read.csv(file_name)
    names(my_data) = c("dates", "rate")
    my_data$date = as.Date(my_data$dates, "%m/%d/%Y")
    my_data = my_data[order(my_data$date), ]
    my_data$ret = -c(NA, 100 * diff(log(as.numeric(my_data$rate))))
    my_data = subset(my_data, (!is.na(ret)) & (date >= "1990-01-04")) 
    plot(my_data$date, my_data$ret, type = "l")
    
    ## set the starting values
    x0 = NULL
    
  } else   if (asset_name == "RUB") {
    ## read the data
    file_name = file.path(folder_name, "FX/USDRUB_042025.csv")
    if (!file.exists(file_name)) throw_error()
    my_data = read.csv(file_name)
    names(my_data) = c("dates", "rate")
    my_data$date = as.Date(my_data$dates, "%m/%d/%Y")
    my_data = my_data[order(my_data$date), ]
    my_data$ret = c(NA, 100 * diff(log(as.numeric(my_data$rate))))
    my_data = subset(my_data, (!is.na(ret)) & (date >= "1996-01-02")) 
    
    ## set the starting values
    x0 = NULL
    
  } else   if (asset_name == "BTC") {
    ## read the data
    file_name = file.path(folder_name, "Crypto/Bitfinex_BTCUSD_1h.csv")
    if (!file.exists(file_name)) throw_error()
    my_data = read.csv(file_name)[ , c("date", "close")]
    names(my_data) = c("dates", "rate")
    my_data$date = as.POSIXct(my_data$dates)
    my_data = my_data[nrow(my_data):1, ]
    my_data$ret = c(NA, 100 * diff(log(as.numeric(my_data$rate))))
    my_data = subset(my_data, (!is.na(ret)) & (date >= "2018-05-15") & (date <= "2023-01-25"))
    
    ## set the starting values
    x0 = NULL
    
  } else   if (asset_name == "ETH") {
    ## read the data
    file_name = file.path(folder_name, "Crypto/Bitfinex_ETHUSD_1h.csv")
    if (!file.exists(file_name)) throw_error()
    my_data = read.csv(file_name)[ , c("date", "close")]
    names(my_data) = c("dates", "rate")
    my_data$date = as.POSIXct(my_data$dates)
    my_data = my_data[nrow(my_data):1, ]
    my_data$ret = c(NA, 100 * diff(log(as.numeric(my_data$rate))))
    my_data = subset(my_data, (!is.na(ret)) & (date >= "2018-05-15") & (date <= "2023-01-28"))
    
    ## set the starting values
    x0 = NULL
    
  } else throw_error(paste0("(wrong asset ", asset_name, ")"))
  
  
  ## return
  return(list(data = my_data, x0 = x0))
}





#######################################################
## Nolde Ziegel loss functions for RIGHT tail
#######################################################
Nolde_Ziegel_2_19 = function(x, r, alfa) return( 
  (1 - alfa) * r + ifelse(x > r, x - r, 0))
Nolde_Ziegel_2_20 = function(x, r, alfa) return( 
  ifelse(r > 0, (1 - alfa) * log(r) + ifelse(x > r, log(x / r), 0), 0) )
Nolde_Ziegel_2_23 = function(x, r1, r2, alfa) return( 
  0.5 * (1 - alfa) * (r1 + r2) / sqrt(r2) + 0.5 * ifelse(x > r1, (x - r1)/sqrt(r2), 0))
Nolde_Ziegel_2_24 = function(x, r1, r2, alfa) return( 
  (1 - alfa) * ((r1 / r2) - 1 + log(r2)) + ifelse(x > r1, (x - r1)/r2, 0))
emp_evaluate_performance = function(my_data, alpha_extreme, corename_VaR = "EVT_VaR_alpha", 
                                    corename_ES = "EVT_EL_alpha", tailindicator = 1) {
  ydata = tailindicator * my_data$y
  EVT_VaR = tailindicator * my_data[ ,  paste0(corename_VaR, alphas_extreme)]
  EVT_ES = tailindicator * my_data[ ,  paste0(corename_ES, alphas_extreme)]
  PZC_VaR = tailindicator * my_data$VaR
  PZC_ES = tailindicator * my_data$EL
  aid1 = as.data.frame(matrix(c(
    mean(ydata >= EVT_VaR),
    mean(Nolde_Ziegel_2_19(ydata, EVT_VaR, alphas_extreme)),
    mean(Nolde_Ziegel_2_20(ydata, EVT_VaR, alphas_extreme)),
    mean(Nolde_Ziegel_2_23(ydata, EVT_VaR, EVT_ES, alphas_extreme)),
    mean(Nolde_Ziegel_2_24(ydata, EVT_VaR, EVT_ES, alphas_extreme)),
    ##
    mean(ydata >= PZC_VaR),
    mean(Nolde_Ziegel_2_19(ydata, PZC_VaR, alphas_extreme)),
    mean(Nolde_Ziegel_2_20(ydata, PZC_VaR, alphas_extreme)),
    mean(Nolde_Ziegel_2_23(ydata, PZC_VaR, PZC_ES, alphas_extreme)),
    mean(Nolde_Ziegel_2_24(ydata, PZC_VaR, PZC_ES, alphas_extreme))
  ), ncol = 2))
  names(aid1) = c("EVT", "PZC")
  rownames(aid1) = c("Hit_rate", "NZ(2.19)", "NZ(2.21)", "NZ(2.23)", "NZ(2.24)")
  return(aid1)
}



