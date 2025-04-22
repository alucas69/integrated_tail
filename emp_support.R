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
    file_name = file.path(folder_name, "FX/EURUSD.csv")
    if (!file.exists(file_name)) throw_error()
    my_data = read.csv(file_name)
    names(my_data) = c("dates", "rate")
    my_data$date = as.Date(my_data$dates, "%m/%d/%y")
    my_data = my_data[order(my_data$date), ]
    my_data$ret = -c(NA, 100 * diff(log(as.numeric(my_data$rate))))
    my_data = subset(my_data, (!is.na(ret)) & (date >= "1999-01-04") & (date <= "2023-12-01"))
    
    ## set the starting values
    x0 = NULL
    
  } else   if (asset_name == "RUB") {
    ## read the data
    file_name = file.path(folder_name, "FX/USDRUB.csv")
    if (!file.exists(file_name)) throw_error()
    my_data = read.csv(file_name)
    names(my_data) = c("dates", "rate")
    my_data$date = as.Date(my_data$dates, "%m/%d/%y")
    my_data = my_data[order(my_data$date), ]
    my_data$ret = c(NA, 100 * diff(log(as.numeric(my_data$rate))))
    my_data = subset(my_data, (!is.na(ret)) & (date >= "1996-01-02") & (date <= "2023-09-06"))
    
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


