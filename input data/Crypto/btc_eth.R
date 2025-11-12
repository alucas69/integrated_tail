library(binancer)
library(ggplot2)
library(lubridate)
library(data.table)

# all_coins = binance_coins_prices()[c(1,11), ]
# 
# binance_coins()
# 
# 
# klines <- rbindlist(lapply(
#   c('BTCUSD'),
#   binance_klines,
#   interval = '1h',
#   limit = 2))
#   # start_time = '2022-01-26 00:00:00',
#   # end_time = '2023-01-26 00:00:00')))
# 
# ggplot(klines, aes(open_time)) +
#   geom_linerange(aes(ymin = open, ymax = close, color = close < open), size = 2) +
#   geom_errorbar(aes(ymin = low, ymax = high), size = 0.25) +
#   theme_bw() + theme('legend.position' = 'none') + xlab('') +
#   ggtitle(paste('Last Updated:', Sys.time())) +
#   scale_color_manual(values = c('#1a9850', '#d73027')) +
#   facet_wrap(~symbol, scales = 'free', nrow = 2)


# Reference: https://www.cryptodatadownload.com/blog/how-to-download-coinbase-price-data-R.html
#Import libraries needed
# Use install.packages("glue") if needed
library(jsonlite)
library(glue)
# Create a function to retrieve daily data
retreive_daily_data <- function(pair) {
  url = glue("https://api.pro.coinbase.com/products/{pair}/candles?granularity=86400")
  columnNames <- c('unix', 'low', 'high', 'open', 'close', glue('{pair} volume'))
  mydata <- fromJSON(url)
  df <- as.data.frame(mydata)
  colnames(df) <- columnNames  # rename the columns
  
  write.csv(df, file = filename)
}

newPair <- "BTC-USD"
fileName <- glue("dailyData{newPair}.csv")
runFunc <- retreive_daily_data(newPair) #, filename = fileName)
runFunc


old_data = read.csv("Bitfinex_BTCUSD_1h.csv")
