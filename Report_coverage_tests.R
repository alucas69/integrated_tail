library(GAS)

alphas_tail = c(0.1,0.05,0.025)
asset_names = c("BTC", "ETH")


alpha_tail = alphas_tail[1]
asset_name = asset_names[1]
for (asset_name in asset_names) {
  for (alpha_tail in alphas_tail) {

    core_filename = paste0("results/", asset_name, " IS ", alpha_tail/10)
    my_data = readxl::read_xlsx(paste0(core_filename, ".xlsx"))
    
    p_violation = round(100 * mean(as.integer(my_data$y > my_data$EVT_VaR)), 2)
    btest = BacktestVaR(-my_data$y, -my_data$EVT_VaR, alpha_tail/10)
    
    # print(round(c(p_violation, btest$LRuc, btest$LRcc, btest$DQ$stat, btest$DQ$pvalue), 3))
    print(round(c(p_violation, btest$LRuc, btest$LRcc), 3))
    
    
  }
}

