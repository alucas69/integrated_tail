results_folder_name = "results"
alpha_tails = c(0.1,0.05,0.025)
asset_names = c("BTC", "ETH")
nr_digits = 4
mmult = 1000

out_tab = matrix("", ncol = 3 * length(alpha_tails), nrow = 2 * length(asset_names))
colnames(out_tab) = paste0(rep(c("a", "o", "T"), length(alpha_tails)), rep(100*alpha_tails, each = 3))
rownames(out_tab) = paste0(rep(asset_names, each = 2), rep(c("", "_T"), length(asset_names)))

for (alpha_tail in alpha_tails) {
  for (asset_name in asset_names) {
    
    load(paste0(results_folder_name, "/", asset_name, " 2025-08-31 ", alpha_tail/10, ".rds"))
    out_tab[ paste0(asset_name, ""  ), paste0("a", 100 * alpha_tail)] = sprintf(paste0("%1.", nr_digits, "f"), exp(EVT_optimizer_outputs$tail.optimizer.output$par[1]))
    out_tab[ paste0(asset_name, "_T"), paste0("a", 100 * alpha_tail)] = paste0("(", sprintf(paste0("%1.", nr_digits, "f"), sqrt(EVT_optimizer_outputs$covmatrices$sandwich['alpha', 'alpha'])), ")")
    out_tab[ paste0(asset_name, ""  ), paste0("o", 100 * alpha_tail)] = sprintf(paste0("%1.", nr_digits, "f"), mmult * exp(EVT_optimizer_outputs$tail.optimizer.output$par[2]))
    out_tab[ paste0(asset_name, "_T"), paste0("o", 100 * alpha_tail)] = paste0("(", sprintf(paste0("%1.", nr_digits, "f"), mmult * sqrt(EVT_optimizer_outputs$covmatrices$sandwich['omega', 'omega'])), ")")
    out_tab[ paste0(asset_name, ""  ), paste0("T", 100 * alpha_tail)] = sprintf("%d", length(which(my_data$y > my_data$EVT_tau)))
    
  }
}


out_tab = as.data.frame(out_tab)
for (i1 in 1:nrow(out_tab)) for (i2 in 1:length(out_tab)) {
  out_tab[i1, i2] = paste0(
    ifelse(i2 == 1, " & ", ""),
    out_tab[i1, i2], 
    ifelse(i2 %% 3 == 0, " & ", ""),
    ifelse(i2 == length(out_tab), " \\", " & "))
}
print(out_tab)
