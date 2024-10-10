library(ggplot2)
library(ggpubr)


NIC_PZC <- function(VaR, EL, alpha, omega, A_mat, B_mat, my_data) {
  ## update using (11)+(12) followed by (9)+(16)
  NIC = rbind(
    -VaR * ((my_data$y < VaR) - alpha),
    (my_data$y < VaR) * my_data$y / alpha - EL
  )
  NIC = outer(omega, rep(1, nrow(my_data))) + (B_mat %*% c(VaR, EL)) %*% rep(1, nrow(my_data)) + A_mat %*% NIC
  my_data$PZC_VaR = NIC[1, ]
  my_data$PZC_EL = NIC[1, ]
  return(my_data)
}



NIC_EVT <- function(tau, f0v, f0e, tail_alpha, extreme_alpha, omega, coef_a, my_data) {
  f = omega + f0v + coef_a * (log(1 + (my_data$y/tau - 1)) - f0v) * (my_data$y < tau)
  my_data$EVT_VaR = tau * ((extreme_alpha / tail_alpha) ^ -f);
  f = omega + f0e + coef_a * (log(1 + (my_data$y/tau - 1)) - f0e) * (my_data$y < tau)
  my_data$EVT_EL = tau * ((extreme_alpha / tail_alpha) ^ -f) / (1 - f)
  return(my_data)
}



## set the benchmark situation
nr_grid = 1000
df_inv = 1 / 3; y_grid = qt(c(0.0001, 0.4) , df = 1/df_inv); 
my_data = data.frame(x = 1:nr_grid, y = sort(y_grid[1] + (y_grid[2] - y_grid[1]) * (1:nr_grid)/nr_grid, decreasing = TRUE))
tail_alpha = 0.05; extreme_alpha = 0.001
tau = qt(tail_alpha, df = 1/df_inv); 
VaR0 = qt(extreme_alpha, df = 1/df_inv)
EL0 = -(1 + df_inv * VaR0^2)/(1 - df_inv) * dt(VaR0, df = 1/df_inv) / extreme_alpha
f0v = -log(VaR0 / tau) / log(extreme_alpha / tail_alpha)
f0e = optim(f0v, function(f) { return ((tau * ((extreme_alpha / tail_alpha) ^ -f) / (1 - f) - EL0)^2)}, method = "BFGS")$par
omega_EVT = 8.002317e-05
a_EVT = 0.01744782
omega_PZC = c(-0.02832264, -0.04838196)
a_PZC = matrix(c(-0.005810349, 0.0008353013, -0.005394566, 0.0013734608), nrow = 2, byrow = TRUE)
b_PZC = diag(c(0.9958568, 0.9955658));


my_data = NIC_PZC(VaR0, EL0, extreme_alpha, omega_PZC, a_PZC, b_PZC, my_data)
my_data = NIC_EVT(tau, f0v, f0e, tail_alpha, extreme_alpha, omega_EVT, a_EVT, my_data)

plot_frame = rbind(
  data.frame(x = my_data$x, y = -my_data$y, VaRv = my_data$PZC_VaR[1]-my_data$PZC_VaR, ELv = my_data$PZC_EL[1]-my_data$PZC_EL, model = 'PZC'),
  data.frame(x = my_data$x, y = -my_data$y, VaRv = my_data$EVT_VaR[1]-my_data$EVT_VaR, ELv = my_data$EVT_EL[1]-my_data$EVT_EL, model = 'EVT')
)

idx = which(plot_frame$x %in% (1:nr_grid))
gg1 = ggplot(data = plot_frame[idx, ], aes(x = y, y = VaRv, color = model)) +
  geom_line() + xlab('-x') + ylab('-VaR')
gg2 = ggplot(data = plot_frame[idx, ], aes(x = y, y = ELv, color = model)) +
  geom_line() + xlab('-x') + ylab('-ES')
plot(ggarrange(gg1, gg2, nrow = 1, common.legend = TRUE))

# write.csv(plot_frame, file = "xin_nic.csv", row.names = FALSE)
