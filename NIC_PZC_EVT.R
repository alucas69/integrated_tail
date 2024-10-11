library(ggplot2)
library(ggpubr)
library(ggrepel)
library(latex2exp)


NIC_PZC <- function(VaR, EL, alpha, omega, A_mat, B_mat, my_data) {
  ## update using (11)+(12) followed by (9)+(16)
  NIC = rbind(
    -VaR * ((my_data$x <= VaR) - alpha),
    (my_data$x <= VaR) * my_data$x / alpha - EL
  )
  NIC = outer(omega, rep(1, nrow(my_data))) + (B_mat %*% c(VaR, EL)) %*% rep(1, nrow(my_data)) + A_mat %*% NIC
  my_data$PZC_VaR = NIC[1, ]
  my_data$PZC_EL = NIC[2, ]
  return(my_data)
}



NIC_EVT <- function(tau, f0v, f0e, tail_alpha, extreme_alpha, omega, coef_a, my_data) {
  f = omega + f0v + coef_a * (log(1 + (my_data$x/tau - 1)) - f0v) * (my_data$x <= tau)
  my_data$EVT_VaR = tau * ((extreme_alpha / tail_alpha) ^ -f);
  f = omega + f0e + coef_a * (log(1 + (my_data$x/tau - 1)) - f0e) * (my_data$x <= tau)
  my_data$EVT_EL = tau * ((extreme_alpha / tail_alpha) ^ -f) / (1 - f)
  return(my_data)
}



## set the benchmark situation
nr_grid = 1000
df_inv = 1 / 3;  tail_alpha = 0.05; extreme_alpha = 0.001
tau = qt(tail_alpha, df = 1/df_inv); 
x_min = qt(0.0001 , df = 1/df_inv); 
x_grid = tau - (tau - x_min) * (1:nr_grid)/nr_grid
y_grid = (x_grid / tau - 1)
my_data = data.frame(x = x_grid, y = y_grid)
VaR0 = qt(extreme_alpha, df = 1/df_inv)
EL0 = -(1 + df_inv * VaR0^2)/(1 - df_inv) * dt(VaR0, df = 1/df_inv) / extreme_alpha
f0v = -log(VaR0 / tau) / log(extreme_alpha / tail_alpha)
f0e = optim(f0v, function(f) { return ((tau * ((extreme_alpha / tail_alpha) ^ -f) / (1 - f) - EL0)^2)}, method = "BFGS")$par
omega_EVT = 5.318383e-05
a_EVT = 0.01462920
omega_PZC = c(-0.02901644, -0.03287773)
a_PZC = matrix(c(0.025583236, 0.0007202833, -0.007727068, 0.0008506198), nrow = 2, byrow = TRUE)
b_PZC = diag(c(0.9956943, 0.9966654));

my_data = NIC_PZC(VaR0, EL0, extreme_alpha, omega_PZC, a_PZC, b_PZC, my_data)
my_data = NIC_EVT(tau, f0v, f0e, tail_alpha, extreme_alpha, omega_EVT, a_EVT, my_data)
my_data = my_data[order(my_data$x, decreasing = TRUE), ]

plot_frame = rbind(
  data.frame(idx = (1:nr_grid), x = -my_data$x, y = my_data$y, VaRv = my_data$PZC_VaR[1]-my_data$PZC_VaR, ELv = my_data$PZC_EL[1]-my_data$PZC_EL, model = 'PZC'),
  data.frame(idx = (1:nr_grid), x = -my_data$x, y = my_data$y, VaRv = my_data$EVT_VaR[1]-my_data$EVT_VaR, ELv = my_data$EVT_EL[1]-my_data$EVT_EL, model = 'EVT')
)

idx = which(plot_frame$idx %in% (1:nr_grid))
gg1 = ggplot(data = plot_frame[idx, ], aes(x = x, y = VaRv, color = model)) +
  geom_line() + xlab('-x') + ylab('-VaR') + 
  annotate('text', x=-tau+0.1, y=-1, label= 'tau', parse = TRUE, hjust = 0) +
  annotate('text', x=-VaR0+0.1, y=-1, label= 'VaR[0.999]^{PZC}', parse = TRUE, hjust = 0) +
  geom_vline(xintercept = -VaR0, linetype = 2) +
  geom_vline(xintercept = -tau, linetype = 2)
gg2 = ggplot(data = plot_frame[idx, ], aes(x = x, y = ELv, color = model)) +
  geom_line() + xlab('-x') + ylab('-ES') +
  annotate('text', x=-tau+0.1, y=-1, label= 'tau', parse = TRUE, hjust = 0) +
  annotate('text', x=-VaR0+0.1, y=-1, label= 'VaR[PZC]', parse = TRUE, hjust = 0) +
  geom_vline(xintercept = -VaR0, linetype = 2) +
  geom_vline(xintercept = -tau, linetype = 2)
plot(ggarrange(gg1, gg2, nrow = 1, common.legend = TRUE))

# write.csv(plot_frame, file = "xin_nic.csv", row.names = FALSE)
