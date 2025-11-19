###################
## Set up functions
###################
library(latex2exp)
library(ggplot2)
library(ggpubr)
library(Rcpp)
library(RcppArmadillo)
library(foreach)
library(doParallel)
library(doRNG)
library(utils)
nr_free_cores = max(1, detectCores() - 1)
registerDoParallel(cores = nr_free_cores)

sourceCpp('tail_llik.cpp')



tail_llik = function(POT_data, omega, alpha, f1) {
  llik = 0
  dim_T = length(POT_data)
  for (i1 in 1:dim_T) {
    llik = llik - log(f1) - (1+1/f1)*log(1+POT_data[i1])
    f1 = omega + f1 + alpha * (log(1+POT_data[i1]) - f1)
  }
  return(llik/dim_T)
}


# set.seed(1234)

dim_T = 5000
f1_true = 0.5
omega_true = 1.5e-5
alpha_true = 0.01
dim_S = 100

POT_data = f_true = rep(0, dim_T)
if (!exists('df_sim')) df_sim = data.frame()

# compute grid 
dim_gr = 50
m_alpha = 1
m_omega = 8
alpha_grid = 0.02 * (1:(m_alpha*dim_gr) - 1)/(m_alpha*dim_gr - 1)
logomega_grid = -7.5 - 12.5 * (1:(m_omega*dim_gr) - 1)/(m_omega*dim_gr - 1)
gr_points = data.frame(
  logomega = rep(logomega_grid, each = length(alpha_grid)),
  alpha = rep(alpha_grid, length(logomega_grid)),
  llik = 0
)


for (i0 in 1:dim_S) {

  # simulate data
  cat('\nSimulation ', i0)
  f1 = f1_true
  for (i1 in 1:dim_T) {
    f_true[i1] = f1
    POT_data[i1] = exp(f1 * rexp(1)) - 1
    f1 = omega_true + f1 + alpha_true * (log(1+POT_data[i1]) - f1)
  }

  # plot(f_true, type = "l", col = "blue")

  # compute likelihood over grid and pick optimum
  f1_hat = mean(log(1+POT_data[1:50]))
  gr_points$llik = foreach(i1=1:nrow(gr_points), .combine = c) %dopar% {
    tail_llik_cpp(POT_data, exp(gr_points[i1, 'logomega']), 
              gr_points[i1, 'alpha'], f1_hat)
  }
  if (dim_S == 1) {
    # make plot and stop
    aid1 = gr_points[which.max(gr_points$llik), ]
    gr_points = subset(gr_points, alpha == aid1$alpha[1])
    gg1 = ggplot(data = gr_points, aes(y=llik, x=logomega)) + geom_line() +
      ylab('profile log-lik') + xlab('log(omega)') + ggtitle(paste0('T=',10*dim_T)) +
      geom_vline(xintercept = log(omega_true), color='red', linetype = 2) +
      geom_vline(xintercept = aid1$logomega, color='blue', linetype = 2)
    plot(gg1)
    # ggsave('simulation5000.2.png')
    stop('intermediate stop to save plot')
  }

  aid1 = gr_points[which.max(gr_points$llik),]
  df_sim[i0, paste0('alpha(T=',10*dim_T,')')] = aid1[1, 'alpha']
  df_sim[i0, paste0('logomega(T=',10*dim_T,')')] = aid1[1, 'logomega']
  df_sim[i0, paste0('llik(T=',10*dim_T,')')] = aid1[1, 'llik']
  
  if (i0 %% 10 == 0) {
    nr_bins = ceiling(dim_S^0.5)
    df_plot = data.frame(
      alpha = df_sim[,paste0('alpha(T=', 10*dim_T, ')')],
      logomega = df_sim[,paste0('logomega(T=', 10*dim_T, ')')],
      llik = df_sim[,paste0('llik(T=', 10*dim_T, ')')]
    )[1:i0, ]
    gg1 = ggplot(data = df_plot, aes(x=logomega)) + geom_histogram(aes(y=after_stat(density)), bins=nr_bins) + geom_density()
    gg2 = ggplot(data = df_plot, aes(x=alpha)) + geom_histogram(aes(y=after_stat(density)), bins=nr_bins) + geom_density()
    gg3 = ggplot(data = df_plot, aes(x=alpha, y=logomega)) + geom_jitter()
    gg = list(3); gg[[1]]=gg1; gg[[2]]=gg2; gg[[3]]=gg3; 
    plot(ggarrange(plotlist=gg,ncol=2,nrow=2))
    
    cat(paste0('\nT=', 10*dim_T, '; % at lower bound: ', 100 * length(which(df_plot[,'logomega'] < min(df_plot[,'logomega'] + 0.1)))/nrow(df_plot), '%'))
  }
}

df_plot = data.frame(
  alpha = df_sim[,paste0('alpha(T=', 10*dim_T, ')')],
  logomega = df_sim[,paste0('logomega(T=', 10*dim_T, ')')],
  llik = df_sim[,paste0('llik(T=', 10*dim_T, ')')]
)
df_plot = df_plot[complete.cases(df_plot), ]
dim_S = nrow(df_plot)
nr_bins = ceiling(dim_S^0.5)
gg1 = ggplot(data = df_plot, aes(x=logomega)) + geom_histogram(aes(y=after_stat(density)), bins=nr_bins) + geom_density()
gg2 = ggplot(data = df_plot, aes(x=alpha)) + geom_histogram(aes(y=after_stat(density)), bins=nr_bins) + geom_density()
gg3 = ggplot(data = df_plot, aes(x=alpha, y=logomega)) + geom_jitter()
gg = list(3); gg[[1]]=gg1; gg[[2]]=gg2; gg[[3]]=gg3; 
plot(ggarrange(plotlist=gg,ncol=2,nrow=2))

cat(paste0('\nT=', 10*dim_T, '; % at lower bound: ', 100 * length(which(df_plot[,'logomega'] < min(df_plot[,'logomega'] + 0.1)))/nrow(df_plot), '%'))

# save_file = 'profile_results.Rds'
# if (file.exists(save_file)) {
#   answer = menu(c('Yes','No'), title = paste0('"', save_file, '" exists: overwrite?'))
#   if (answer == 1) save(df_sim, file = save_file)
# }

# df_plot$logomega = -df_plot$logomega
# gg1 = ggplot(data = df_plot, aes(x = logomega)) + 
#   geom_histogram(aes(y = after_stat(density)), bins = nr_bins, color = 'lightgray', fill = 'lightgray') + 
#   geom_density(color = 'blue') + 
#   geom_vline(xintercept = -log(1.5e-5)) +
#   theme_bw() +   # (i) white background
#   labs(title = "T=500,000", x = TeX("$-\\log(\\omega)$"), y = "Density") +   # (iii) title and labels
#   theme(
#     axis.title.x = element_text(size = 14),   # (ii) larger x label
#     axis.title.y = element_text(size = 14),   # (ii) larger y label
#     axis.text = element_text(size=14),
#     plot.title = element_text(size = 16, face = "bold")  # (iii) larger title
#   )
# plot(gg1)
