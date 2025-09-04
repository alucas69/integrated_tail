library(ggplot2)

alpha_tail = c(0.1,0.05,0.025)[3]
asset_name = c("BTC", "ETH")[2]
flip_back_to_left = -1
scale_factor = 2
which_models = c("GARCH", "PZC", "EVT", "TVGPD")
which_model_names = paste0(c("GARCH", "PZC", "EVT", "TVGPD"), "_VaR")
which_models = c("EVT", "TVGPD", "PZC")
which_model_names = paste0(c("EVT", "TVGPD", "PZC"), "_VaR")

core_filename = paste0("results/", asset_name, " IS ", alpha_tail/10)
my_data = readxl::read_xlsx(paste0(core_filename, ".xlsx"))
    
idx = which( (as.Date(my_data$dates) < "2024-02-01") & (as.Date(my_data$dates) >= "2022-06-01") )
my_data = my_data[idx , ]

gg =
  ggplot(
    data = data.frame(
      x = rep(my_data$dates, length(which_models)),
      y = c(as.matrix(flip_back_to_left * my_data[ , which_model_names])),
      VaR = factor(rep(which_models, each = nrow(my_data)), levels = which_models, ordered = TRUE)
    ),
    aes(x = x, y = y, color = VaR)
  ) +
  geom_line(data = my_data, aes(x = dates, y = flip_back_to_left * y), color = 'gray') +
  geom_line() +
  coord_cartesian(ylim = c(NA, 0)) +
  scale_color_manual(
    values = c("PZC" = "blue", "EVT" = "red", "TVGPD" = "aquamarine4", "GARCH" = "darkgoldenrod3"),
    # labels = c("PZC" = expression(tau[t]), "VaR", "ES"),
    name = ""
  ) + xlab("") + ylab("") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),   # remove major gridlines
    panel.grid.minor = element_blank(),   # remove minor gridlines
    # axis.title = element_text(size = 26),     # axis label size
    axis.text  = element_text(size = 24),     # tick label size
    legend.title = element_text(size = 24),   # legend title size
    legend.text  = element_text(size = 24),   # legend labels size
    legend.position = c(0.8, 0.3)           # top-left inside plot (x, y in [0,1])
  )
plot(gg)
base_size = 10; ggsave(filename = paste0(core_filename, " pzc comp.png"), plot = gg, width = base_size, height = base_size / scale_factor)
