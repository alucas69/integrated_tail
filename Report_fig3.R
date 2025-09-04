library(ggplot2)
library(latex2exp)

asset_name = "BTC"
core_filename = paste0("results/", asset_name, " IS 0.005")
my_data = readxl::read_xlsx(paste0(core_filename, ".xlsx"))
flip_back_to_left = -1
scale_factor = 2

gg =
  ggplot(
    data = data.frame(
      x = rep(my_data$dates, 3),
      y = flip_back_to_left * c(my_data$EVT_tau, my_data$EVT_VaR, my_data$EVT_ES),
      VaR = factor(rep(c("tau", "VaR", "ES"), each = nrow(my_data)), levels = c("tau", "VaR","ES"), ordered = TRUE)
    ),
    aes(x = x, y = y, color = VaR)
  ) +
  geom_line(data = my_data, aes(x = dates, y = flip_back_to_left * y), color = 'gray') +
  geom_line() +
  coord_cartesian(ylim = c(-25,20)) +
  scale_color_manual(
    values = c("tau" = "red", "VaR" = "blue", "ES" = "aquamarine4"),
    labels = c("tau" = expression(tau[t]), "VaR", "ES"),
    name = ""
  ) + xlab("") + ylab("") +
  scale_x_datetime(
    breaks = as.POSIXct(paste0(2019:2025,"-01-01"), tz = "UTC"),
    date_labels = "%Y"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),   # remove major gridlines
    panel.grid.minor = element_blank(),   # remove minor gridlines
    # axis.title = element_text(size = 26),     # axis label size
    axis.text  = element_text(size = 24),     # tick label size
    legend.title = element_text(size = 24),   # legend title size
    legend.text  = element_text(size = 24),   # legend labels size
    legend.position = c(0.9, 0.85)           # top-left inside plot (x, y in [0,1])
  )
plot(gg)
base_size = 10; ggsave(filename = paste0(core_filename, " full.png"), plot = gg, width = base_size, height = base_size / scale_factor)


gg =
  ggplot(
    data = data.frame(
      x = rep(my_data$dates, 3),
      y = flip_back_to_left * c(my_data$EVT_tau, my_data$EVT_VaR, my_data$EVT_ES),
      VaR = factor(rep(c("tau", "VaR", "ES"), each = nrow(my_data)), levels = c("tau", "VaR","ES"), ordered = TRUE)
    ),
    aes(x = x, y = y, color = VaR)
  ) +
  geom_line(data = my_data, aes(x = dates, y = flip_back_to_left * y), color = 'gray') +
  geom_line() +
  scale_color_manual(
    values = c("tau" = "red", "VaR" = "blue", "ES" = "aquamarine4"),
    labels = c("tau" = expression(tau[t]), "VaR", "ES"),
    name = ""
  ) + xlab("") + ylab("") + 
  coord_cartesian(
    xlim = c(as.POSIXct("2021-10-01", tz = "UTC"), as.POSIXct("2023-02-28", tz = "UTC")),
    ylim = c(-9, 15)
  ) +
  geom_vline(xintercept = as.POSIXct("2022-05-10", tz = "UTC")) +
  annotate(
    "text",
    x = as.POSIXct("2022-05-10") + 20,   # shift a bit to the right of the line
    y = 12,
    label = "(1)",
    size = 10,
    hjust = 0
  ) +
  geom_vline(xintercept = as.POSIXct("2022-06-13", tz = "UTC")) +
  annotate(
    "text",
    x = as.POSIXct("2022-06-13") + 20,   # shift a bit to the right of the line
    y = 12,
    label = "(2)",
    size = 10,
    hjust = 0
  ) +
  geom_vline(xintercept = as.POSIXct("2022-11-11", tz = "UTC")) +
  annotate(
    "text",
    x = as.POSIXct("2022-11-11") + 20,   # shift a bit to the right of the line
    y = 12,
    label = "(3)",
    size = 10,
    hjust = 0
  ) +
  scale_x_datetime(
      breaks = as.POSIXct(c("2022-01-01","2022-04-01","2022-07-01","2022-10-01","2023-01-01"), tz = "UTC"),
      date_labels = c("%Y","Apr","Jul","Oct","%Y"),
    ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),   # remove major gridlines
    panel.grid.minor = element_blank(),   # remove minor gridlines
    # axis.title = element_text(size = 26),     # axis label size
    axis.text  = element_text(size = 24),     # tick label size
    legend.title = element_text(size = 24),   # legend title size
    legend.text  = element_text(size = 24),   # legend labels size
    legend.position = c(0.9, 0.85)           # top-left inside plot (x, y in [0,1])
  )
plot(gg)
base_size = 10; ggsave(filename = paste0(core_filename, " zoom.png"), plot = gg, width = base_size, height = base_size / scale_factor)



gg = ggplot(my_data, aes(x = dates)) +
  geom_ribbon(aes(ymin = EVT_bandL, ymax = EVT_bandU), fill = "lightblue", alpha = 0.8) +
  geom_line(aes(y = EVT_ft), color = "blue") + #, size = 1.2) +
  xlab("") +
  ylab("") +
  scale_x_datetime(
    breaks = as.POSIXct(paste0(2019:2025,"-01-01"), tz = "UTC"),
    date_labels = "%Y"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 24)
  )
plot(gg)
base_size = 10; ggsave(filename = paste0(core_filename, " fi.png"), plot = gg, width = base_size, height = base_size / scale_factor)



#################################################################
#################################################################
#################################################################
#################################################################


asset_name = "ETH"
core_filename = paste0("results/", asset_name, " IS 0.005")
my_data = readxl::read_xlsx(paste0(core_filename, ".xlsx"))
flip_back_to_left = -1

gg =
  ggplot(
    data = data.frame(
      x = rep(my_data$dates, 3),
      y = flip_back_to_left * c(my_data$EVT_tau, my_data$EVT_VaR, my_data$EVT_ES),
      VaR = factor(rep(c("tau", "VaR", "ES"), each = nrow(my_data)), levels = c("tau", "VaR","ES"), ordered = TRUE)
    ),
    aes(x = x, y = y, color = VaR)
  ) +
  geom_line(data = my_data, aes(x = dates, y = flip_back_to_left * y), color = 'gray') +
  geom_line() +
  coord_cartesian(ylim = c(-25,20)) +
  scale_color_manual(
    values = c("tau" = "red", "VaR" = "blue", "ES" = "aquamarine4"),
    labels = c("tau" = expression(tau[t]), "VaR", "ES"),
    name = ""
  ) + xlab("") + ylab("") +
  scale_x_datetime(
    breaks = as.POSIXct(paste0(2019:2025,"-01-01"), tz = "UTC"),
    date_labels = "%Y"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),   # remove major gridlines
    panel.grid.minor = element_blank(),   # remove minor gridlines
    # axis.title = element_text(size = 26),     # axis label size
    axis.text  = element_text(size = 24),     # tick label size
    legend.title = element_text(size = 24),   # legend title size
    legend.text  = element_text(size = 24),   # legend labels size
    legend.position = c(0.9, 0.85)           # top-left inside plot (x, y in [0,1])
  )
plot(gg)
base_size = 10; ggsave(filename = paste0(core_filename, " full.png"), plot = gg, width = base_size, height = base_size / scale_factor)


gg =
  ggplot(
    data = data.frame(
      x = rep(my_data$dates, 3),
      y = flip_back_to_left * c(my_data$EVT_tau, my_data$EVT_VaR, my_data$EVT_ES),
      VaR = factor(rep(c("tau", "VaR", "ES"), each = nrow(my_data)), levels = c("tau", "VaR","ES"), ordered = TRUE)
    ),
    aes(x = x, y = y, color = VaR)
  ) +
  geom_line(data = my_data, aes(x = dates, y = flip_back_to_left * y), color = 'gray') +
  geom_line() +
  scale_color_manual(
    values = c("tau" = "red", "VaR" = "blue", "ES" = "aquamarine4"),
    labels = c("tau" = expression(tau[t]), "VaR", "ES"),
    name = ""
  ) + xlab("") + ylab("") + 
  coord_cartesian(
    xlim = c(as.POSIXct("2021-10-01", tz = "UTC"), as.POSIXct("2023-02-28", tz = "UTC")),
    ylim = c(-9, 15)
  ) +
  geom_vline(xintercept = as.POSIXct("2022-05-10", tz = "UTC")) +
  annotate(
    "text",
    x = as.POSIXct("2022-05-10") + 20,   # shift a bit to the right of the line
    y = 12,
    label = "(1)",
    size = 10,
    hjust = 0
  ) +
  geom_vline(xintercept = as.POSIXct("2022-06-13", tz = "UTC")) +
  annotate(
    "text",
    x = as.POSIXct("2022-06-13") + 20,   # shift a bit to the right of the line
    y = 12,
    label = "(2)",
    size = 10,
    hjust = 0
  ) +
  geom_vline(xintercept = as.POSIXct("2022-11-11", tz = "UTC")) +
  annotate(
    "text",
    x = as.POSIXct("2022-11-11") + 20,   # shift a bit to the right of the line
    y = 12,
    label = "(3)",
    size = 10,
    hjust = 0
  ) +
  scale_x_datetime(
    breaks = as.POSIXct(c("2022-01-01","2022-04-01","2022-07-01","2022-10-01","2023-01-01"), tz = "UTC"),
    date_labels = c("%Y","Apr","Jul","Oct","%Y"),
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),   # remove major gridlines
    panel.grid.minor = element_blank(),   # remove minor gridlines
    # axis.title = element_text(size = 26),     # axis label size
    axis.text  = element_text(size = 24),     # tick label size
    legend.title = element_text(size = 24),   # legend title size
    legend.text  = element_text(size = 24),   # legend labels size
    legend.position = c(0.9, 0.85)           # top-left inside plot (x, y in [0,1])
  )
plot(gg)
base_size = 10; ggsave(filename = paste0(core_filename, " zoom.png"), plot = gg, width = base_size, height = base_size / scale_factor)


gg = ggplot(my_data, aes(x = dates)) +
  geom_ribbon(aes(ymin = EVT_bandL, ymax = EVT_bandU), fill = "lightblue", alpha = 0.8) +
  geom_line(aes(y = EVT_ft), color = "blue") + #, size = 1.2) +
  xlab("") +
  ylab("") +
  scale_x_datetime(
    breaks = as.POSIXct(paste0(2019:2025,"-01-01"), tz = "UTC"),
    date_labels = "%Y"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 24)
  )
plot(gg)
base_size = 10; ggsave(filename = paste0(core_filename, " fi.png"), plot = gg, width = base_size, height = base_size / scale_factor)
