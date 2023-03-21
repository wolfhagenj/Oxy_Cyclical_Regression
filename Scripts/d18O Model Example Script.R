library("data.table")
library("cmdstanr")
library("rstan")
library("ggplot2")
#
blanz_etal_data <- fread("./Data/Blanz, et al. 2020 Sheep CO Data (Cleaned).csv")
rousay_data <- blanz_etal_data[Site %in% "ROU"]
#
rousay_data[, Specimen_No := as.numeric(as.factor(Specimen))]
rousay_data[, Sample_No := 1:.N, Specimen_No]

rousay_standata <- list(
  N_Specimens = rousay_data[!is.na(d18O), .N, Specimen_No][, .N],
  N_Samples = rousay_data[!is.na(d18O), .N],
  Specimen = rousay_data[!is.na(d18O), as.numeric(as.factor(Specimen_No))],
  #observations
  dist_obs = rousay_data[!is.na(d18O), Dist],
  oxy_obs = rousay_data[!is.na(d18O), d18O],
  dist_sd = rep(0.5, rousay_data[!is.na(d18O), .N]),
  oxy_sd = rousay_data[!is.na(d18O), O.Error],
  #hyper-parameters
  prior_period = c(30, 3.5),
  period_low = 10,
  period_high = 70
)
oxy_slope_stanmodel <- cmdstan_model("./Scripts/d18O_multilevel_slopemodel.stan")
rousay_samples <- oxy_slope_stanmodel$sample(
  data = rousay_standata,
  chains = 4,
  parallel_chains = 4,
  refresh = 250,
  adapt_delta = 0.95,
  max_treedepth = 15
)

rousay_samples$summary(c("mu_period_oxy", "mu_alpha_oxy", "mu_beta_oxy", "mu_Obeta1", "mu_Obeta2", "sigma_specimen", "mu_mean_oxy", "mu_amplitude_oxy"))
rousay_samples$summary(c("period_oxy", "alpha_oxy", "beta_oxy", "Obeta1", "Obeta2"))
rousay_samples$summary("sigmasample_O")
#
rousay_samples$summary(c("sampled_period_oxy", "sampled_alpha_oxy", "sampled_beta_oxy"))
rousay_samples$summary(c("mu_max_oxy", "sampled_max_oxy"))

rousay_stanfit <- rstan::read_stan_csv(rousay_samples$output_files())

#Visualizing some of the model fit using traceplots and pair plots
traceplot(rousay_stanfit, pars = "alpha_oxy")
pairs(rousay_stanfit, pars = c("mu_alpha_oxy", "mu_beta_oxy", "mu_Obeta1", "mu_Obeta2", "sigmasample_O", "lp__"))

#extracting the posterior samples
rousay_post <- extract(rousay_stanfit)

####Functions####
#This function calculates the sinusoidal curve for a series of parameters and distance values
oxyfx <- function(params, dist) {
  period <- params[1]
  intercept <- params[2]
  slope <- params[3]
  beta1 <- params[4]
  beta2 <- params[5]
  (intercept + slope * dist) + (beta1 * cos(2 * pi * (dist / period))) + (beta2 * sin(2 * pi * (dist / period)))
}

#adds confidence intervals for observed d18O values
add_d18O_CIs <- function(data, posterior, specimen_no, col = "black") {
  #collect which observations from the posterior model to include
  observations <- which(data[, Specimen_No] %in% specimen_no)
  CI_data <- rbindlist(lapply(observations, function(x) data.table(Dist = rep(data[x, Dist], 2), CI_low = quantile(posterior[, x], 0.025), CI_high = quantile(posterior[, x], 0.975))))
  lapply(1:CI_data[, .N], function(x) geom_segment(data = CI_data[x], aes(x = Dist, xend = Dist, y = CI_low, yend = CI_high), col = col))
}

#adds isotopic profile curves (with uncertainty) and linear trend
add_isotopic_profile <- function(posterior, specimen_no, col = "black", linear_trend = T) {
  isotopic_profile_df <- rbindlist(lapply(1:nrow(posterior$alpha_oxy), function(x) data.table(Iteration = x, Dist = seq(40, 0, length.out = 1000), m_d18O = oxyfx(params = c(posterior$period_oxy[x, specimen_no], posterior$alpha_oxy[x, specimen_no], posterior$beta_oxy[x, specimen_no], posterior$Obeta1[x, specimen_no], posterior$Obeta2[x, specimen_no]), dist = seq(40, 0, length.out = 1000)))))
  linear_trend_df <- rbindlist(lapply(1:nrow(posterior$alpha_oxy), function(x) data.table(Iteration = x, Dist = seq(40, 0, length.out = 1000), m_d18O_line = posterior$alpha_oxy[x, specimen_no] + posterior$beta_oxy[x, specimen_no] * seq(40, 0, length.out = 1000))))
  if(linear_trend) {
    output_list <- list(
      #add the isotopic profile (and CI)
      geom_line(data = isotopic_profile_df[, .(m_d18O = mean(m_d18O)), .(Dist)], aes(x = Dist, y = m_d18O, linetype = "Isotopic Profile"), col = col),
      geom_ribbon(data = isotopic_profile_df[, .(CI_low = quantile(m_d18O, 0.025), CI_high = quantile(m_d18O, 0.975)), Dist], aes(x = Dist, ymin = CI_low, ymax = CI_high, linetype = "Isotopic Profile", alpha = "Isotopic Profile"), fill = col, col = col),
      #add the linear trend (and CI)
      geom_line(data = linear_trend_df[, .(m_d18O_line = mean(m_d18O_line)), Dist], aes(x = Dist, y = m_d18O_line, linetype = "Linear Trend"), col = col),
      geom_ribbon(data = linear_trend_df[, .(CI_low = quantile(m_d18O_line, 0.025), CI_high = quantile(m_d18O_line, 0.975)), Dist], aes(x = Dist, ymin = CI_low, ymax = CI_high, linetype = "Linear Trend", alpha = "Linear Trend"), fill = col, col = col)
    )
  } else {
    output_list <-   list(
      #add the isotopic profile (and CI)
      geom_line(data = isotopic_profile_df[, .(m_d18O = mean(m_d18O)), .(Dist)], aes(x = Dist, y = m_d18O, linetype = "Isotopic Profile"), col = col),
      geom_ribbon(data = isotopic_profile_df[, .(CI_low = quantile(m_d18O, 0.025), CI_high = quantile(m_d18O, 0.975)), Dist], aes(x = Dist, ymin = CI_low, ymax = CI_high, linetype = "Isotopic Profile", alpha = "Isotopic Profile"), fill = col, col = col)
    )
  }
  output_list
}

#Adds an isotopic profile using site-average parameter values
add_site_isotopic_profile <- function(posterior, col = "black") {
  isotopic_profile_df <- rbindlist(lapply(1:nrow(posterior$mu_alpha_oxy), function(x) data.table(Iteration = x, Dist = seq(40, 0, length.out = 1000), m_d18O = oxyfx(params = c(posterior$mu_period_oxy[x], posterior$mu_alpha_oxy[x], posterior$mu_beta_oxy[x], posterior$mu_Obeta1[x], posterior$mu_Obeta2[x]), dist = seq(40, 0, length.out = 1000)))))
  linear_trend_df <- rbindlist(lapply(1:nrow(posterior$mu_alpha_oxy), function(x) data.table(Iteration = x, Dist = seq(40, 0, length.out = 1000), m_d18O_line = posterior$mu_alpha_oxy[x] + posterior$mu_beta_oxy[x] * seq(40, 0, length.out = 1000))))
  list(
    #add the isotopic profile (and CI)
    geom_line(data = isotopic_profile_df[, .(m_d18O = mean(m_d18O)), .(Dist)], aes(x = Dist, y = m_d18O, linetype = "Isotopic Profile"), col = col),
    geom_polygon(data = isotopic_profile_df[, .(CI_low = quantile(m_d18O, 0.025), CI_high = quantile(m_d18O, 0.975)), Dist][, .(Dist = c(Dist[1], Dist, Dist[.N], Dist[.N], Dist[.N:1], Dist[1]), m_d18O_CI = c(CI_high[1], CI_low, CI_low[.N], CI_high[.N], CI_high[.N:1], CI_high[1]))], aes(x = Dist, y = m_d18O_CI, linetype = "Isotopic Profile", alpha = "Isotopic Profile"), fill = col, col = col),
    #add the linear trend (and CI)
    geom_line(data = linear_trend_df[, .(m_d18O_line = mean(m_d18O_line)), Dist], aes(x = Dist, y = m_d18O_line, linetype = "Linear Trend"), col = col),
    geom_polygon(data = linear_trend_df[, .(CI_low = quantile(m_d18O_line, 0.025), CI_high = quantile(m_d18O_line, 0.975)), Dist][, .(Dist = c(Dist[1], Dist, Dist[.N], Dist[.N], Dist[.N:1], Dist[1]), m_d18O_line_CI = c(CI_high[1], CI_low, CI_low[.N], CI_high[.N], CI_high[.N:1], CI_high[1]))], aes(x = Dist, y = m_d18O_line_CI, linetype = "Linear Trend", alpha = "Linear Trend"), fill = col, col = col)
  )
}

#These functions are used to create estimates of the relative position of curve peaks
#(relative to the periodicity), which is used as a relative estimate of birth seasonality.
#If values are spread out across an assemblage, then the birth seasonality was wider.

#Calculates the quantiles (default is 95%) of the relative position of the maximum value
#This is necessary because these could range across the 0 point in circular data
get_quantile_arc <- function(theta_vals, ci = 0.95, n = 1000) {
  #data are set up with theta (in radians)
  quantile_vals <- as.numeric(quantile(circular::circular(theta_vals), c((1 - ci) / 2, 1 - ((1 - ci) / 2))))
  #create the range of theta values, but first check if the values span the 0 point
  if(quantile_vals[1] <= quantile_vals[2]) {
    #doesn't span 0
    arc_vals <- seq(quantile_vals[1], quantile_vals[2], length.out = n)
  }
  if(quantile_vals[1] > quantile_vals[2]) {
    #spans 0
    arc_vals <- seq(quantile_vals[1], quantile_vals[2] + 2 * pi, length.out = n)
  }
  arc_vals
}

#Adds quantile curve of a specific quantile length to a circular plot
max_oxy_quantiles <- function(theta_vals, ci = 0.95) {
  #data are set up with theta (in radians)
  quantile_vals <- as.numeric(quantile(circular::circular(theta_vals), c((1 - ci) / 2, 1 - ((1 - ci) / 2))))
  #check if the values span the 0 point (if so, add 2*pi to the second value)
  if(quantile_vals[1] <= quantile_vals[2]) {
    #doesn't span 0
    arc_vals <- data.table(CI = ci, Min = quantile_vals[1], Max = quantile_vals[2])
  }
  if(quantile_vals[1] > quantile_vals[2]) {
    #spans 0
    arc_vals <- data.table(CI = ci, Min = quantile_vals[1], Max = quantile_vals[2] + (2 * pi))
  }
  arc_vals
}


#Adds 80% and 95% quantiles (along with median point) of the relative position of a curve's maximum in circular form
add_max_oxy_line <- function(data, radius, col = "black") {
  plot_data <- rbindlist(lapply(c(0.95, 0.80), function(x) max_oxy_quantiles(data, ci = x)))
  plot_median <- data.table(Median = as.numeric(median(circular::circular(data))))
  list(
    geom_arc(data = plot_data[CI %in% 0.95], aes(x0 = 0, y0 = 0, r = radius, start = Min, end = Max), col = col, size = 0.75),
    geom_arc(data = plot_data[CI %in% 0.80], aes(x0 = 0, y0 = 0, r = radius, start = Min, end = Max), col = col, size = 2),
    geom_point(data = plot_median, aes(x = radius * sin(Median), y = radius * cos(Median)), col = col, pch = 16, size = 4)
  )
}

#Calculates the location of all of the curves (and average values) for an assemblage to help with labeling
plot_data_creator <- function(data, posterior) {
  #collect summary information on each specimen
  plot_data <- rbindlist(lapply(1:data[, .N, Specimen_No][, .N], function(x) data.table(Specimen_No = x, Specimen = data[, .N, .(Specimen_No, Specimen)][order(Specimen_No), Specimen][x], Median = as.numeric(median(circular::circular(posterior$max_oxy[, x]))), max_oxy_quantiles(posterior$max_oxy[, x], ci = 0.95))))
  #add in the overall mean
  plot_data <- rbind(plot_data,
                     data.table(Specimen_No = (data[, .N, Specimen_No][, .N] + 1), Specimen = "Herd Average", Median = as.numeric(median(circular::circular(posterior$mu_max_oxy))), max_oxy_quantiles(posterior$mu_max_oxy, ci = 0.95)))
  
  #convert median values to x/y to organize plot, as if radius = 1 (have to kind of go by quandrant--goal is to go from closest-to-zero [vertical] to furthest-from-zero [which wraps all the way around, awkwardly, to vertical])
  plot_data[, c("Median_X", "Median_Y") := .(1 * sin(Median), 1 * cos(Median))]
  #first quadrant: X negative, Y positive: rank by decreasing X (0 = max)
  plot_data[Median_X < 0 & Median_Y >= 0, Median_Rank := frankv(plot_data[Median_X < 0 & Median_Y >= 0, Median_X], order = -1)]
  #second quadrant: X negative, Y negative: rank by decreasing Y (-1 = max)
  plot_data[Median_X < 0 & Median_Y < 0, Median_Rank := plot_data[, .N, Median_Rank][!is.na(Median_Rank), .N] + frankv(plot_data[Median_X <= 0 & Median_Y < 0, Median_X], order = -1)]
  #third quadrant: X positive, Y negative: rank by increasing X (1 = max)
  plot_data[Median_X >= 0 & Median_Y < 0, Median_Rank := plot_data[, .N, Median_Rank][!is.na(Median_Rank), .N] + frankv(plot_data[Median_X > 0 & Median_Y < 0, Median_X], order = 1)]
  #fourth quadrant: X positive, Y positive: rank by increasing Y (1 = max)
  plot_data[Median_X >= 0 & Median_Y >= 0, Median_Rank := plot_data[, .N, Median_Rank][!is.na(Median_Rank), .N] + frankv(plot_data[Median_X >= 0 & Median_Y >= 0, Median_X], order = 1)]
  
  #add in plotting radius based on Median_Rank
  plot_data[, Radius := seq(0.75, 1.25, 0.05)[Median_Rank]]
  #add "invisible" points (1000 between each Min and Max) to ensure that the labels do not overlap the lines#
  plot_data_with_points <- plot_data[, .(Specimen_No, Specimen, Radius, Min, Max)]
  for(i in 1:plot_data[!is.na(Specimen), .N, Specimen][, .N]) {
    plot_data_with_points <- rbind(
      plot_data_with_points,
      data.table(Specimen_No = NA, Specimen = "", Radius = plot_data[Specimen_No %in% i, Radius], Min = seq(plot_data[Specimen_No %in% i, Min], plot_data[Specimen_No %in% i, Max], length.out = 1000), Max = seq(plot_data[Specimen_No %in% i, Min], plot_data[Specimen_No %in% i, Max], length.out = 1000))
    )
  }
  list(Summary = plot_data, Label_Plotting = plot_data_with_points)
}

#Plotting the isotopic profile of a single specimen####
specimen_no <- 2
title_message <- paste0("Rousay Sheep M2: ", rousay_data[, .N, .(Specimen_No, Specimen)][Specimen_No %in% specimen_no, Specimen])
#
ggplot(rousay_data[Specimen_No %in% specimen_no]) +
  add_isotopic_profile(posterior = rousay_post, specimen_no = specimen_no, col = "black") +
  geom_point(aes(x = Dist, y = d18O), col = "black") +
  add_d18O_CIs(data = rousay_data, posterior = rousay_post$d18O, specimen_no = specimen_no, col = "black") +
  coord_cartesian(xlim = c(40, 0), ylim = c(-10, 0), expand = F) +
  scale_linetype_manual(name = "", values = c("Isotopic Profile" = "solid", "Linear Trend" = "dotted"), labels = c("Isotopic Profile", "Linear Trend")) +
  scale_alpha_manual(name = "", values = c("Isotopic Profile" = 0.3, "Linear Trend" = 0.1), labels = c("Isotopic Profile", "Linear Trend")) +
  labs(title = title_message, x = "Distance from REJ (mm)", y = expression(paste(delta^{18}, "O (\u2030)"))) +
  theme_classic() + theme(panel.border = element_rect(fill = NA), legend.text = element_text(size = 10), legend.position = "bottom", legend.justification = c("center"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6), plot.title = element_text(size = 16), plot.subtitle = element_text(size = 12), axis.title = element_text(size = 12), axis.text = element_text(size = 10))

#Average profile for the whole assemblage####
site_name <- "Rousay"
ggplot(rousay_data[Specimen_No %in% specimen_no, .(Dist, d18O = fourth_d18O)][!is.na(d18O)]) +
  add_site_isotopic_profile(posterior = rousay_post, col = "black") +
  coord_cartesian(xlim = c(40, 0), ylim = c(-7, -1), expand = F) +
  scale_linetype_manual(name = "", values = c("Isotopic Profile" = "solid", "Linear Trend" = "dotted"), labels = c("Isotopic Profile", "Linear Trend")) +
  scale_alpha_manual(name = "", values = c("Isotopic Profile" = 0.3, "Linear Trend" = 0.1), labels = c("Isotopic Profile", "Linear Trend")) +
  labs(title = paste(site_name, "(Average Profile)", sep = " "), x = "Distance from REJ (mm)", y = expression(paste(delta^{18}, "O (\u2030)"))) +
  theme_classic() + theme(panel.border = element_rect(fill = NA), legend.text = element_text(size = 10), legend.position = "bottom", legend.justification = c("center"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6), plot.title = element_text(size = 16), plot.subtitle = element_text(size = 12), axis.title = element_text(size = 12), axis.text = element_text(size = 10))

#Seasonality estimate for every sample with average####
full_plot_data <- plot_data_creator(data = rousay_data, posterior = rousay_post)
ggplot() +
  geom_arc(aes(x0 = 0, y0 = 0, r = 1, start = 0, end = 2 * pi), linetype = "dotted") +
  geom_segment(aes(x = -1, xend = 1, y = 0, yend = 0), linetype = "dotted") +
  geom_segment(aes(x = 0, xend = 0, y = -1, yend = 1), linetype = "dotted") +
  sapply(1:10, function(x) add_max_oxy_line(data = rousay_post$max_oxy[, x], radius = seq(0.75, 1.25, 0.05)[full_plot_data$Summary[order(Specimen_No), Median_Rank]][x], col = "grey20")) +
  add_max_oxy_line(data = rousay_post$mu_max_oxy, radius = seq(0.75, 1.25, 0.05)[full_plot_data$Summary[order(Specimen_No), Median_Rank]][11], col = "black") +
  #  geom_text_repel(data = full_plot_data$Label_Plotting[is.na(Specimen_No) | (Specimen_No %in% 1:10)], aes(x = Radius * sin(Max), y = Radius * cos(Max), label = Specimen), max.overlaps = Inf, min.segment.length = 0, size = 3) +
  geom_text_repel(data = full_plot_data$Label_Plotting[is.na(Specimen_No) | (Specimen_No %in% 11)], aes(x = Radius * sin(Min), y = Radius * cos(Min), label = Specimen), max.overlaps = Inf, col = "black", min.segment.length = 0) +
  labs(title = "Estimated Birth Seasonality: Rousay Sheep M2", x = "", y = "") +
  scale_color_manual(labels = c("Specimen", "Herd Average"), values = c("grey", "black")) +
  theme_classic() + theme(legend.text = element_text(size = 10), legend.position = "bottom", legend.justification = c("center"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6), axis.title = element_blank(), axis.text = element_blank())

#Average seasonality estimate (and expected range) for the site####
ggplot() +
  geom_arc(aes(x0 = 0, y0 = 0, r = 1, start = 0, end = 2 * pi), linetype = "dotted") +
  geom_segment(aes(x = -1, xend = 1, y = 0, yend = 0), linetype = "dotted") +
  geom_segment(aes(x = 0, xend = 0, y = -1, yend = 1), linetype = "dotted") +
  add_max_oxy_line(data = rousay_post$mu_max_oxy, radius = 0.95, col = "black") +
  add_max_oxy_line(data = rousay_post$sampled_max_oxy, radius = 1.05, col = "blue") +
  geom_arc_bar(data = max_oxy_quantiles(rousay_post$mu_max_oxy, ci = 0.95), aes(x0 = 0, y0 = 0, r0 = 0, r = 1.05, start = Min, end = Max), alpha = 0.05, fill = "black", col = "black", linetype = "dashed") +
  labs(title = "Herd Average (Black) and Range (Blue)", subtitle = site_name, x = "", y = "") +
  theme_classic() + theme(legend.text = element_text(size = 10), legend.position = "bottom", legend.justification = c("center"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6), axis.title = element_blank(), axis.text = element_blank())
