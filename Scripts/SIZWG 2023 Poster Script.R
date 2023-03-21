####SIZWG 2023 Poster Script####
#This script produces the analyses and figures used in the poster
#Examining Variability in Seasonal Behaviors with Bayesian Cyclical d18O Regression Models
#by Jesse Wolfhagen

#For strict replicability, random seeds are set for each Bayesian model,
#but these can be removed to provide conceptual replicability.

####Import libraries####
#For data analysis and manipulation
library("data.table")
library("circular")
#For Bayesian analysis
library("cmdstanr")
library("rstan")
#For visualizations
library("ggplot2")
library("ggdist")
library("ggforce")
library("ggpubr")
library("ggrepel")

####Functions####
#These functions are used to turn model parameters into isotopic profile curves
#and to create visualizations that incorporate uncertainty in the location of these curves

#This functions calculates model error (sigma_d18O) for MLE-based fits
balasse_model_fit <- function(data, par) {
  d18O_predict <- par[, amplitude_oxy] * cos(2 * pi * ((data[, Dist] - par[, max_point]) / par[, period_oxy])) + par[, mean_oxy]
  #  sum((d18O_predict - data[, d18O])^2)
  data.table(data[, .(Specimen, Dist, d18O)], d18O_predict = d18O_predict)
}

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

####Importing and setting up the data####
blanz_etal_data <- fread("./Data/Blanz, et al. 2020 Sheep CO Data (Cleaned).csv")
rousay_data <- blanz_etal_data[Site %in% "ROU"]
#
rousay_data[, Specimen_No := as.numeric(as.factor(Specimen))]
rousay_data[, Sample_No := 1:.N, Specimen_No]

#Saving every third (or fourth) sample per tooth
rousay_data[, c("third_d18O", "fourth_d18O") := .(d18O, d18O)]
rousay_data[Sample_No %% 3 != 1, third_d18O := NA]
rousay_data[Sample_No %% 4 != 1, fourth_d18O := NA]
#
#lower half of the samples (as if the tooth were more heavily worn)
rousay_data[, c("first_d18O") := .(d18O)]
rousay_data[, Max_Sample := .N, Specimen]
rousay_data[Sample_No / Max_Sample < 0.5, first_d18O := NA]

#Summary of sample counts per tooth (table in poster)
rousay_sample_table <- data.table(Sampling = c("Full", "Third", "Fourth", "First Half"), rousay_data[, .(N_Samples = .N, Third_Samples = sum(!is.na(third_d18O)), Fourth_Samples = sum(!is.na(fourth_d18O)), First_Samples = sum(!is.na(first_d18O))), .(Specimen_No, Specimen)][, .(Mean = lapply(.SD, mean), Median = lapply(.SD, median), Min = lapply(.SD, min), Max = lapply(.SD, max)), .SDcols = c("N_Samples", "Third_Samples", "Fourth_Samples", "First_Samples")])

#Calculating MLE fit (for MLE error)
published_rousay_fits <- fread("./Data/Published Rousay d18O Cyclical Regression Results.csv")
rousay_d18O_data <- rousay_data[, .(Site, Taxon = "Ovis aries", Specimen, Tooth, Dist, d18O)]
previous_d18O_fits <- published_rousay_fits[rousay_d18O_data[, .N, .(Specimen, Tooth)], on = c("Specimen", "Tooth")][, .(Site, Specimen, Tooth, period_oxy = Periodicity, amplitude_oxy = Amplitude, mean_oxy = Mean_d18O, max_point = Critical_Value )][!is.na(Site)]
#
previous_d18O_mle_fits <- rbindlist(lapply(previous_d18O_fits[, Specimen], function(x) data.table(Site = previous_d18O_fits[Specimen %in% x, Site], balasse_model_fit(rousay_d18O_data[Specimen %in% x], previous_d18O_fits[Specimen %in% x]))))
previous_d18O_mle_fits[, d18O_error := (d18O_predict - d18O)]

####Running the Bayesian analyses for the Rousay data####

#Model every sample
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
  seed = 575829583,
  max_treedepth = 15
)
rousay_stanfit <- rstan::read_stan_csv(rousay_samples$output_files())
rousay_post <- extract(rousay_stanfit)

#Model only the every-third sample, create estimates for the values of "non" sampled areas
#Evaluate true value to estimate
rousay_third_missing_standata <- list(
  N_Specimens = rousay_data[!is.na(third_d18O), .N, Specimen_No][, .N],
  N_Samples = rousay_data[!is.na(third_d18O), .N],
  N_Samples_Missing = rousay_data[is.na(third_d18O), .N],
  Specimen = rousay_data[!is.na(third_d18O), as.numeric(as.factor(Specimen_No))],
  Specimen_Missing = rousay_data[is.na(third_d18O), as.numeric(as.factor(Specimen_No))],
  #observations
  dist_obs = rousay_data[!is.na(third_d18O), Dist],
  oxy_obs = rousay_data[!is.na(third_d18O), d18O],
  dist_sd = rep(0.5, rousay_data[!is.na(third_d18O), .N]),
  oxy_sd = rousay_data[!is.na(third_d18O), O.Error],
  dist_obs_missing = rousay_data[is.na(third_d18O), Dist],
  dist_sd_missing = rep(0.5, rousay_data[is.na(third_d18O), .N]),
  #hyper-parameters
  # prior_log_period = c(3.3, 0.1)
  prior_period = c(30, 3.5),
  period_low = 10,
  period_high = 70,
  prior_sigmasample_O = c(0.08, 0.025)
)
oxy_slope_missing_stanmodel <- cmdstan_model("./Scripts/d18O_multilevel_slopemodel_missing.stan")
rousay_third_missing_samples <- oxy_slope_missing_stanmodel$sample(
  data = rousay_third_missing_standata,
  chains = 4,
  parallel_chains = 4,
  refresh = 250,
  adapt_delta = 0.95,
  seed = 911075373,
  max_treedepth = 15
)
rousay_third_missing_stanfit <- rstan::read_stan_csv(rousay_third_missing_samples$output_files())
rousay_third_missing_post <- extract(rousay_third_missing_stanfit)

#Model only the every fourth sample, create estimates for the values of "non" sampled areas
#Evaluate true value to estimate
rousay_fourth_missing_standata <- list(
  N_Specimens = rousay_data[!is.na(fourth_d18O), .N, Specimen_No][, .N],
  N_Samples = rousay_data[!is.na(fourth_d18O), .N],
  N_Samples_Missing = rousay_data[is.na(fourth_d18O), .N],
  Specimen = rousay_data[!is.na(fourth_d18O), as.numeric(as.factor(Specimen_No))],
  Specimen_Missing = rousay_data[is.na(fourth_d18O), as.numeric(as.factor(Specimen_No))],
  #observations
  dist_obs = rousay_data[!is.na(fourth_d18O), Dist],
  oxy_obs = rousay_data[!is.na(fourth_d18O), d18O],
  dist_sd = rep(0.5, rousay_data[!is.na(fourth_d18O), .N]),
  oxy_sd = rousay_data[!is.na(fourth_d18O), O.Error],
  dist_obs_missing = rousay_data[is.na(fourth_d18O), Dist],
  dist_sd_missing = rep(0.5, rousay_data[is.na(fourth_d18O), .N]),
  #hyper-parameters
  prior_period = c(30, 3.5),
  period_low = 10,
  period_high = 70,
  prior_sigmasample_O = c(0.08, 0.025)
)
oxy_slope_missing_stanmodel <- cmdstan_model("./Scripts/d18O_multilevel_slopemodel_missing.stan")
rousay_fourth_missing_samples <- oxy_slope_missing_stanmodel$sample(
  data = rousay_fourth_missing_standata,
  chains = 4,
  parallel_chains = 4,
  refresh = 250,
  adapt_delta = 0.95,
  seed = 1386218260,
  max_treedepth = 15
)
rousay_fourth_missing_stanfit <- rstan::read_stan_csv(rousay_fourth_missing_samples$output_files())
rousay_fourth_missing_post <- extract(rousay_fourth_missing_stanfit)

#Model only the lower half of the samples, create estimates for the values of "non" sampled areas
#Evaluate true value to estimate
rousay_first_missing_standata <- list(
  N_Specimens = rousay_data[!is.na(first_d18O), .N, Specimen_No][, .N],
  N_Samples = rousay_data[!is.na(first_d18O), .N],
  N_Samples_Missing = rousay_data[is.na(first_d18O), .N],
  Specimen = rousay_data[!is.na(first_d18O), as.numeric(as.factor(Specimen_No))],
  Specimen_Missing = rousay_data[is.na(first_d18O), as.numeric(as.factor(Specimen_No))],
  #observations
  dist_obs = rousay_data[!is.na(first_d18O), Dist],
  oxy_obs = rousay_data[!is.na(first_d18O), d18O],
  dist_sd = rep(0.5, rousay_data[!is.na(first_d18O), .N]),
  oxy_sd = rousay_data[!is.na(first_d18O), O.Error],
  dist_obs_missing = rousay_data[is.na(first_d18O), Dist],
  dist_sd_missing = rep(0.5, rousay_data[is.na(first_d18O), .N]),
  #hyper-parameters
  # prior_log_period = c(3.3, 0.1)
  prior_period = c(30, 3.5),
  period_low = 10,
  period_high = 70,
  prior_sigmasample_O = c(0.08, 0.025)
)
oxy_slope_missing_stanmodel <- cmdstan_model("./Scripts/d18O_multilevel_slopemodel_missing.stan")
rousay_first_missing_samples <- oxy_slope_missing_stanmodel$sample(
  data = rousay_first_missing_standata,
  chains = 4,
  parallel_chains = 4,
  refresh = 250,
  adapt_delta = 0.99,
  seed = 678709198,
  max_treedepth = 15
)
rousay_first_missing_stanfit <- rstan::read_stan_csv(rousay_first_missing_samples$output_files())
rousay_first_missing_post <- extract(rousay_first_missing_stanfit)

####Creating visualizations####

#make plots to get the legend for the different sampling intensities
sampling_intensity_legend_dt <- data.table(
  Iteration = rep(rep(1:4000), 4),
  `Sampling Intensity` = factor(rep(c("Every Sample", "Every Third", "Every Fourth", "Lower Half"), each = 4000), levels = c("Every Sample", "Every Third", "Every Fourth", "Lower Half")),
  Value = c(rousay_post$mu_alpha_oxy, rousay_third_missing_post$mu_alpha_oxy, rousay_fourth_missing_post$mu_alpha_oxy, rousay_first_missing_post$mu_alpha_oxy)
)
legend_plot <-
  ggplot(sampling_intensity_legend_dt) + aes(x = Value) +
  geom_histogram(aes(fill = `Sampling Intensity`)) +
  scale_fill_manual(name = "Sampling Intensity", values = c("black", "blue", "red", "orange"), na.translate = F) +
  guides(fill = guide_legend(ncol = 2, nrow = 2, title.position = "top", title.hjust = 0.5)) +
  theme_classic() + theme(legend.text = element_text(size = 10), legend.position = "bottom", legend.justification = c("center"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6), axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5), axis.title = element_text(size = 16), axis.text.y = element_text(angle = 0, size = 10, hjust = 0.5))
one_row_legend_plot <-
  ggplot(sampling_intensity_legend_dt) + aes(x = Value) +
  geom_histogram(aes(fill = `Sampling Intensity`)) +
  scale_fill_manual(name = "Sampling Intensity", values = c("black", "blue", "red", "orange"), na.translate = F) +
  guides(fill = guide_legend(nrow = 1, title.position = "top", title.hjust = 0.5)) +
  theme_classic() + theme(legend.text = element_text(size = 10), legend.position = "bottom", legend.justification = c("center"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6), axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5), axis.title = element_text(size = 16), axis.text.y = element_text(angle = 0, size = 10, hjust = 0.5))

#Example isotopic profile (with parameter annotations)####
regression_example_dt <- data.table(
  Dist = seq(50, -10, length.out = 1000),
  mu_d18O = oxyfx(params = c(30, -3.5, -0.03, -0.6, 1.2), dist = seq(50, -10, length.out = 1000))
)
#calculating key locations on the curve
example_oxy_max <- (((atan(1.2 / -0.6) * (30 / (2 * pi))) + (30 / 2)))
example_oxy_min <- (((atan(1.2 / -0.6) * (30 / (2 * pi))) + (30)))
#
example_regression_plot <- 
  ggplot(regression_example_dt) + aes(x = Dist, y = mu_d18O) +
  #the curve
  geom_line(lwd = 1.2) +
  #mean line (beta)
  geom_line(data = data.table(Dist = seq(50, -10, length.out = 1000), mu_d18O = -3.5 - 0.03 * seq(50, -10, length.out = 1000)), linetype = "solid") +
  geom_segment(data = data.table(Dist = c(40.5, 37), mu_d18O = rep((-3.5 - 0.03 * 40), 2)), aes(x = Dist[1], xend = Dist[2], y = mu_d18O[1], yend = mu_d18O[2]), linetype = "solid", col = "black", alpha = 0.2) +
  geom_segment(data = data.table(Dist = c(37, 37), mu_d18O = c(-3.5 - 0.03 * 40, -3.5 - 0.03 * 37)), aes(x = Dist[1], xend = Dist[2], y = mu_d18O[1], yend = mu_d18O[2]), linetype = "solid", col = "black", alpha = 0.2) +
  annotate("text", x = 37, y = (-3.5 - 0.03 * 37) - 0.15, label = expression(paste("Slope (", beta[0], ")"))) +
  #alpha
  geom_point(data = data.table(Dist = 0, mu_d18O = -3.5), col = "blue", size = 3) +
  geom_segment(data = data.table(Dist = c(0, 0), mu_d18O = c(-3.5, -10)), aes(x = Dist[1], xend = Dist[2], y = mu_d18O[1], yend = mu_d18O[2]), col = "blue", linetype = "dotdash", alpha = 0.2) +
  annotate("text", x = 0 + 0.45, y = (-3.5) + 0.15, label = expression(paste("Intercept (", alpha, ")")), col = "blue") +
  #periodicity
  #  geom_segment(data = data.table(Dist = c(example_oxy_max, example_oxy_max + 30), mu_d18O = oxyfx(params = c(30, -3.5, -0.03, -0.6, 1.2), dist = c(example_oxy_max, example_oxy_max + 30))), aes(x = Dist[1], xend = Dist[2], y = mu_d18O[1], yend = mu_d18O[2]), linetype = "dotdash", arrow = arrow(angle = 90, length = unit(0.1, "inches"), ends = "both"), col = "red") +
  geom_segment(data = data.table(Dist = c(example_oxy_max, example_oxy_max + 30), mu_d18O = ((oxyfx(params = c(30, -3.5, -0.03, -0.6, 1.2), dist = c(example_oxy_max, example_oxy_max + 0))))), aes(x = Dist[1], xend = Dist[2], y = mu_d18O[1], yend = mu_d18O[2]), linetype = "dashed", arrow = arrow(angle = 90, length = unit(0.1, "inches"), ends = "both"), col = "red") +
  annotate("text", x = example_oxy_max + 15, y = ((oxyfx(params = c(30, -3.5, -0.03, -0.6, 1.2), dist = c(example_oxy_max))) - 0.15), label = expression(paste("Periodicity (", lambda, ")")), col = "red") +
  #amplitude
  geom_segment(data = data.table(Dist = c(example_oxy_min, example_oxy_min), mu_d18O = c(oxyfx(params = c(30, -3.5, -0.03, -0.6, 1.2), dist = c(example_oxy_min)), -3.5 - 0.03 * example_oxy_min)), aes(x = Dist[1], xend = Dist[2], y = mu_d18O[1], yend = mu_d18O[2]), linetype = "dotted", arrow = arrow(angle = 90, length = unit(0.1, "inches"), ends = "both"), col = "orange", lwd = 0.8) +
  geom_segment(data = data.table(Dist = c(example_oxy_max, example_oxy_max), mu_d18O = c(oxyfx(params = c(30, -3.5, -0.03, -0.6, 1.2), dist = c(example_oxy_max)), -3.5 - 0.03 * example_oxy_max)), aes(x = Dist[1], xend = Dist[2], y = mu_d18O[1], yend = mu_d18O[2]), linetype = "dotted", arrow = arrow(angle = 90, length = unit(0.1, "inches"), ends = "both"), col = "orange", lwd = 0.8) +
  annotate("text", x = (example_oxy_min + 3), y = (mean(c(oxyfx(params = c(30, -3.5, -0.03, -0.6, 1.2), dist = c(example_oxy_min)), -3.5 - 0.03 * example_oxy_min)) + 0.15), label = expression(paste("Amplitude (", Alpha, ")")), col = "orange") +
  annotate("text", x = (example_oxy_max - 3), y = (mean(c(oxyfx(params = c(30, -3.5, -0.03, -0.6, 1.2), dist = c(example_oxy_max)), -3.5 - 0.03 * example_oxy_max)) - 0.15), label = expression(paste("Amplitude (", Alpha, ")")), col = "orange") +
  #the equations
  annotate("text", x = 20, y = -2, label = "delta^{18}~O~{}(x) == alpha + (beta[0] * x) + beta[1] * cos(2 * pi * frac(x, lambda)) + beta[2] * sin(2 * pi * frac(x, lambda))", parse = T) +
  annotate("text", x = 20, y = -2.25, label = "Alpha == sqrt(beta[1]^{2} + beta[2]^{2})", parse = T) +
  coord_cartesian(xlim = c(40, 0), ylim = c(-6, -2), expand = T) +
  labs(x = "Distance from REJ (x) (mm)", y = expression(paste(delta^{18}, "O (\u2030)"))) +
  theme_classic() + theme(panel.border = element_rect(fill = NA), legend.text = element_text(size = 10), legend.position = "bottom", legend.justification = c("center"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6), plot.title = element_text(size = 16), plot.subtitle = element_text(size = 12), axis.title = element_text(size = 12), axis.text = element_text(size = 10))
example_regression_plot

#Tooth-specific parameter comparison####
tooth_alpha_oxy_df <- data.table(
  Iteration = rep(rep(1:4000, 10), 4),
  `Sampling Intensity` = factor(rep(c("Every Sample", "Every Third", "Every Fourth", "Lower Half"), each = 4000 * rousay_standata$N_Specimens), levels = c("Every Sample", "Every Third", "Every Fourth", "Lower Half")),
  Specimen = rep(factor(rep(rousay_data[, .N, .(Specimen_No, Specimen)][order(Specimen_No), Specimen], each = 4000), levels = rousay_data[, .N, .(Specimen_No, Specimen)][order(Specimen_No), Specimen]), 4),
  Value = c(rousay_post$alpha_oxy, rousay_third_missing_post$alpha_oxy, rousay_fourth_missing_post$alpha_oxy, rousay_first_missing_post$alpha_oxy)
)
tooth_alpha_plot <-
  ggplot(tooth_alpha_oxy_df[, .(`Sampling Intensity`, Specimen, Value, Plot_Group = paste(Specimen, as.numeric(`Sampling Intensity`), sep = "."))]) + aes(y = Value, x = Plot_Group) +
  stat_slab(normalize = "groups", aes(fill = `Sampling Intensity`, fill_ramp = stat(cut_cdf_qi(cdf, .width = c(1, 0.95, 0.80), labels = scales::percent_format())))) +
  stat_pointinterval(.width = c(0.8, 0.95), aes(color = `Sampling Intensity`), position = position_dodge(width = 0.2, preserve = "single")) +
  scale_fill_ramp_discrete(name = "Interval", range = c(0.8, 0.3), na.translate = F) +
  scale_color_manual(name = "Sampling Intensity", values = c("black", "blue", "red", "orange"), na.translate = F) +
  scale_fill_manual(name = "Sampling Intensity", values = c("black", "blue", "red", "orange"), na.translate = F) +
  scale_x_discrete(name = "", breaks = tooth_period_oxy_df[, .N, .(Specimen)][order(Specimen), paste(Specimen, ".2", sep = "")], labels = tooth_period_oxy_df[, .N, .(Specimen)][order(Specimen), Specimen]) +
  guides(color = guide_legend(nrow = 2, ncol = 2, title.position = "top", title.hjust = 0.5), fill = guide_legend(nrow = 2, ncol = 2, title.position = "top", title.hjust = 0.5)) +
  labs(y = expression(paste("Intercept: ", alpha))) +
  theme_classic() + theme(legend.text = element_text(size = 10), legend.position = "bottom", legend.justification = c("center"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6), axis.ticks.x = element_blank(), axis.text.x = element_text(size = 10, angle = -90, vjust = 0), axis.title = element_text(size = 16), axis.text.y = element_text(angle = 0, size = 10, hjust = 0.5))
#
tooth_beta_oxy_df <- data.table(
  Iteration = rep(rep(1:4000, 10), 4),
  `Sampling Intensity` = factor(rep(c("Every Sample", "Every Third", "Every Fourth", "Lower Half"), each = 4000 * rousay_standata$N_Specimens), levels = c("Every Sample", "Every Third", "Every Fourth", "Lower Half")),
  Specimen = rep(factor(rep(rousay_data[, .N, .(Specimen_No, Specimen)][order(Specimen_No), Specimen], each = 4000), levels = rousay_data[, .N, .(Specimen_No, Specimen)][order(Specimen_No), Specimen]), 4),
  Value = c(rousay_post$beta_oxy, rousay_third_missing_post$beta_oxy, rousay_fourth_missing_post$beta_oxy, rousay_first_missing_post$beta_oxy)
)
tooth_beta_plot <-
  ggplot(tooth_beta_oxy_df[, .(`Sampling Intensity`, Specimen, Value, Plot_Group = paste(Specimen, as.numeric(`Sampling Intensity`), sep = "."))]) + aes(y = Value, x = Plot_Group) +
  stat_slab(normalize = "groups", aes(fill = `Sampling Intensity`, fill_ramp = stat(cut_cdf_qi(cdf, .width = c(1, 0.95, 0.80), labels = scales::percent_format())))) +
  stat_pointinterval(.width = c(0.8, 0.95), aes(color = `Sampling Intensity`), position = position_dodge(width = 0.2, preserve = "single")) +
  scale_fill_ramp_discrete(name = "Interval", range = c(0.8, 0.3), na.translate = F) +
  scale_color_manual(name = "Sampling Intensity", values = c("black", "blue", "red", "orange"), na.translate = F) +
  scale_fill_manual(name = "Sampling Intensity", values = c("black", "blue", "red", "orange"), na.translate = F) +
  scale_x_discrete(name = "", breaks = tooth_period_oxy_df[, .N, .(Specimen)][order(Specimen), paste(Specimen, ".2", sep = "")], labels = tooth_period_oxy_df[, .N, .(Specimen)][order(Specimen), Specimen]) +
  guides(color = guide_legend(nrow = 2, ncol = 2, title.position = "top", title.hjust = 0.5), fill = guide_legend(nrow = 2, ncol = 2, title.position = "top", title.hjust = 0.5)) +
  labs(y = expression(paste("Slope: ", beta))) +
  theme_classic() + theme(legend.text = element_text(size = 10), legend.position = "bottom", legend.justification = c("center"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6), axis.ticks.x = element_blank(), axis.text.x = element_text(size = 10, angle = -90, vjust = 0), axis.title = element_text(size = 16), axis.text.y = element_text(angle = 0, size = 10, hjust = 0.5))
#
tooth_amplitude_oxy_df <- data.table(
  Iteration = rep(rep(1:4000, 10), 4),
  `Sampling Intensity` = factor(rep(c("Every Sample", "Every Third", "Every Fourth", "Lower Half"), each = 4000 * rousay_standata$N_Specimens), levels = c("Every Sample", "Every Third", "Every Fourth", "Lower Half")),
  Specimen = rep(factor(rep(rousay_data[, .N, .(Specimen_No, Specimen)][order(Specimen_No), Specimen], each = 4000), levels = rousay_data[, .N, .(Specimen_No, Specimen)][order(Specimen_No), Specimen]), 4),
  Value = c(rousay_post$amplitude_oxy, rousay_third_missing_post$amplitude_oxy, rousay_fourth_missing_post$amplitude_oxy, rousay_first_missing_post$amplitude_oxy)
)
tooth_amplitude_plot <-
  ggplot(tooth_amplitude_oxy_df[, .(`Sampling Intensity`, Specimen, Value, Plot_Group = paste(Specimen, as.numeric(`Sampling Intensity`), sep = "."))]) + aes(y = Value, x = Plot_Group) +
  stat_slab(normalize = "groups", aes(fill = `Sampling Intensity`, fill_ramp = stat(cut_cdf_qi(cdf, .width = c(1, 0.95, 0.80), labels = scales::percent_format())))) +
  stat_pointinterval(.width = c(0.8, 0.95), aes(color = `Sampling Intensity`), position = position_dodge(width = 0.2, preserve = "single")) +
  scale_fill_ramp_discrete(name = "Interval", range = c(0.8, 0.3), na.translate = F) +
  scale_color_manual(name = "Sampling Intensity", values = c("black", "blue", "red", "orange"), na.translate = F) +
  scale_fill_manual(name = "Sampling Intensity", values = c("black", "blue", "red", "orange"), na.translate = F) +
  scale_x_discrete(name = "", breaks = tooth_period_oxy_df[, .N, .(Specimen)][order(Specimen), paste(Specimen, ".2", sep = "")], labels = tooth_period_oxy_df[, .N, .(Specimen)][order(Specimen), Specimen]) +
  guides(color = guide_legend(nrow = 2, ncol = 2, title.position = "top", title.hjust = 0.5), fill = guide_legend(nrow = 2, ncol = 2, title.position = "top", title.hjust = 0.5)) +
  labs(y = expression(paste("Amplitude: ", Alpha))) +
  theme_classic() + theme(legend.text = element_text(size = 10), legend.position = "bottom", legend.justification = c("center"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6), axis.ticks.x = element_blank(), axis.text.x = element_text(size = 10, angle = -90, vjust = 0), axis.title = element_text(size = 16), axis.text.y = element_text(angle = 0, size = 10, hjust = 0.5))
#
tooth_period_oxy_df <- data.table(
  Iteration = rep(rep(1:4000, 10), 4),
  `Sampling Intensity` = factor(rep(c("Every Sample", "Every Third", "Every Fourth", "Lower Half"), each = 4000 * rousay_standata$N_Specimens), levels = c("Every Sample", "Every Third", "Every Fourth", "Lower Half")),
  Specimen = rep(factor(rep(rousay_data[, .N, .(Specimen_No, Specimen)][order(Specimen_No), Specimen], each = 4000), levels = rousay_data[, .N, .(Specimen_No, Specimen)][order(Specimen_No), Specimen]), 4),
  Value = c(rousay_post$period_oxy, rousay_third_missing_post$period_oxy, rousay_fourth_missing_post$period_oxy, rousay_first_missing_post$period_oxy)
)
tooth_period_plot <-
  ggplot(tooth_period_oxy_df[, .(`Sampling Intensity`, Specimen, Value, Plot_Group = paste(Specimen, as.numeric(`Sampling Intensity`), sep = "."))]) + aes(y = Value, x = Plot_Group) +
  stat_slab(normalize = "groups", aes(fill = `Sampling Intensity`, fill_ramp = stat(cut_cdf_qi(cdf, .width = c(1, 0.95, 0.80), labels = scales::percent_format())))) +
  stat_pointinterval(.width = c(0.8, 0.95), aes(color = `Sampling Intensity`), position = position_dodge(width = 0.2, preserve = "single")) +
  scale_fill_ramp_discrete(name = "Interval", range = c(0.8, 0.3), na.translate = F) +
  scale_color_manual(name = "Sampling Intensity", values = c("black", "blue", "red", "orange"), na.translate = F) +
  scale_fill_manual(name = "Sampling Intensity", values = c("black", "blue", "red", "orange"), na.translate = F) +
  scale_x_discrete(name = "", breaks = tooth_period_oxy_df[, .N, .(Specimen)][order(Specimen), paste(Specimen, ".2", sep = "")], labels = tooth_period_oxy_df[, .N, .(Specimen)][order(Specimen), Specimen]) +
  labs(y = expression(paste("Periodicity: ", lambda))) +
  guides(color = guide_legend(nrow = 2, ncol = 2, title.position = "top", title.hjust = 0.5), fill = guide_legend(nrow = 2, ncol = 2, title.position = "top", title.hjust = 0.5)) +
  theme_classic() + theme(legend.text = element_text(size = 10), legend.position = "bottom", legend.justification = c("center"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6), axis.ticks.x = element_blank(), axis.text.x = element_text(size = 10, angle = -90, vjust = 0), axis.title = element_text(size = 16), axis.text.y = element_text(angle = 0, size = 10, hjust = 0.5))
#
tooth_parameters <- ggarrange(tooth_alpha_plot, tooth_beta_plot, tooth_period_plot, tooth_amplitude_plot, nrow = 2, ncol = 2, common.legend = T, legend = "bottom")
tooth_parameter_plot <- annotate_figure(tooth_parameters, top = "Specimen-Specific Parameter Estimates")
tooth_parameter_plot

#Error plot comparison####
full_d18O_obs_error <- rbindlist(lapply(1:rousay_data[, .N], function(x) data.table(`Sampling Intensity` = "Every Sample", Type = "Observation", rousay_data[!is.na(d18O)][x, .(Specimen_No, Sample_No)], d18O_Error = rousay_data[!is.na(d18O)][x, d18O] - rousay_post$d18O[, x])))
third_d18O_error <- rbindlist(lapply(1:rousay_data[is.na(third_d18O), .N], function(x) data.table(`Sampling Intensity` = "Every Third", Type = "Estimation", rousay_data[is.na(third_d18O)][x, .(Specimen_No, Sample_No)], d18O_Error = rousay_data[is.na(third_d18O)][x, d18O] - rousay_third_missing_post$d18O_missing[, x])))
third_d18O_obs_error <- rbindlist(lapply(1:rousay_data[!is.na(third_d18O), .N], function(x) data.table(`Sampling Intensity` = "Every Third", Type = "Observation", rousay_data[!is.na(third_d18O)][x, .(Specimen_No, Sample_No)], d18O_Error = rousay_data[!is.na(third_d18O)][x, d18O] - rousay_third_missing_post$d18O[, x])))
fourth_d18O_error <- rbindlist(lapply(1:rousay_data[is.na(fourth_d18O), .N], function(x) data.table(`Sampling Intensity` = "Every Fourth", Type = "Estimation", rousay_data[is.na(fourth_d18O)][x, .(Specimen_No, Sample_No)], d18O_Error = rousay_data[is.na(fourth_d18O)][x, d18O] - rousay_fourth_missing_post$d18O_missing[, x])))
fourth_d18O_obs_error <- rbindlist(lapply(1:rousay_data[!is.na(fourth_d18O), .N], function(x) data.table(`Sampling Intensity` = "Every Fourth", Type = "Observation", rousay_data[!is.na(fourth_d18O)][x, .(Specimen_No, Sample_No)], d18O_Error = rousay_data[!is.na(fourth_d18O)][x, d18O] - rousay_fourth_missing_post$d18O[, x])))
first_d18O_error <- rbindlist(lapply(1:rousay_data[is.na(first_d18O), .N], function(x) data.table(`Sampling Intensity` = "Lower Half", Type = "Estimation", rousay_data[is.na(first_d18O)][x, .(Specimen_No, Sample_No)], d18O_Error = rousay_data[is.na(first_d18O)][x, d18O] - rousay_first_missing_post$d18O_missing[, x])))
first_d18O_obs_error <- rbindlist(lapply(1:rousay_data[!is.na(first_d18O), .N], function(x) data.table(`Sampling Intensity` = "Lower Half", Type = "Observation", rousay_data[!is.na(first_d18O)][x, .(Specimen_No, Sample_No)], d18O_Error = rousay_data[!is.na(first_d18O)][x, d18O] - rousay_first_missing_post$d18O[, x])))

d18O_error_data <- rbind(full_d18O_obs_error, third_d18O_error, third_d18O_obs_error, fourth_d18O_error, fourth_d18O_obs_error, first_d18O_error, first_d18O_obs_error)
d18O_error_data <- rbind(d18O_error_data, previous_d18O_mle_fits[Site %in% "Rousay", .(`Sampling Intensity` = "MLE Fit", Type = "Observation", Specimen_No = Specimen, Sample_No = NA, d18O_Error = d18O_error)])
d18O_error_data[, `Sampling Intensity` := factor(`Sampling Intensity`, levels = c("Every Sample", "Every Third", "Every Fourth", "Lower Half", "MLE Fit"))]
d18O_error_data[`Sampling Intensity` %in% "Lower Half", Max_Observed_Sample := min(Sample_No[Type %in% "Observation"]), Specimen_No]
d18O_error_data[`Sampling Intensity` %in% "Lower Half" & Type %in% "Estimation", Estimated_Group := ifelse(Max_Observed_Sample - Sample_No <= 5, "Nearest 5 Samples", ifelse(Sample_No %in% 1:5, "Furthest 5 Samples", "Intermediate Samples"))]
d18O_error_data[, Estimated_Group := factor(Estimated_Group, levels = c("Nearest 5 Samples", "Intermediate Samples", "Furthest 5 Samples"))]
#observed data points
observed_error_plot <- 
  ggplot(data = d18O_error_data[Type %in% "Observation", .(`Sampling Intensity`, d18O_Error = abs(d18O_Error))]) + aes(y = `Sampling Intensity`, x = d18O_Error) +
  stat_slab(normalize = "groups", slab_type = "ccdf", aes(fill = `Sampling Intensity`, fill_ramp = stat(cut_cdf_qi(cdf, .width = c(1, 0.95, 0.80), labels = scales::percent_format())))) +
  stat_pointinterval(.width = c(0.8, 0.95), slab_type = "ccdf", aes(color = `Sampling Intensity`), position = position_dodge(width = 0.2, preserve = "single")) +
  geom_vline(data = d18O_error_data[`Sampling Intensity` %in% "Every Sample", .(d18O_Error = quantile(abs(d18O_Error), 0.975))], aes(xintercept = d18O_Error), linetype = "dashed") +
  scale_fill_ramp_discrete(name = "Interval", range = c(0.8, 0.3), na.translate = F) +
  scale_color_manual(name = "Sampling Intensity", values = c("Every Sample" = "black", "Every Third" = "blue", "Every Fourth" = "red", "Lower Half" = "orange"), na.translate = T) +
  scale_fill_manual(name = "Sampling Intensity", values = c("Every Sample" = "black", "Every Third" = "blue", "Every Fourth" = "red", "Lower Half" = "orange"), na.translate = T) +
  scale_y_discrete(name = "") +
  scale_x_continuous(breaks = seq(0, 10, 0.5)) + coord_cartesian(expand = FALSE) +
  labs(subtitle = "Error on Observed Data Points", x = expression(paste(delta^{18}, "O (\u2030) Error"))) +
  guides(color = guide_legend(nrow = 2, ncol = 2, title.position = "top", title.hjust = 0.5), fill = guide_legend(nrow = 2, ncol = 2, title.position = "top", title.hjust = 0.5)) +
  theme_classic() + theme(legend.text = element_text(size = 10), legend.position = "bottom", legend.justification = c("center"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6), axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5), axis.title = element_text(size = 16), axis.text.y = element_text(angle = 0, size = 10, hjust = 1))
#omitted data points
estimated_error_plot <-
  ggplot(data = d18O_error_data[Type %in% "Estimation", .(`Sampling Intensity`, d18O_Error = abs(d18O_Error))]) + aes(y = `Sampling Intensity`, x = d18O_Error) +
  stat_slab(normalize = "groups", slab_type = "ccdf", aes(fill = `Sampling Intensity`, fill_ramp = stat(cut_cdf_qi(cdf, .width = c(1, 0.95, 0.80), labels = scales::percent_format())))) +
  stat_pointinterval(.width = c(0.8, 0.95), slab_type = "ccdf", aes(color = `Sampling Intensity`), position = position_dodge(width = 0.2, preserve = "single")) +
  geom_vline(data = d18O_error_data[`Sampling Intensity` %in% "Every Sample", .(d18O_Error = quantile(abs(d18O_Error), 0.95))], aes(xintercept = d18O_Error), linetype = "dashed") +
  scale_fill_ramp_discrete(name = "Interval", range = c(0.8, 0.3), na.translate = F) +
  scale_color_manual(name = "Sampling Intensity", values = c("Every Sample" = "black", "Every Third" = "blue", "Every Fourth" = "red", "Lower Half" = "orange"), na.translate = T) +
  scale_fill_manual(name = "Sampling Intensity", values = c("Every Sample" = "black", "Every Third" = "blue", "Every Fourth" = "red", "Lower Half" = "orange"), na.translate = T) +
  scale_y_discrete(name = "") +
  scale_x_continuous(breaks = seq(0, 10, 1)) + coord_cartesian(expand = FALSE) +
  labs(subtitle = "Error on Omitted Data Points", x = expression(paste(delta^{18}, "O (\u2030) Error"))) +
  guides(color = guide_legend(nrow = 2, ncol = 2, title.position = "top", title.hjust = 0.5), fill = guide_legend(nrow = 2, ncol = 2, title.position = "top", title.hjust = 0.5)) +
  theme_classic() + theme(legend.text = element_text(size = 10), legend.position = "bottom", legend.justification = c("center"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6), axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5), axis.title = element_text(size = 16), axis.text.y = element_text(angle = 0, size = 10, hjust = 1))
#lower half (by distance from observation)
lowerhalf_error_plot <-
  ggplot(data = d18O_error_data[!is.na(Estimated_Group), .(`Sampling Intensity`, d18O_Error = abs(d18O_Error), Estimated_Group)]) + aes(y = `Estimated_Group`, x = d18O_Error) +
  stat_slab(normalize = "groups", slab_type = "ccdf", aes(fill = `Sampling Intensity`, fill_ramp = stat(cut_cdf_qi(cdf, .width = c(1, 0.95, 0.80), labels = scales::percent_format())))) +
  stat_pointinterval(.width = c(0.8, 0.95), slab_type = "ccdf", aes(color = `Sampling Intensity`), position = position_dodge(width = 0.2, preserve = "single")) +
  geom_vline(data = d18O_error_data[`Sampling Intensity` %in% "Every Sample", .(d18O_Error = quantile(abs(d18O_Error), 0.95))], aes(xintercept = d18O_Error), linetype = "dashed") +
  scale_fill_ramp_discrete(name = "Interval", range = c(0.8, 0.3), na.translate = F) +
  scale_color_manual(name = "Sampling Intensity", values = c("Every Sample" = "black", "Every Third" = "blue", "Every Fourth" = "red", "Lower Half" = "orange"), na.translate = T) +
  scale_fill_manual(name = "Sampling Intensity", values = c("Every Sample" = "black", "Every Third" = "blue", "Every Fourth" = "red", "Lower Half" = "orange"), na.translate = T) +
  scale_y_discrete(name = "") +
  scale_x_continuous(breaks = seq(0, 10, 1)) + coord_cartesian(expand = FALSE) +
  labs(subtitle = "Error by Distance from Observations\n('Lower Half' Strategy)", x = expression(paste(delta^{18}, "O (\u2030) Error"))) +
  guides(color = guide_legend(nrow = 2, ncol = 2, title.position = "top", title.hjust = 0.5), fill = guide_legend(nrow = 2, ncol = 2, title.position = "top", title.hjust = 0.5)) +
  theme_classic() + theme(legend.text = element_text(size = 10), legend.position = "bottom", legend.justification = c("center"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6), axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5), axis.title = element_text(size = 16), axis.text.y = element_text(angle = 0, size = 10, hjust = 1))
#
error_legend_plot <- tooth_alpha_plot + theme(legend.position = "right", legend.justification = "center", legend.box.just = "center") +
  guides(fill_ramp = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 1))
error_legend_plot <- as_ggplot(get_legend(error_legend_plot))
observed_error_plot_nolegend <- observed_error_plot + rremove("legend")
estimated_error_plot_nolegend <- estimated_error_plot + rremove("legend")
lowerhalf_error_plot_nolegend <- lowerhalf_error_plot + rremove("legend")
error_summary_plot <- ggarrange(observed_error_plot_nolegend, estimated_error_plot_nolegend, lowerhalf_error_plot_nolegend, error_legend_plot, nrow = 2, ncol = 2)
error_summary_plot

#Isotopic profile comparisons####
#Note: you can change the "specimen_no" variable to look at a different specific tooth (range is 1-10)
specimen_no <- 2
title_message <- paste0("Rousay Sheep M2: ", rousay_data[, .N, .(Specimen_No, Specimen)][Specimen_No %in% specimen_no, Specimen])
#
full_subtitle_message <- paste0("Every Sample (N = ", rousay_data[, .(N_Samples = sum(!is.na(d18O))), .(Specimen_No)][Specimen_No %in% specimen_no, N_Samples], ")")
full_profile <- 
  ggplot(rousay_data[Specimen_No %in% specimen_no]) +
  add_isotopic_profile(posterior = rousay_post, specimen_no = specimen_no, col = "black") +
  geom_point(aes(x = Dist, y = d18O), col = "black") +
  add_d18O_CIs(data = rousay_data, posterior = rousay_post$d18O, specimen_no = specimen_no, col = "black") +
  coord_cartesian(xlim = c(40, 0), ylim = c(-10, 0), expand = F) +
  scale_linetype_manual(name = "", values = c("Isotopic Profile" = "solid", "Linear Trend" = "dotted"), labels = c("Isotopic Profile", "Linear Trend")) +
  scale_alpha_manual(name = "", values = c("Isotopic Profile" = 0.3, "Linear Trend" = 0.1), labels = c("Isotopic Profile", "Linear Trend")) +
  labs(title = NULL, subtitle = full_subtitle_message, x = "Distance from REJ (mm)", y = expression(paste(delta^{18}, "O (\u2030)"))) +
  theme_classic() + theme(panel.border = element_rect(fill = NA), legend.text = element_text(size = 10), legend.position = "bottom", legend.justification = c("center"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6), plot.title = element_text(size = 16), plot.subtitle = element_text(size = 12), axis.title = element_text(size = 12), axis.text = element_text(size = 10))
#
third_subtitle_message <- paste0("Every Third Sample (N = ", rousay_data[, .(N_Samples = sum(!is.na(third_d18O))), .(Specimen_No)][Specimen_No %in% specimen_no, N_Samples], ")")
third_profile <- 
  ggplot(rousay_data[Specimen_No %in% specimen_no, .(Dist, d18O = third_d18O)][!is.na(d18O)]) +
  add_isotopic_profile(posterior = rousay_third_missing_post, specimen_no = specimen_no, col = "blue") +
  add_isotopic_profile(posterior = rousay_post, specimen_no = specimen_no, col = "black", linear_trend = FALSE) +
  geom_point(aes(x = Dist, y = d18O, col = "Observed")) +
  add_d18O_CIs(data = rousay_data[, .(Specimen_No, Dist, d18O = third_d18O)][!is.na(d18O)], posterior = rousay_third_missing_post$d18O, specimen_no = specimen_no) +
  #add in the missing values
  geom_point(data = rousay_data[Specimen_No %in% specimen_no, .(Dist, d18O, third_d18O)][is.na(third_d18O)], aes(x = Dist, y = d18O, col = "Missing")) +
  add_d18O_CIs(data = rousay_data[, .(Specimen_No, Dist, d18O = third_d18O)][is.na(d18O)], posterior = rousay_third_missing_post$d18O_missing, specimen_no = specimen_no, col = "blue") +
  coord_cartesian(xlim = c(40, 0), ylim = c(-10, 0), expand = F) +
  scale_linetype_manual(name = "", values = c("Isotopic Profile" = "solid", "Linear Trend" = "dotted"), labels = c("Isotopic Profile", "Linear Trend")) +
  scale_alpha_manual(name = "", values = c("Isotopic Profile" = 0.3, "Linear Trend" = 0.1), labels = c("Isotopic Profile", "Linear Trend")) +
  scale_color_manual(name = "", values = c("Observed" = "black", "Missing" = "blue"), labels = c("Observed", "Missing")) +
  labs(title = NULL, subtitle = third_subtitle_message, x = "Distance from REJ (mm)", y = expression(paste(delta^{18}, "O (\u2030)"))) +
  theme_classic() + theme(panel.border = element_rect(fill = NA), legend.text = element_text(size = 10), legend.position = "bottom", legend.justification = c("center"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6), plot.title = element_text(size = 16), plot.subtitle = element_text(size = 12), axis.title = element_text(size = 12), axis.text = element_text(size = 10))
#
fourth_subtitle_message <- paste0("Every Fourth Sample (N = ", rousay_data[, .(N_Samples = sum(!is.na(fourth_d18O))), .(Specimen_No)][Specimen_No %in% specimen_no, N_Samples], ")")
fourth_profile <- 
  ggplot(rousay_data[Specimen_No %in% specimen_no, .(Dist, d18O = fourth_d18O)][!is.na(d18O)]) +
  add_isotopic_profile(posterior = rousay_fourth_missing_post, specimen_no = specimen_no, col = "red") +
  add_isotopic_profile(posterior = rousay_post, specimen_no = specimen_no, col = "black", linear_trend = FALSE) +
  geom_point(aes(x = Dist, y = d18O, col = "Observed")) +
  add_d18O_CIs(data = rousay_data[, .(Specimen_No, Dist, d18O = fourth_d18O)][!is.na(d18O)], posterior = rousay_fourth_missing_post$d18O, specimen_no = specimen_no) +
  #add in the missing values
  geom_point(data = rousay_data[Specimen_No %in% specimen_no, .(Dist, d18O, fourth_d18O)][is.na(fourth_d18O)], aes(x = Dist, y = d18O, col = "Missing"), col = "red") +
  add_d18O_CIs(data = rousay_data[, .(Specimen_No, Dist, d18O = fourth_d18O)][is.na(d18O)], posterior = rousay_fourth_missing_post$d18O_missing, specimen_no = specimen_no, col = "red") +
  coord_cartesian(xlim = c(40, 0), ylim = c(-10, 0), expand = F) +
  scale_linetype_manual(name = "", values = c("Isotopic Profile" = "solid", "Linear Trend" = "dotted"), labels = c("Isotopic Profile", "Linear Trend")) +
  scale_alpha_manual(name = "", values = c("Isotopic Profile" = 0.3, "Linear Trend" = 0.1), labels = c("Isotopic Profile", "Linear Trend")) +
  scale_color_manual(name = "", values = c("Observed" = "black", "Missing" = "red"), labels = c("Observed", "Missing")) +
  labs(title = NULL, subtitle = fourth_subtitle_message, x = "Distance from REJ (mm)", y = expression(paste(delta^{18}, "O (\u2030)"))) +
  theme_classic() + theme(panel.border = element_rect(fill = NA), legend.text = element_text(size = 10), legend.position = "bottom", legend.justification = c("center"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6), plot.title = element_text(size = 16), plot.subtitle = element_text(size = 12), axis.title = element_text(size = 12), axis.text = element_text(size = 10))
#
firsthalf_subtitle_message <- paste0("First Half Samples (N = ", rousay_data[, .(N_Samples = sum(!is.na(first_d18O))), .(Specimen_No)][Specimen_No %in% specimen_no, N_Samples], ")")
firsthalf_profile <- 
  ggplot(rousay_data[Specimen_No %in% specimen_no, .(Dist, d18O = first_d18O)][!is.na(d18O)]) +
  add_isotopic_profile(posterior = rousay_first_missing_post, specimen_no = specimen_no, col = "orange") +
  add_isotopic_profile(posterior = rousay_post, specimen_no = specimen_no, linear_trend = FALSE, col = "black") +
  geom_point(aes(x = Dist, y = d18O, col = "Observed"), col = "black") +
  add_d18O_CIs(data = rousay_data[, .(Specimen_No, Dist, d18O = first_d18O)][!is.na(d18O)], posterior = rousay_first_missing_post$d18O, specimen_no = specimen_no, col = "black") +
  #add in the missing values
  geom_point(data = rousay_data[Specimen_No %in% specimen_no, .(Dist, d18O, first_d18O)][is.na(first_d18O)], aes(x = Dist, y = d18O, col = "Missing")) +
  add_d18O_CIs(data = rousay_data[, .(Specimen_No, Dist, d18O = first_d18O)][is.na(d18O)], posterior = rousay_first_missing_post$d18O_missing, specimen_no = specimen_no, col = "orange") +
  coord_cartesian(xlim = c(40, 0), ylim = c(-10, 0), expand = F) +
  scale_linetype_manual(name = "", values = c("Isotopic Profile" = "solid", "Linear Trend" = "dotted"), labels = c("Isotopic Profile", "Linear Trend")) +
  scale_alpha_manual(name = "", values = c("Isotopic Profile" = 0.3, "Linear Trend" = 0.1), labels = c("Isotopic Profile", "Linear Trend")) +
  scale_color_manual(name = "", values = c("Observed" = "black", "Missing" = "orange"), labels = c("Observed", "Missing")) +
  labs(title = NULL, subtitle = firsthalf_subtitle_message, x = "Distance from REJ (mm)", y = expression(paste(delta^{18}, "O (\u2030)"))) +
  theme_classic() + theme(panel.border = element_rect(fill = NA), legend.text = element_text(size = 10), legend.position = "bottom", legend.justification = c("center"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6), plot.title = element_text(size = 16), plot.subtitle = element_text(size = 12), axis.title = element_text(size = 12), axis.text = element_text(size = 10))
#
#remove the legend to put in the sampling intensity
profile_legend_plot <- as_ggplot(get_legend(one_row_legend_plot))
full_profile_nolegend <- full_profile + rremove("legend")
third_profile_nolegend <- third_profile + rremove("legend")
fourth_profile_nolegend <- fourth_profile + rremove("legend")
firsthalf_profile_nolegend <- firsthalf_profile + rremove("legend")
profile_plots <- ggarrange(full_profile_nolegend, third_profile_nolegend, fourth_profile_nolegend, firsthalf_profile_nolegend, nrow = 2, ncol = 2)
profile_plots <- ggarrange(profile_plots, profile_legend_plot, ncol = 1, nrow = 2, heights = c(0.9, 0.1))
profile_comparison <- annotate_figure(profile_plots, top = title_message)
profile_comparison

#Overall summary comparison####
full_mean_profile <- 
  ggplot(rousay_data[Specimen_No %in% specimen_no, .(Dist, d18O = fourth_d18O)][!is.na(d18O)]) +
  add_site_isotopic_profile(posterior = rousay_post, col = "black") +
  coord_cartesian(xlim = c(40, 0), ylim = c(-7, -1), expand = F) +
  scale_linetype_manual(name = "", values = c("Isotopic Profile" = "solid", "Linear Trend" = "dotted"), labels = c("Isotopic Profile", "Linear Trend")) +
  scale_alpha_manual(name = "", values = c("Isotopic Profile" = 0.3, "Linear Trend" = 0.1), labels = c("Isotopic Profile", "Linear Trend")) +
  labs(title = NULL, subtitle = "Every Sample\n(Average)", x = "Distance from REJ (mm)", y = expression(paste(delta^{18}, "O (\u2030)"))) +
  theme_classic() + theme(panel.border = element_rect(fill = NA), legend.text = element_text(size = 10), legend.position = "bottom", legend.justification = c("center"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6), plot.title = element_text(size = 16), plot.subtitle = element_text(size = 12), axis.title = element_text(size = 12), axis.text = element_text(size = 10))
third_mean_profile <- 
  ggplot(rousay_data[Specimen_No %in% specimen_no, .(Dist, d18O = fourth_d18O)][!is.na(d18O)]) +
  add_site_isotopic_profile(posterior = rousay_third_missing_post, col = "blue") +
  coord_cartesian(xlim = c(40, 0), ylim = c(-7, -1), expand = F) +
  scale_linetype_manual(name = "", values = c("Isotopic Profile" = "solid", "Linear Trend" = "dotted"), labels = c("Isotopic Profile", "Linear Trend")) +
  scale_alpha_manual(name = "", values = c("Isotopic Profile" = 0.3, "Linear Trend" = 0.1), labels = c("Isotopic Profile", "Linear Trend")) +
  labs(title = NULL, subtitle = "Every Third Sample\n(Average)", x = "Distance from REJ (mm)", y = expression(paste(delta^{18}, "O (\u2030)"))) +
  theme_classic() + theme(panel.border = element_rect(fill = NA), legend.text = element_text(size = 10), legend.position = "bottom", legend.justification = c("center"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6), plot.title = element_text(size = 16), plot.subtitle = element_text(size = 12), axis.title = element_text(size = 12), axis.text = element_text(size = 10))
fourth_mean_profile <- 
  ggplot(rousay_data[Specimen_No %in% specimen_no, .(Dist, d18O = fourth_d18O)][!is.na(d18O)]) +
  add_site_isotopic_profile(posterior = rousay_fourth_missing_post, col = "red") +
  coord_cartesian(xlim = c(40, 0), ylim = c(-7, -1), expand = F) +
  scale_linetype_manual(name = "", values = c("Isotopic Profile" = "solid", "Linear Trend" = "dotted"), labels = c("Isotopic Profile", "Linear Trend")) +
  scale_alpha_manual(name = "", values = c("Isotopic Profile" = 0.3, "Linear Trend" = 0.1), labels = c("Isotopic Profile", "Linear Trend")) +
  labs(title = NULL, subtitle = "Every Fourth Sample\n(Average)", x = "Distance from REJ (mm)", y = expression(paste(delta^{18}, "O (\u2030)"))) +
  theme_classic() + theme(panel.border = element_rect(fill = NA), legend.text = element_text(size = 10), legend.position = "bottom", legend.justification = c("center"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6), plot.title = element_text(size = 16), plot.subtitle = element_text(size = 12), axis.title = element_text(size = 12), axis.text = element_text(size = 10))
first_mean_profile <- 
  ggplot(rousay_data[Specimen_No %in% specimen_no, .(Dist, d18O = fourth_d18O)][!is.na(d18O)]) +
  add_site_isotopic_profile(posterior = rousay_first_missing_post, col = "orange") +
  coord_cartesian(xlim = c(40, 0), ylim = c(-7, -1), expand = F) +
  scale_linetype_manual(name = "", values = c("Isotopic Profile" = "solid", "Linear Trend" = "dotted"), labels = c("Isotopic Profile", "Linear Trend")) +
  scale_alpha_manual(name = "", values = c("Isotopic Profile" = 0.3, "Linear Trend" = 0.1), labels = c("Isotopic Profile", "Linear Trend")) +
  labs(title = NULL, subtitle = "Lower Half Samples\n(Average)", x = "Distance from REJ (mm)", y = expression(paste(delta^{18}, "O (\u2030)"))) +
  theme_classic() + theme(panel.border = element_rect(fill = NA), legend.text = element_text(size = 10), legend.position = "bottom", legend.justification = c("center"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6), plot.title = element_text(size = 16), plot.subtitle = element_text(size = 12), axis.title = element_text(size = 12), axis.text = element_text(size = 10))
mean_profile_comparison <- ggarrange(full_mean_profile, third_mean_profile, fourth_mean_profile, first_mean_profile, nrow = 1, ncol = 4, common.legend = T, legend = "none")
mean_profile_comparison

#Birth seasonality plot####
full_plot_data <- plot_data_creator(data = rousay_data, posterior = rousay_post)
rousay_seasonality_full_plot <-
  ggplot() +
  geom_arc(aes(x0 = 0, y0 = 0, r = 1, start = 0, end = 2 * pi), linetype = "dotted") +
  geom_segment(aes(x = -1, xend = 1, y = 0, yend = 0), linetype = "dotted") +
  geom_segment(aes(x = 0, xend = 0, y = -1, yend = 1), linetype = "dotted") +
  sapply(1:10, function(x) add_max_oxy_line(data = rousay_post$max_oxy[, x], radius = seq(0.75, 1.25, 0.05)[full_plot_data$Summary[order(Specimen_No), Median_Rank]][x], col = "grey20")) +
  add_max_oxy_line(data = rousay_post$mu_max_oxy, radius = seq(0.75, 1.25, 0.05)[full_plot_data$Summary[order(Specimen_No), Median_Rank]][11], col = "black") +
  #  geom_text_repel(data = full_plot_data$Label_Plotting[is.na(Specimen_No) | (Specimen_No %in% 1:10)], aes(x = Radius * sin(Max), y = Radius * cos(Max), label = Specimen), max.overlaps = Inf, min.segment.length = 0, size = 3) +
  geom_text_repel(data = full_plot_data$Label_Plotting[is.na(Specimen_No) | (Specimen_No %in% 11)], aes(x = Radius * sin(Min), y = Radius * cos(Min), label = Specimen), max.overlaps = Inf, col = "black", min.segment.length = 0) +
  labs(title = "Estimated Birth Seasonality: Rousay Sheep M2", subtitle = paste0("Every Sample (Average = ", rousay_sample_table[Sampling %in% "Full", Mean], ")"), x = "", y = "") +
  scale_color_manual(labels = c("Specimen", "Herd Average"), values = c("grey", "black")) +
  theme_classic() + theme(legend.text = element_text(size = 10), legend.position = "bottom", legend.justification = c("center"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6), axis.title = element_blank(), axis.text = element_blank())
#
rousay_seasonality_comparison_plot <-
  ggplot() +
  geom_arc(aes(x0 = 0, y0 = 0, r = 1, start = 0, end = 2 * pi), linetype = "dotted") +
  geom_segment(aes(x = -1, xend = 1, y = 0, yend = 0), linetype = "dotted") +
  geom_segment(aes(x = 0, xend = 0, y = -1, yend = 1), linetype = "dotted") +
  add_max_oxy_line(data = rousay_post$mu_max_oxy, radius = 0.95, col = "black") +
  add_max_oxy_line(data = rousay_third_missing_post$mu_max_oxy, radius = 1.00, col = "blue") +
  add_max_oxy_line(data = rousay_fourth_missing_post$mu_max_oxy, radius = 1.05, col = "red") +
  add_max_oxy_line(data = rousay_first_missing_post$mu_max_oxy, radius = 1.10, col = "orange") +
  #  geom_text_repel(data = data.table(Sampling = c("Every Sample", "Every Third", "Every Fourth", "First Half"), Radius = c(0.95, 1, 1.05, 1.10), rbind(max_oxy_quantiles(rousay_post$mu_max_oxy, ci = 0.95), max_oxy_quantiles(rousay_third_missing_post$mu_max_oxy, ci = 0.95), max_oxy_quantiles(rousay_fourth_missing_post$mu_max_oxy, ci = 0.95), max_oxy_quantiles(rousay_first_missing_post$mu_max_oxy, ci = 0.95))), aes(x = Radius * sin(Min), y = Radius * cos(Min), label = Sampling), max.overlaps = Inf, col = "black", min.segment.length = 0, nudge_x = -0.2, nudge_y = 0.25) +
  geom_arc_bar(data = max_oxy_quantiles(rousay_post$mu_max_oxy, ci = 0.95), aes(x0 = 0, y0 = 0, r0 = 0, r = 1.15, start = Min, end = Max), alpha = 0.05, fill = "black", col = "black", linetype = "dashed") +
  labs(title = "Comparison of Herd Averages", subtitle = "Variation in Sampling Intensity", x = "", y = "") +
  theme_classic() + theme(legend.text = element_text(size = 10), legend.position = "bottom", legend.justification = c("center"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6), axis.title = element_blank(), axis.text = element_blank())
#
rousay_seasonality_sampled_comparison_plot <-
  ggplot() +
  geom_arc(aes(x0 = 0, y0 = 0, r = 1, start = 0, end = 2 * pi), linetype = "dotted") +
  geom_segment(aes(x = -1, xend = 1, y = 0, yend = 0), linetype = "dotted") +
  geom_segment(aes(x = 0, xend = 0, y = -1, yend = 1), linetype = "dotted") +
  add_max_oxy_line(data = rousay_post$sampled_max_oxy, radius = 0.95, col = "black") +
  add_max_oxy_line(data = rousay_third_missing_post$sampled_max_oxy, radius = 1.00, col = "blue") +
  add_max_oxy_line(data = rousay_fourth_missing_post$sampled_max_oxy, radius = 1.05, col = "red") +
  add_max_oxy_line(data = rousay_first_missing_post$sampled_max_oxy, radius = 1.10, col = "orange") +
  #  geom_text_repel(data = data.table(Sampling = c("Every Sample", "Every Third", "Every Fourth", "First Half"), Radius = c(0.95, 1, 1.05, 1.10), rbind(max_oxy_quantiles(rousay_post$sampled_max_oxy, ci = 0.95), max_oxy_quantiles(rousay_third_missing_post$sampled_max_oxy, ci = 0.95), max_oxy_quantiles(rousay_fourth_missing_post$sampled_max_oxy, ci = 0.95), max_oxy_quantiles(rousay_first_missing_post$sampled_max_oxy, ci = 0.95))), aes(x = Radius * sin(Min), y = Radius * cos(Min), label = Sampling), max.overlaps = Inf, col = "black", min.segment.length = 0, nudge_x = -0.2, nudge_y = 0.25) +
  geom_arc_bar(data = max_oxy_quantiles(rousay_post$sampled_max_oxy, ci = 0.95), aes(x0 = 0, y0 = 0, r0 = 0, r = 1.15, start = Min, end = Max), alpha = 0.05, fill = "black", col = "black", linetype = "dashed") +
  labs(title = "Comparison of Herd Ranges", subtitle = "Variation in Sampling Intensity", x = "", y = "") +
  theme_classic() + theme(legend.text = element_text(size = 10), legend.position = "bottom", legend.justification = c("center"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6), axis.title = element_blank(), axis.text = element_blank())
#
rousay_seasonality_legend_plot <- as_ggplot(get_legend(legend_plot))
#
rousay_seasonality_plot <- ggarrange(rousay_seasonality_full_plot, rousay_seasonality_comparison_plot, rousay_seasonality_sampled_comparison_plot, rousay_seasonality_legend_plot, ncol = 2, nrow = 2)
rousay_seasonality_plot

#####Copyright Statement#####
#Copyright (c) <2023> <Jesse Wolfhagen>

#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included
#in all copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.
#####