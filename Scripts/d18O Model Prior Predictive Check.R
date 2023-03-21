library("boot")
library("circular")
library("data.table")
library("doParallel")
library("ggdist")
library("ggforce")
library("ggplot2")
library("ggpubr")
library("ggrepel")
library("parallel")
#The goal is to simulate an assemblage from the prior distributions to get a sense of their "reasonability"
#setting resources to do simulations with multicore processing (on Windows)
num_cores <- detectCores(logical = T)
cl <- makeCluster(num_cores - 4) #creating a virtual cluster
registerDoParallel(cl) #registering the cluster

#
blanz_etal_data <- fread("./Data/Blanz, et al. 2020 Sheep CO Data (Cleaned).csv")
rousay_data <- blanz_etal_data[Site %in% "ROU"]
#
rousay_data[, Specimen_No := as.numeric(as.factor(Specimen))]
rousay_data[, Sample_No := tstrsplit(Sample, "M2-", keep = 2)]
rousay_data[, Sample_No := as.numeric(Sample_No)]
#
rousay_data[, c("third_d18O", "fourth_d18O") := .(d18O, d18O)]
rousay_data[Sample_No %% 3 != 1, third_d18O := NA]
rousay_data[Sample_No %% 4 != 1, fourth_d18O := NA]
#
#first half of the samples (as if the tooth were more heavily worn)
rousay_data[, c("first_d18O") := .(d18O)]
rousay_data[, Max_Sample := .N, Specimen]
rousay_data[Sample_No / Max_Sample < 0.5, first_d18O := NA]
#
rousay_sample_summary <- data.table(Sampling = c("Full", "Third", "Fourth", "First Half"), rousay_data[, .(N_Samples = .N, Third_Samples = sum(!is.na(third_d18O)), Fourth_Samples = sum(!is.na(fourth_d18O)), First_Samples = sum(!is.na(first_d18O))), .(Specimen_No, Specimen)][, .(Mean = lapply(.SD, mean), Median = lapply(.SD, median), Min = lapply(.SD, min), Max = lapply(.SD, max)), .SDcols = c("N_Samples", "Third_Samples", "Fourth_Samples", "First_Samples")])

####Bring in published MLE fits to help inform prior distributions####
published_fits <- fread("./Published d18O Cyclical Regression Results.csv")
published_fits[Taxon %in% c("Ovis aries"), .(Mean = mean(Periodicity), SD = sd(Periodicity)), .(Tooth)]
published_fits[, .(Amplitude = quantile(Amplitude, c(0.025, 0.975)))]
#

#other observed data (to check model error)
previous_d18O_data <- rbind(
  rousay_data[, .(Site, Taxon = "Ovis aries", Specimen, Tooth, Dist, d18O)],
  fread("./Data/Blaise and Balasse 2011 - Modern Reference Sheep Data.csv")[, .(Site, Taxon, Specimen, Tooth, Dist, d18O)],
  fread("./Data/Balasse et al 2017 d18O Data.csv")[, .(Site, Taxon, Specimen, Tooth, Dist, d18O)],
  fread("./Data/Balasse, et al. 2012 - Bercy.csv")[, .(Site, Taxon, Specimen, Tooth, Dist, d18O)]
)
#
previous_d18O_fits <- published_fits[previous_d18O_data[, .N, .(Specimen, Tooth)], on = c("Specimen", "Tooth")][, .(Site, Specimen, Tooth, period_oxy = Periodicity, amplitude_oxy = Amplitude, mean_oxy = Mean_d18O, max_point = Critical_Value )][!is.na(Site)]

####Figure out parameters for model error (sigma_d18O) using the Rousay MLE fits####
balasse_model_fit <- function(data, par) {
  d18O_predict <- par[, amplitude_oxy] * cos(2 * pi * ((data[, Dist] - par[, max_point]) / par[, period_oxy])) + par[, mean_oxy]
  #  sum((d18O_predict - data[, d18O])^2)
  data.table(data[, .(Specimen, Dist, d18O)], d18O_predict = d18O_predict)
}
#

previous_d18O_mle_fits <- rbindlist(lapply(previous_d18O_fits[, Specimen], function(x) data.table(Site = previous_d18O_fits[Specimen %in% x, Site], balasse_model_fit(previous_d18O_data[Specimen %in% x], previous_d18O_fits[Specimen %in% x]))))
previous_d18O_mle_fits[, d18O_error := (d18O_predict - d18O)]
previous_d18O_mle_fits[Site %in% "Carmejane farm" == F, sd(d18O_error)]
#

#bootstrapping to get uncertainty on the estimate
mean(sapply(1:1000, function(x) previous_d18O_mle_fits[abs(d18O_error) <= 1][sample(1:.N, .N, replace = T)][, sd(d18O_error)]))
sd(sapply(1:1000, function(x) previous_d18O_mle_fits[abs(d18O_error) <= 1][sample(1:.N, .N, replace = T)][, sd(d18O_error)]))
#visualize
hist(rnorm(1e5, 0.25, 0.01), col = "white", border = NULL, density = 10, freq = F, ylim = c(0, 60), main = "Estimate of Sinusoidal Regression Model Error")
hist(sapply(1:1000, function(x) sd(previous_d18O_mle_fits[abs(d18O_error) <= 1][sample(1:.N, .N, replace = T)][, d18O_error])), freq = F, add = T)
hist(rnorm(1e5, 0.25, 0.01), col = "red", density = 10, freq = F, add = T)
legend(x = "topright", fill = c("gray", "red"), legend = c("Derived from Literature (Bootstrapped)", "Proposed Prior"), bty = "n")

####Prior Predictive Checking for the Sinsusoidal Regression Model####


singlesite_N_Specimens <- 10
singlesite_observed_dist <- seq(5, 30, 3)
singlesite_period_boundaries <- c(10, 70)
singlesite_d18O_error <- 0.1
singlesite_grand_prior_list <- list(
  prior_period_oxy = c(30, 3.5),
  prior_logperiod_oxy = c(3.3, 0.275),
  prior_alpha_oxy = c(-5, 2.5),
  prior_beta_oxy = c(0, 0.05),
  prior_Obeta1 = c(0, 1.5),
  prior_Obeta2 = c(0, 1.5)
)
singlesite_specimen_sigma_values <- c(
  0.25, #transformed_period_oxy
#  0.1, #period_log_oxy
#  2, #period_oxy
  1, #alpha_oxy
  0.01, #beta_oxy
  1, #Obeta1
  1 #Obeta2
)
singlesite_prior_sigmasample <- 0.15

#Functions
rnorm_bounded <- function(n, mu, sigma, lower, upper) {
  #samples from a defined normal distribution until it finds a positive value
  accept = FALSE
  while(!accept) {
    value <- rnorm(n, mean = mu, sd = sigma)
    if(sum(value > lower & value < upper) == length(value)) accept <- TRUE
  }
  value
}
rposnorm <- function(n, mu, sigma) {
  #samples from a defined normal distribution until it finds a positive value
  accept = FALSE
  while(!accept) {
    value <- rnorm(n, mean = mu, sd = sigma)
    if(sum(value > 0) == length(value)) accept <- TRUE
  }
  value
}
rho_matrix_maker <- function(lkjcorr, K) {
  Rho_matrix <- matrix(NA, ncol = K, nrow = K)
  for(i in 1:K) {
    Rho_matrix[i,] <- c(rep(1, K - i), lkjcorr, rep(1, i - 1))
  }
  Rho_matrix
}

singlesite_d18O_regression_prior_predictive_check <- function(N_Specimens = singlesite_N_Specimens, grand_prior = singlesite_grand_prior_list, specimen_sigma_vals = singlesite_specimen_sigma_values, period_bounds = singlesite_period_boundaries, prior_sigmasample = singlesite_prior_sigmasample, d18O_error = singlesite_d18O_error, observed_dist = singlesite_observed_dist) {
  #Simulating from the priors, start at the top
  mu_period_oxy = rnorm_bounded(1, grand_prior$prior_period_oxy[1], grand_prior$prior_period_oxy[2], period_bounds[1], period_bounds[2])
  mu_alpha_oxy = rnorm(1, grand_prior$prior_alpha_oxy[1], grand_prior$prior_alpha_oxy[2])
  mu_beta_oxy = rnorm(1, grand_prior$prior_beta_oxy[1], grand_prior$prior_beta_oxy[2])
  mu_Obeta1 = rnorm(1, grand_prior$prior_Obeta1[1], grand_prior$prior_Obeta1[2])
  mu_Obeta2 = rnorm(1, grand_prior$prior_Obeta2[1], grand_prior$prior_Obeta2[2])

  #Internal priors for the multilevel modeling
  specimen_sigma <- c(rposnorm(1, 0, specimen_sigma_vals[1]),
                      rposnorm(1, 0, specimen_sigma_vals[2]),
                      rposnorm(1, 0, specimen_sigma_vals[3]),
                      rposnorm(1, 0, specimen_sigma_vals[4]),
                      rposnorm(1, 0, specimen_sigma_vals[5]))
  #
  z_specimen <- matrix(rnorm(5 * N_Specimens, 0, 1), ncol = N_Specimens, nrow = 5)
  Rho_specimen <- rlkjcorr_marginal(1, K = 5, eta = 2)

  #Sampling variability (deviation from regression line)
  sigmasample_O = rposnorm(1, 0.08, 0.025)
  
  #Calculate a single set of the parameter values
  v_specimen <- t(diag(specimen_sigma) %*% rho_matrix_maker(Rho_specimen, 5) %*% diag(specimen_sigma) %*% z_specimen)
  
  #calculating the output variables for each specimen (creates the curves)
  #transforming period_oxy to be constrained within the limits (lower, upper)
  transformed_mu_period_oxy = boot::logit((mu_period_oxy - period_bounds[1]) / (period_bounds[2] - period_bounds[1]))
  period_oxy = period_bounds[1] + ((period_bounds[2] - period_bounds[1]) * boot::inv.logit(transformed_mu_period_oxy + v_specimen[, 1]))
  alpha_oxy = mu_alpha_oxy + v_specimen[, 2]
  beta_oxy = mu_beta_oxy + v_specimen[, 3]
  Obeta1 = mu_Obeta1 + v_specimen[, 4]
  Obeta2 = mu_Obeta2 + v_specimen[, 5]

  #use the curves to estimate values for a set of observed distance values
  estimated_d18O_values <- data.table(Specimen_No = rep(1:N_Specimens, each = length(observed_dist)), Modeled_Dist = rnorm(length(rep(observed_dist, N_Specimens)), rep(observed_dist, N_Specimens), 0.5))
  estimated_d18O_values[, Mu_d18O := (alpha_oxy[Specimen_No] + beta_oxy[Specimen_No] * Modeled_Dist) + Obeta1[Specimen_No] * cos(2 * pi * (Modeled_Dist / period_oxy[Specimen_No])) + Obeta2[Specimen_No] * sin(2 * pi * (Modeled_Dist / period_oxy[Specimen_No]))]

  #Sampling error (deviation from mean)
  estimated_d18O_values[, Modeled_d18O := rnorm(.N, Mu_d18O, rep(sigmasample_O, .N))]
  
  #Observation error (deviation from modeled value)
  estimated_d18O_values[, Observed_d18O := rnorm(.N, Modeled_d18O, d18O_error)]

  #
  list(Overall_Parameters = data.table(period_oxy = mu_period_oxy, alpha_oxy = mu_alpha_oxy, beta_oxy = mu_beta_oxy, Obeta1 = mu_Obeta1, Obeta2 = mu_Obeta2, sigmasample_O = sigmasample_O), Specimen_Parameters = cbind(period_oxy, alpha_oxy, beta_oxy, Obeta1, Obeta2), Observations = data.table(estimated_d18O_values))
}

singlesite_d18O_regression_prior_predictive_check()$Specimen_Parameters
#First, need to details about the assemblage

clusterExport(cl, list('data.table', 'rbindlist', 'rnorm_bounded', 'rposnorm', 'rho_matrix_maker', 'rlkjcorr_marginal', 'inv.logit', 'singlesite_d18O_regression_prior_predictive_check', 'singlesite_N_Specimens', 'singlesite_grand_prior_list', 'singlesite_specimen_sigma_values', 'singlesite_prior_sigmasample', 'singlesite_observed_dist', 'singlesite_d18O_error', 'singlesite_period_boundaries'))

singlesite_predicted_data <- parLapply(cl = cl, 1:1000, fun = function(x) singlesite_d18O_regression_prior_predictive_check(N_Specimens = singlesite_N_Specimens))
beepr::beep()

oxyfx <- function(params, dist) {
  period <- params[1]
  intercept <- params[2]
  slope <- params[3]
  beta1 <- params[4]
  beta2 <- params[5]
  (intercept + slope * dist) + (beta1 * cos(2 * pi * (dist / period))) + (beta2 * sin(2 * pi * (dist / period)))
}

#
prior_predict_parameters <- rbindlist(lapply(singlesite_predicted_data, function(x) data.table(x$Overall_Parameters)))
prior_predict_specimen_parameters <- rbindlist(lapply(singlesite_predicted_data, function(x) data.table(x$Specimen_Parameters)))
prior_predict_d18O <- rbindlist(lapply(singlesite_predicted_data, function(x) data.table(x$Observations)))

#
add_prior_predict_isotopic_profile <- function(posterior, col = "black") {
  isotopic_profile_df <- rbindlist(lapply(1:length(posterior$alpha_oxy), function(x) data.table(Iteration = x, Dist = seq(40, 0, length.out = 1000), m_d18O = oxyfx(params = c(posterior$period_oxy[x], posterior$alpha_oxy[x], posterior$beta_oxy[x], posterior$Obeta1[x], posterior$Obeta2[x]), dist = seq(40, 0, length.out = 1000)))))
  linear_trend_df <- rbindlist(lapply(1:length(posterior$alpha_oxy), function(x) data.table(Iteration = x, Dist = seq(40, 0, length.out = 1000), m_d18O_line = posterior$alpha_oxy[x] + posterior$beta_oxy[x] * seq(40, 0, length.out = 1000))))
  list(
    #add the isotopic profile (and CI)
    geom_line(data = isotopic_profile_df[, .(m_d18O = mean(m_d18O)), .(Dist)], aes(x = Dist, y = m_d18O, linetype = "Isotopic Profile"), col = col),
    geom_ribbon(data = isotopic_profile_df[, .(CI_low = quantile(m_d18O, 0.025), CI_high = quantile(m_d18O, 0.975)), Dist], aes(x = Dist, ymin = CI_low, ymax = CI_high, linetype = "Isotopic Profile", alpha = "Isotopic Profile"), fill = col, col = col),
    #add the linear trend (and CI)
    geom_line(data = linear_trend_df[, .(m_d18O_line = mean(m_d18O_line)), Dist], aes(x = Dist, y = m_d18O_line, linetype = "Linear Trend"), col = col),
    geom_ribbon(data = linear_trend_df[, .(CI_low = quantile(m_d18O_line, 0.025), CI_high = quantile(m_d18O_line, 0.975)), Dist], aes(x = Dist, ymin = CI_low, ymax = CI_high, linetype = "Linear Trend", alpha = "Linear Trend"), fill = col, col = col)
  )
}
add_prior_predict_isotopic_lines <- function(posterior, col = "black") {
  isotopic_profile_df <- rbindlist(lapply(1:length(posterior$alpha_oxy), function(x) data.table(Iteration = x, Dist = seq(40, 0, length.out = 1000), m_d18O = oxyfx(params = c(posterior$period_oxy[x], posterior$alpha_oxy[x], posterior$beta_oxy[x], posterior$Obeta1[x], posterior$Obeta2[x]), dist = seq(40, 0, length.out = 1000)))))
  list(
    #add the isotopic profile
    geom_line(data = isotopic_profile_df[, .(m_d18O = mean(m_d18O)), .(Dist)], aes(x = Dist, y = m_d18O, linetype = "Isotopic Profile"), col = col)
  )
}

#Overall Parameters
ggplot() +
  add_prior_predict_isotopic_profile(posterior = prior_predict_parameters) +
  coord_cartesian(xlim = c(40, 0), ylim = c(-20, 10), expand = F) +
  scale_linetype_manual(name = "", values = c("Isotopic Profile" = "solid", "Linear Trend" = "dotted"), labels = c("Isotopic Profile", "Linear Trend")) +
  scale_alpha_manual(name = "", values = c("Isotopic Profile" = 0.3, "Linear Trend" = 0.1), labels = c("Isotopic Profile", "Linear Trend")) +
  labs(title = "Prior Predictive Check", subtitle = "Averaged Parameters", x = "Distance from REJ (mm)", y = expression(paste(delta^{18}, "O (\u2030)"))) +
  theme_classic() + theme(legend.text = element_text(size = 10), legend.position = "bottom", legend.justification = c("center"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6), plot.title = element_text(size = 16), plot.subtitle = element_text(size = 12), axis.title = element_text(size = 12), axis.text = element_text(size = 10))
hist(prior_predict_parameters$period_oxy)
#
ggplot() +
  add_prior_predict_isotopic_profile(posterior = prior_predict_specimen_parameters) +
  add_prior_predict_isotopic_lines(posterior = prior_predict_specimen_parameters[1]) +
  add_prior_predict_isotopic_lines(posterior = prior_predict_specimen_parameters[((3 - 1) * singlesite_N_Specimens) + 1]) +
  add_prior_predict_isotopic_lines(posterior = prior_predict_specimen_parameters[((3 - 1) * singlesite_N_Specimens) + 1]) +
  add_prior_predict_isotopic_lines(posterior = prior_predict_specimen_parameters[((1000 - 1) * singlesite_N_Specimens) + 1]) +
  coord_cartesian(xlim = c(40, 0), ylim = c(-20, 10), expand = F) +
  scale_linetype_manual(name = "", values = c("Isotopic Profile" = "solid", "Linear Trend" = "dotted"), labels = c("Isotopic Profile", "Linear Trend")) +
  scale_alpha_manual(name = "", values = c("Isotopic Profile" = 0.3, "Linear Trend" = 0.1), labels = c("Isotopic Profile", "Linear Trend")) +
  labs(title = "Prior Predictive Check", subtitle = "Specimen-Specific Parameters", x = "Distance from REJ (mm)", y = expression(paste(delta^{18}, "O (\u2030)"))) +
  theme_classic() + theme(legend.text = element_text(size = 10), legend.position = "bottom", legend.justification = c("center"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6), plot.title = element_text(size = 16), plot.subtitle = element_text(size = 12), axis.title = element_text(size = 12), axis.text = element_text(size = 10))
hist(prior_predict_specimen_parameters$period_oxy)
