data {
    int<lower = 1> N_Specimens;
    int<lower = 1> N_Samples;
    array[N_Samples] int<lower = 0> Specimen;
    //observations
    array[N_Samples] real<lower = 0> dist_obs;
    array[N_Samples] real oxy_obs;
    array[N_Samples] real<lower = 0> dist_sd;
    array[N_Samples] real<lower = 0> oxy_sd;
    //hyper-parameters
    array[2] real<lower = 0> prior_period;
    real<lower = 0> period_low;
    real<lower = period_low> period_high;
}
parameters {
    //observations
    vector[N_Samples] d18O;
    vector<lower = 0>[N_Samples] dist;
    //hyper-parameters (among all teeth)
    real<lower = period_low, upper = period_high> mu_period_oxy;
    real mu_alpha_oxy;
    real mu_beta_oxy;
    real mu_Obeta1;
    real mu_Obeta2;
    //variation between teeth
    vector<lower = 0>[5] sigma_specimen;
    matrix[5, N_Specimens] z_specimen;
    cholesky_factor_corr[5] L_Rho_specimen;
    //sampling variability
    real<lower = 0> sigmasample_O;
}
transformed parameters {
  matrix[N_Specimens, 5] v_specimen;
  matrix[5, 5] Rho_specimen;
  real transformed_mu_period_oxy;
  //tooth-level parameters
  vector[N_Specimens] period_oxy;
  vector[N_Specimens] alpha_oxy;
  vector[N_Specimens] beta_oxy;
  vector[N_Specimens] Obeta1;
  vector[N_Specimens] Obeta2;
  //between-tooth variability
  v_specimen = (diag_pre_multiply(sigma_specimen, L_Rho_specimen) * z_specimen)';
  Rho_specimen = L_Rho_specimen * L_Rho_specimen';
  //tooth-level parameters
  transformed_mu_period_oxy = logit((mu_period_oxy - period_low) / (period_high - period_low));
  period_oxy = period_low + ((period_high - period_low) * inv_logit(transformed_mu_period_oxy + col(v_specimen, 1)));
  alpha_oxy = mu_alpha_oxy + col(v_specimen, 2);
  beta_oxy = mu_beta_oxy + col(v_specimen, 3);
  Obeta1 = mu_Obeta1 + col(v_specimen, 4);
  Obeta2 = mu_Obeta2 + col(v_specimen, 5);
}
model {
    real d18O_mu;
    mu_period_oxy ~ normal(prior_period[1], prior_period[2]) T[period_low, period_high];
    mu_alpha_oxy ~ normal(-5, 2.5);
    mu_beta_oxy ~ normal(0, 0.05);
    mu_Obeta1 ~ normal(0, 1.5);
    mu_Obeta2 ~ normal(0, 1.5);
    //variation between specimens
    sigma_specimen[1] ~ normal(0, 0.5); //bounded_period_oxy
    sigma_specimen[2] ~ normal(0, 1); //alpha_oxy
    sigma_specimen[3] ~ normal(0, 0.01); //beta_oxy
    sigma_specimen[4:5] ~ normal(0, 1); //Obeta1, Obeta2
    to_vector(z_specimen) ~ normal(0, 1);
    L_Rho_specimen ~ lkj_corr_cholesky(2);
    //sampling variability
    sigmasample_O ~ normal(0.08, 0.025);
    //observation model
    dist_obs ~ normal(dist, dist_sd);
    for(i in 1:N_Samples) {
      oxy_obs[i] ~ normal(d18O[i], oxy_sd[i]);
    }
    for(i in 1:N_Samples) {
      //matches formula in Chazin, et al. (2018) Formula 3.2 (with slope term added) rather than Stolwijk, et al. (1999)
    	d18O_mu = (alpha_oxy[Specimen[i]] + beta_oxy[Specimen[i]] * dist[i]) + Obeta1[Specimen[i]] * cos(2 * pi() * (dist[i] / period_oxy[Specimen[i]])) + Obeta2[Specimen[i]] * sin(2 * pi() * (dist[i] / period_oxy[Specimen[i]]));
    	//sampling variability
    	d18O[i] ~ normal(d18O_mu, sigmasample_O);
    }
}
generated quantities {
  //herd-level parameters
  real mu_mean_oxy;
  real mu_amplitude_oxy;
  real mu_max_oxy_dist;
  real mu_min_oxy_dist;
  real mu_max_oxy_radian;
  real mu_max_oxy_d18O; //d18O value at max_oxy
  real mu_min_oxy_d18O; //d18O value at min_oxy point (max_oxy + 0.5 * periodicity)
  real mu_mid_oxy_d18O; //d18O value at mid-poitn between max_oxy and min_oxy (max_oxy + 0.25 * periodicity)
  //specimen-level parameters
  vector[N_Specimens] amplitude_oxy;
  vector[N_Specimens] mean_oxy;
  vector[N_Specimens] max_oxy_dist;
  vector[N_Specimens] min_oxy_dist;
  vector[N_Specimens] max_oxy_radian;
  vector[N_Specimens] max_oxy_d18O; //d18O value at max_oxy
  vector[N_Specimens] min_oxy_d18O; //d18O value at min_oxy point (max_oxy + 0.5 * periodicity)
  vector[N_Specimens] mid_oxy_d18O; //d18O value at mid-poitn between max_oxy and min_oxy (max_oxy + 0.25 * periodicity)
  //sampling the parameters to get a set of uncertainty about the population
  matrix[5, 1] z_sampled;
  matrix[1, 5] v_sampled;
  real sampled_period_oxy;
  real sampled_alpha_oxy;
  real sampled_beta_oxy;
  real sampled_Obeta1;
  real sampled_Obeta2;
  real sampled_max_oxy_dist;
  real sampled_min_oxy_dist;
  real sampled_max_oxy_radian;
  real sampled_max_oxy_d18O;
  real sampled_min_oxy_d18O;
  real sampled_mid_oxy_d18O;
  
  //herd-level parameters for average values (mean_oxy, amplitude_oxy, max_oxy_dist, min_oxy_dist, max_oxy_radian)
  mu_mean_oxy = ((mu_alpha_oxy + (mu_beta_oxy * mu_period_oxy)) + (mu_alpha_oxy + (mu_beta_oxy * 0))) / 2;
  mu_amplitude_oxy = sqrt(pow(mu_Obeta1, 2) + pow(mu_Obeta2, 2));
  if(mu_Obeta2 / mu_Obeta1 > 0 && mu_Obeta2 > 0) {
    mu_max_oxy_dist = (atan(mu_Obeta2 / mu_Obeta1) * (mu_period_oxy / (2 * pi()))) + 0;
    mu_min_oxy_dist = mu_max_oxy_dist + (mu_period_oxy / 2);
  }
  if(mu_Obeta2 / mu_Obeta1 > 0 && mu_Obeta2 <= 0) {
    mu_max_oxy_dist = (atan(mu_Obeta2 / mu_Obeta1) * (mu_period_oxy / (2 * pi()))) + (mu_period_oxy / 2);
    mu_min_oxy_dist = mu_max_oxy_dist - (mu_period_oxy / 2);
  }
  if(mu_Obeta2 / mu_Obeta1 <= 0 && mu_Obeta2 > 0) {
    mu_max_oxy_dist = (atan(mu_Obeta2 / mu_Obeta1) * (mu_period_oxy / (2 * pi()))) + (mu_period_oxy / 2);
    mu_min_oxy_dist = mu_max_oxy_dist + (mu_period_oxy / 2);
  }
  if(mu_Obeta2 / mu_Obeta1 <= 0 && mu_Obeta2 <= 0) {
    mu_max_oxy_dist = (atan(mu_Obeta2 / mu_Obeta1) * (mu_period_oxy / (2 * pi()))) + (mu_period_oxy);
    mu_min_oxy_dist = mu_max_oxy_dist - (mu_period_oxy / 2);
  }
  mu_max_oxy_radian = (mu_max_oxy_dist / mu_period_oxy) * 2 * pi();
  //calculations for curve-related d18O values (max_d18O, min_d180, mean_d18O)
  mu_max_oxy_d18O = mu_alpha_oxy + (mu_beta_oxy * mu_max_oxy_dist) + mu_Obeta1 * cos(2 * pi() * (mu_max_oxy_dist / mu_period_oxy)) + mu_Obeta2 * sin(2 * pi() * (mu_max_oxy_dist / mu_period_oxy));
  mu_min_oxy_d18O = mu_alpha_oxy + (mu_beta_oxy * mu_min_oxy_dist) + mu_Obeta1 * cos(2 * pi() * (mu_min_oxy_dist / mu_period_oxy)) + mu_Obeta2 * sin(2 * pi() * (mu_min_oxy_dist / mu_period_oxy));
  mu_mid_oxy_d18O = mu_alpha_oxy + (mu_beta_oxy * (mu_max_oxy_dist + mu_min_oxy_dist) / 2);

  //specimen-level parameters (mean_oxy, amplitude_oxy, max_oxy_dist, min_oxy_dist, max_oxy_radian)
  for(i in 1:N_Specimens) {
    mean_oxy[i] = ((alpha_oxy[i] + beta_oxy[i] * period_oxy[i]) + (alpha_oxy[i] + beta_oxy[i] * 0)) / 2;
    amplitude_oxy[i] = sqrt(pow(Obeta1[i], 2) + pow(Obeta2[i], 2));
    if(Obeta2[i] / Obeta1[i] > 0 && Obeta2[i] > 0) {
      max_oxy_dist[i] = (atan(Obeta2[i] / Obeta1[i]) * (period_oxy[i] / (2 * pi()))) + 0;
      min_oxy_dist[i] = max_oxy_dist[i] + (period_oxy[i] / 2);
    }
    if(Obeta2[i] / Obeta1[i] > 0 && Obeta2[i] <= 0) {
      max_oxy_dist[i] = (atan(Obeta2[i] / Obeta1[i]) * (period_oxy[i] / (2 * pi()))) + (period_oxy[i] / 2);
      min_oxy_dist[i] = max_oxy_dist[i] - (period_oxy[i] / 2);
    }
    if(Obeta2[i] / Obeta1[i] <= 0 && Obeta2[i] > 0) {
      max_oxy_dist[i] = (atan(Obeta2[i] / Obeta1[i]) * (period_oxy[i] / (2 * pi()))) + (period_oxy[i] / 2);
      min_oxy_dist[i] = max_oxy_dist[i] + (period_oxy[i] / 2);
    }
    if(Obeta2[i] / Obeta1[i] <= 0 && Obeta2[i] <= 0) {
      max_oxy_dist[i] = (atan(Obeta2[i] / Obeta1[i]) * (period_oxy[i] / (2 * pi()))) + (period_oxy[i]);
      min_oxy_dist[i] = max_oxy_dist[i] - (period_oxy[i] / 2);
    }
    max_oxy_radian[i] = (max_oxy_dist[i] / period_oxy[i]) * 2 * pi();
    //specimen-level calculations for curve-related d18O values (max_d18O, min_d180, mean_d18O)
    max_oxy_d18O[i] = alpha_oxy[i] + (beta_oxy[i] * max_oxy_dist[i]) + Obeta1[i] * cos(2 * pi() * (max_oxy_dist[i] / period_oxy[i])) + Obeta2[i] * sin(2 * pi() * (max_oxy_dist[i] / period_oxy[i]));
    min_oxy_d18O[i] = alpha_oxy[i] + (beta_oxy[i] * min_oxy_dist[i]) + Obeta1[i] * cos(2 * pi() * (min_oxy_dist[i] / period_oxy[i])) + Obeta2[i] * sin(2 * pi() * (min_oxy_dist[i] / period_oxy[i]));
    mid_oxy_d18O[i] = alpha_oxy[i] + (beta_oxy[i] * (max_oxy_dist[i] + min_oxy_dist[i]) / 2);
  }
  
  //values for a sampled specimen (uncertainty: sampled_max_oxy)
  for(i in 1:5) {
    z_sampled[i, 1] = normal_rng(0, 1);
  }
  v_sampled = (diag_pre_multiply(sigma_specimen, L_Rho_specimen) * z_sampled)';
  sampled_period_oxy = period_low + ((period_high - period_low) * inv_logit(transformed_mu_period_oxy + v_sampled[1, 1]));
  sampled_alpha_oxy = mu_alpha_oxy + v_sampled[1, 2];
  sampled_beta_oxy = mu_beta_oxy + v_sampled[1, 3];
  sampled_Obeta1 = mu_Obeta1 + v_sampled[1, 4];
  sampled_Obeta2 = mu_Obeta2 + v_sampled[1, 5];
  //calculate the sampled max_oxy value
  if(sampled_Obeta2 / sampled_Obeta1 > 0 && sampled_Obeta2 > 0) {
    sampled_max_oxy_dist = (atan(sampled_Obeta2 / sampled_Obeta1) * (sampled_period_oxy / (2 * pi()))) + 0;
    sampled_min_oxy_dist = sampled_max_oxy_dist + (sampled_period_oxy / 2);
  }
  if(sampled_Obeta2 / sampled_Obeta1 > 0 && sampled_Obeta2 <= 0) {
    sampled_max_oxy_dist = (atan(sampled_Obeta2 / sampled_Obeta1) * (sampled_period_oxy / (2 * pi()))) + (sampled_period_oxy / 2);
    sampled_min_oxy_dist = sampled_max_oxy_dist - (sampled_period_oxy / 2);
  }
  if(sampled_Obeta2 / sampled_Obeta1 <= 0 && sampled_Obeta2 > 0) {
    sampled_max_oxy_dist = (atan(sampled_Obeta2 / sampled_Obeta1) * (sampled_period_oxy / (2 * pi()))) + (sampled_period_oxy / 2);
    sampled_min_oxy_dist = sampled_max_oxy_dist + (sampled_period_oxy / 2);
  }
  if(sampled_Obeta2 / sampled_Obeta1 <= 0 && sampled_Obeta2 <= 0) {
    sampled_max_oxy_dist = (atan(sampled_Obeta2 / sampled_Obeta1) * (sampled_period_oxy / (2 * pi()))) + (sampled_period_oxy);
    sampled_min_oxy_dist = sampled_max_oxy_dist - (sampled_period_oxy / 2);
  }
  sampled_max_oxy_radian = (sampled_max_oxy_dist / sampled_period_oxy) * 2 * pi();
  sampled_max_oxy_d18O = sampled_alpha_oxy + (sampled_beta_oxy * sampled_max_oxy_dist) + sampled_Obeta1 * cos(2 * pi() * (sampled_max_oxy_dist / sampled_period_oxy)) + sampled_Obeta2 * sin(2 * pi() * (sampled_max_oxy_dist / sampled_period_oxy));
  sampled_min_oxy_d18O = sampled_alpha_oxy + (sampled_beta_oxy * sampled_min_oxy_dist) + sampled_Obeta1 * cos(2 * pi() * (sampled_min_oxy_dist / sampled_period_oxy)) + sampled_Obeta2 * sin(2 * pi() * (sampled_min_oxy_dist / sampled_period_oxy));
  sampled_mid_oxy_d18O = sampled_alpha_oxy + (sampled_beta_oxy * (sampled_max_oxy_dist + sampled_min_oxy_dist) / 2);
}
