functions {
#include hghp_funs.stan
}
data {
  int<lower=0>                    n;                      // number of observations
  int<lower=0>                    p;                      // total number of covariates (incl. intercept)
  int<lower=0>                    G;                      // number of factors to perform variable selection
  int<lower=0>                    psel;                   // number of the p covariates to perform variable selection
  int<lower=0>                    pforce;                 // number of the p covariates to perform variable selection
  int<lower=0>                    nlevelparms;            // sum of (number of levels in each group - 1)
  array[psel] int<lower=0>        indx_select;            // indices pertaining to variable selection
  array[pforce] int<lower=0>      indx_force;             // indices pertaining to covariates forced into the model
  array[p,G] int<lower=0>         lvl_indx;               // Gxp array consisting of 0s (doesn't belong to group) and (1, ..., Lg-1)'s (levels of the group where Lg = number of levels)
  array[nlevelparms] int<lower=1> lvlgrp_indx;            // nlevelparms array: lvlgrp_indx[j] is the group to which the level belongs. Used for shrinkage at the level-level.
  vector[n]                       y;                      // response vector (binary responses)
  matrix[n,p]                     X;                      // desgin matrix
  real<lower=0>                   beta_force_sd;          // SD for normal prior on forced-in betas
  real<lower=0>                   sigma_mean;             // mean parameter for normal prior on standard deviation
  real<lower=0>                   sigma_sd;               // std. dev. parameter for normal prior on standard deviation
}
transformed data {
  vector[pforce] sd_force = rep_vector(beta_force_sd, pforce);
  matrix[n,p] X_std;
  vector[p] sd_X;
  sd_X[indx_force] = rep_vector(1, pforce);  // set forced-terms SD = 1
  X_std[,indx_force] = X[,indx_force];     // Keep forced-terms as-is
  for ( j in indx_select ) {
    sd_X[j] = sd(X[, j]);
    X_std[, j] = X[, j] / sd_X[j];
  }
}
parameters {
  vector[p]                                 beta_raw_std; // raw regression coefficients: beta_raw_std ~ N(0, 1)
  real<lower=0,upper=pi()/2>                tau_raw;      // raw global shrinkage param: tau_raw ~ U(0,pi/2)
  vector<lower=0,upper=pi()/2>[psel]        lambda_raw;   // raw individual var selection params: lambda_raw ~ U(0,pi/2)
  vector<lower=0,upper=pi()/2>[G]           deltagrp_raw; // raw group var selection params: deltagrp_raw ~ U(0,pi/2)
  vector<lower=0,upper=pi()/2>[nlevelparms] deltalvl_raw; // raw level var selection params: deltalvl_raw ~ U(0,pi/2)
  real<lower=0>                             sigma;        // standard deviation
}
transformed parameters {
  real                tau       = tan(tau_raw);         // tau       ~ HC(0, 1)
  vector[psel]        lambda    = tan(lambda_raw);      // lambda    ~ HC(0, 1)
  vector[G]           delta_grp = tan(deltagrp_raw);    // delta_grp ~ HC(0, 1)
  
  // delta_lvl ~ HC(0, delta_grp)
  vector[nlevelparms] delta_lvl = get_delta_lvl(
    deltalvl_raw, delta_grp, lvlgrp_indx, nlevelparms
  );
  // Obtain prior SD for regression coefficients
  vector[p] sd_beta = get_beta_sd(
    delta_lvl, tau, lambda, lvl_indx, indx_select, indx_force, beta_force_sd
    , p, sd_force, G
  );
  // Normal prior on beta_std: 
  //   beta_raw_std ~ N(0, 1) => beta_std ~ N(0, sd_beta^2)
  vector[p] beta_std = sd_beta .* beta_raw_std; 
}
model {
  beta_raw_std  ~ std_normal();
  sigma ~ normal(sigma_mean, sigma_sd);
  
  // likelihood
  y ~ normal_id_glm(X_std, 0, beta_std, sigma);
}
generated quantities {
  vector[p] beta = inv(sd_X) .* beta_std;
}
