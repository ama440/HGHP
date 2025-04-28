// This file is a modification to the bayesian_lasso_mod.stan file to standardize covariates first
// I also leave the intercept unshrunk now
data {
  int<lower=0>                            n;                      // number of observations
  int<lower=0>                            p;                      // total number of covariates (incl. intercept)
  array[n] int<lower=0, upper=1>          y;                      // response vector (binary responses)
  matrix[n,p]                             X;                      // design matrix
  int<lower=0, upper=1>                   force_trt;              // indicator for forcing in treatment  
  array[force_trt+1] int<lower=1, upper=p>  index_force;          // Indices for noninformative prior on intercept and optionally treatment
  array[force_trt ? p-2: p-1] int<lower=1, upper=p>  index_sel;   // Index of covariates being shrunk
}
transformed data {
  matrix[n,p] X_std;
  vector[p] sd_X;
  
  sd_X[index_force] = rep_vector(1, force_trt+1);   // set forced-terms SD = 1
  X_std[,index_force] = X[,index_force];            // Keep forced-terms as-is
  for (j in index_sel) {
    if (j > 1) {
      sd_X[j] = sd(X[, j]);
      X_std[, j] = X[, j] / sd_X[j];
    }
  }
}
parameters {
  vector[p] beta_std;                                 // Coefficients
  real<lower=0,upper=pi()/2> lasso_scale_uniform;     // Shrinkage parameter
}
transformed parameters {
  real lasso_scale = tan(lasso_scale_uniform);
}
model {
  // prior
  if (force_trt == 1) {
    // Noninformative priors on intercept and optionally beta_trt
    beta_std[index_force] ~ normal(0, 10);            
    
    // Rest of parameters get double exponential prior
    target += double_exponential_lpdf(beta_std[index_sel] | 0, lasso_scale);
  } else {
    target += double_exponential_lpdf(beta_std | 0, lasso_scale);                   // double exponential prior on all covariates
  }
  lasso_scale ~ cauchy(0,1);
  
  // likelihood
  y ~ bernoulli_logit_glm(X_std, 0, beta_std);
}
generated quantities {
  vector[p] beta = inv(sd_X) .* beta_std;
}
