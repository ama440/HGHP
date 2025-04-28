data {
  int<lower=0> n;                      // number of observations
  int<lower=0> d;                      // number of predictors
  array[n] int<lower=0, upper=1> y;    // outputs
  matrix[n,d] x;                       // inputs, not including intercept column
  real<lower=0> scale_icept;           // prior std for the intercept
  real<lower=0> scale_global;          // scale for the half-t prior for tau
  real<lower=1> nu_global;             // degrees of freedom for the half-t priors for tau
  real<lower=1> nu_local;              // degrees of freedom for the half-t priors for lambdas
  real<lower=0> slab_scale;            // slab scale for the regularized horseshoe
  real<lower=0> slab_df;               // slab degrees of freedom for the regularized horseshoe
}
transformed data {
  matrix[n,d] x_std;
  vector[d] sd_x;
  
  for (j in 1:d) {
    sd_x[j] = sd(x[, j]);
    x_std[, j] = x[, j] / sd_x[j];
  }
}
parameters {
  real beta0;
  vector[d] z;
  real<lower=0> aux1_global;
  real<lower=0> aux2_global;
  vector<lower=0>[d] aux1_local;
  vector<lower=0>[d] aux2_local;
  real<lower=0> caux;
}

transformed parameters {
  real<lower=0> tau;                   // global shrinkage parameter
  vector<lower=0>[d] lambda;           // local shrinkage parameter
  vector<lower=0>[d] lambda_tilde;     // 'truncated' local shrinkage parameter
  real<lower=0> c;                     // slab scale
  vector[d] beta_std;                  // regression coefficients
  vector[n] f;                         // latent function values

  lambda = aux1_local .* sqrt(aux2_local);
  tau = aux1_global * sqrt(aux2_global) * scale_global;
  c = slab_scale * sqrt(caux);
  lambda_tilde = sqrt(c^2 * square(lambda) ./ (c^2 + tau^2 * square(lambda)));
  beta_std = z .* lambda_tilde * tau;
  f = beta0 + x_std * beta_std;
}

model {
  // half-t priors for lambdas and tau, and inverse-gamma for c^2
  z ~ normal(0, 1);
  aux1_local ~ normal(0, 1);
  aux2_local ~ inv_gamma(0.5 * nu_local, 0.5 * nu_local);
  aux1_global ~ normal(0, 1);
  aux2_global ~ inv_gamma(0.5 * nu_global, 0.5 * nu_global);
  caux ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df);
  beta0 ~ normal(0, scale_icept);
  y ~ bernoulli_logit(f);
}
generated quantities {
  vector[d] beta = inv(sd_x) .* beta_std;
}
