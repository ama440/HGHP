data {
  int<lower=0>                    n;            // number of observations
  int<lower=0>                    p;            // total number of covariates (incl. intercept)
  array[n] int<lower=0, upper=1>  y;            // response vector (binary responses)
  matrix[n,p]                     X;            // design matrix
}
// transformed data {
//   matrix[n,p] X_std;
//   vector[p] sd_X;
//   sd_X[1] = 1;           // Set intercept sd to 1
//   X_std[,1] = X[,1];     // Keep intercept as is
//   for (j in 2:p) {
//     sd_X[j] = sd(X[,j]);
//     X_std[,j] = X[,j] / sd(X[,j]);
//   }
// }
parameters {
  vector[p] beta;              // coefficients
}
model {
  // prior
  beta[1] ~ normal(0, 1.645);
  beta[2:p] ~ normal(0, 10);

  // likelihood
  y ~ bernoulli_logit_glm(X, 0, beta);
}





