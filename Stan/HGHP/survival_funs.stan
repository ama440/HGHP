

//' PWE log likelihood
//' @param y               failure / censoring time
//' @param eta             linear predictor
//' @param lambda          vector of baseline hazards
//' @param log_lambda      vector of log(lambda)
//' @param breaks          (J+1)-dim vector giving intervals
//' @param j               index of the interval for which the individual failed / was censored: j \in {1, ..., J}
//' @param event_ind       integer giving whether individual died (event_ind == 1) or was censored (event_ind == 0)
//' @param lambda_d_breaks (J-1)-dim vector giving lambda[j] * (s[j] - s[j-1]), j = 1, ..., J-1
real pwe_loglik(real y, real eta, vector lambda, vector log_lambda, vector breaks, int j, int event_ind, vector cumblhaz ) {
  real loglik;
  // Initialize cumhaz to lambda[j] ( y - s[j] )
  real cumhaz = lambda[j] * ( y - breaks[j] );
  // add on (sum_{g=1}^{j-1} lambda[j] * ( s[j] - s[j-1] )
  if ( j > 1 )
    cumhaz += cumblhaz[j-1];
  // Initialize loglik = - cumhaz * exp(x'beta) = log(Survival)
  loglik = -cumhaz * exp(eta);
  // If individual experienced event, add log hazard: log(lambda[j]) + x'beta
  if ( event_ind == 1 )
    loglik += log_lambda[j] + eta;
  return(loglik);
}

//' Weibull AFT log likelihood
//' Likelihood for the Weibull AFT model
//'
//' @param logT logarithm of event time
//' @param event_ind indicator for whether someone failed (event_ind == 1) or was censored `event_ind == 0`
//' @param eta linear predictor
//' @param inv_sigma inverse of scale parameter in Weibull model
//' @param log_sigma logarithm of scale parameter in Weibull model
real weibull_aft_loglik(real logT, int event_ind, real eta, real inv_sigma, real log_sigma) {
  real y      = (logT - eta) * inv_sigma;
  real loglik = -exp(y);   // initialize to log of survival
  if ( event_ind == 1 )
    loglik += -log_sigma + y;  // add on log hazard
  return loglik;
}

//' Log-logistic AFT log likelihood
//' Likelihood for the Log-logistic AFT model
//'
//' @param logT logarithm of event time
//' @param event_ind indicator for whether someone failed (event_ind == 1) or was censored `event_ind == 0`
//' @param eta linear predictor
//' @param inv_sigma inverse of scale parameter in log-logistic model
//' @param log_sigma logarithm of scale parameter in log-logistic model
real loglogistic_aft_loglik(real logT, int event_ind, real eta, real inv_sigma, real log_sigma) {
  real y      = (logT - eta) * inv_sigma;
  real z      = log1p_exp(-y);   // = log(1 + exp(-y))
  real loglik = -y - z;  // initialize to log of survival
  if ( event_ind == 1 )
    loglik += -log_sigma - z;
  return loglik;
}

//' Log-logistic AFT log likelihood
//' Likelihood for the Log-logistic AFT model
//'
//' @param logT logarithm of event time
//' @param event_ind indicator for whether someone failed (event_ind == 1) or was censored `event_ind == 0`
//' @param eta linear predictor
//' @param inv_sigma inverse of scale parameter in log-logistic model
//' @param log_sigma logarithm of scale parameter in log-logistic model
//' @param dist_ind indicator of which distribution (1 = weibull; 2 = log-logistic)
real aft_loglik(real logT, int event_ind, real eta, real inv_sigma, real log_sigma, int dist_ind) {
  if ( dist_ind == 1 )
    return weibull_aft_loglik(logT, event_ind, eta, inv_sigma, log_sigma);
  else if ( dist_ind == 2 )
    return loglogistic_aft_loglik(logT, event_ind, eta, inv_sigma, log_sigma);
  else reject("incorrect distribution");
  return negative_infinity(); // never reached
}






