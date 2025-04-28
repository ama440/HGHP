library(cmdstanr)

file_noninf <- file.path("Stan/Comparison/noninf.stan")
noninf <- cmdstan_model(file_noninf)

#' Get stan data for noninformative prior
#' 
#' Obtains stan data for the Bayesian Lasso shrinkage method based on user input. This
#' is an internal function.
#' 
#' @param formula two-sided formula explaining how the outcome relates to treatment, covariates, and subgroups.
#' @param data data set
#' 
#' @return a list giving data to be passed onto `rstan` or `cmdstanr` models.
#' @examples 
get.standata.noninf <- function(formula, data) {
  ##' Create design matrix
  X <- model.matrix(formula, data)
  
  ##' Obtain response vector
  formula.lhs <- as.character( formula[[2]] )
  y <- data[[formula.lhs[1]]]
  
  data_list <- list(
    'n'               = nrow(X)
    , 'p'             = ncol(X)
    , 'y'             = y
    , 'X'             = X
  )
  
  return(data_list)
}

#' Sample from posterior distribution of non-informative prior model
#' 
#' This is a wrapper function to obtain posterior samples for non-informative prior model.
#' 
#' @param formula two-sided formula explaining how the outcome relates to treatment, covariates, and subgroups.
#' @param data data set
#' @param ... other arguments to pass onto `cmdstanr::sample`
#' 
#' @return an object of class `draws_df` giving posterior draws. See the `posterior` package for more details.
#' @examples 
#' 
noninf.mcmc <- function(formula, data, ...) {
  # Get standata to pass to model
  standata <- get.standata.noninf(formula, data)
  
  # Obtain posterior draws
  fit_noninf <- noninf$sample(
    data = standata, 
    # seed = 123, 
    # chains = 1, 
    # parallel_chains = 1, 
    # iter_warmup = 2000, 
    # iter_sampling = 10000,
    ...
  )
  
  # Reformat into dataframe
  fit_noninf <- fit_noninf$draws(format = 'draws_df')
  
  # Rename regression coefficients
  index <- grep('^beta\\[', colnames(fit_noninf))
  colnames(fit_noninf)[index] <- colnames(standata$X)
  
  index_std <- grep('beta_std\\[', colnames(fit_noninf))
  colnames(fit_noninf)[index_std] <- paste0(colnames(standata$X), "_std")
  
  return(fit_noninf)
}

rm(file_noninf)