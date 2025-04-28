library(cmdstanr)

file_blasso <- file.path("Stan/Comparison/bayesian_lasso.stan")
blasso <- cmdstan_model(file_blasso)

#' Get stan data for Bayesian Lasso
#' 
#' Obtains stan data for the Bayesian Lasso shrinkage method based on user input. This
#' is an internal function.
#' 
#' @param formula two-sided formula explaining how the outcome relates to treatment, covariates, and subgroups.
#' @param data data set
#' @param index_force the index of the column of the design matrix corresponding to treatment, if no shrinkage on treatment is desired
#' 
#' @return a list giving data to be passed onto `rstan` or `cmdstanr` models.
#' @examples 
get.standata.blasso <- function(formula, data, index_force = NULL) {
  ##' Create design matrix
  X <- model.matrix(formula, data)
  
  ##' Obtain response vector
  formula.lhs <- as.character( formula[[2]] )
  y <- data[[formula.lhs[1]]]
  
  if (length(index_force) > 0) {
    index_force <- unique(c(1, index_force))
    index_sel <- c(1:ncol(X))[-index_force]
  } else {
    index_force <- 1
    index_sel <- c(1:ncol(X))[-index_force]
  }
  
  data_list <- list(
    'n'               = nrow(X)
    , 'p'             = ncol(X)
    , 'y'             = y
    , 'X'             = X
    , 'force_trt'     = ifelse(length(index_force) > 0, 1, 0)
    , 'index_force'   = index_force
    , 'index_sel'     = index_sel
  )
  
  return(data_list)
}

#' Sample from posterior distribution of Bayesian Lasso model
#' 
#' This is a wrapper function to obtain posterior samples for Bayesian Lasso
#' shrinkage.
#' 
#' @param formula two-sided formula explaining how the outcome relates to treatment, covariates, and subgroups.
#' @param data data set
#' @param index_force the index of the column of the design matrix corresponding to treatment, if no shrinkage on treatment is desired
#' @param ... other arguments to pass onto `cmdstanr::sample`
#' 
#' @return an object of class `draws_df` giving posterior draws. See the `posterior` package for more details.
#' @examples 
#' 
blasso.mcmc <- function(formula, data, index_force = NULL, ...) {
  # Get standata to pass to model
  standata <- get.standata.blasso(formula, data, index_force)
  
  # Obtain posterior draws
  fit_blasso <- blasso$sample(
    data = standata, 
    # seed = 123, 
    # chains = 1, 
    # parallel_chains = 1, 
    # iter_warmup = 2000, 
    # iter_sampling = 10000,
    ...
  )
  
  # Reformat into dataframe
  fit_blasso <- fit_blasso$draws(format = 'draws_df')
  
  # Rename regression coefficients
  index <- grep('^beta\\[', colnames(fit_blasso))
  colnames(fit_blasso)[index] <- colnames(standata$X)
  
  index_std <- grep('beta_std\\[', colnames(fit_blasso))
  colnames(fit_blasso)[index_std] <- paste0(colnames(standata$X), "_std")
  
  return(fit_blasso)
}

rm(file_blasso)
