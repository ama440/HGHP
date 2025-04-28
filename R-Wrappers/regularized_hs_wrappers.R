library(cmdstanr)

file_hs_reg <- file.path("Stan/Comparison/regularized_horseshoe.stan")
hs_reg <- cmdstan_model(file_hs_reg)

#' Get stan data for Regularized Horshoe
#' 
#' Obtains stan data for the Regularized Horseshoe shrinkage method based on user input. This
#' is an internal function.
#' 
#' @param formula two-sided formula explaining how the outcome relates to treatment, covariates, and subgroups.
#' @param data data set
#' @param scale_icept prior std for the intercept
#' @param p_0 number of parameters anticipated to be nonzero
#' @param nu_global degrees of freedom for the half-t priors for tau
#' @param nu_local degrees of freedom for the half-t priors for lambdas
#' @param slab_scale slab scale for the regularized horseshoe
#' @param slab_df slab degrees of freedom for the regularized horseshoe
#' 
#' @return a list giving data to be passed onto `rstan` or `cmdstanr` models.
#' @examples 
get.standata.hs_reg <- function(formula, data, 
                                scale_icept = 10, p_0 = 5,
                                nu_global = 1, nu_local = 1, 
                                slab_scale = 10, slab_df = 2) {
  # Create design matrix
  X <- model.matrix(formula, data)
  
  # Obtain response vector
  formula.lhs <- as.character( formula[[2]] )
  y <- data[[formula.lhs[1]]]
  
  # Obtain scale_global
  scale_global <- p_0/ ((ncol(X) - 1 - p_0) * sqrt(nrow(X)))
  
  data_list <- list(
    'n'               = nrow(X)
    , 'd'             = ncol(X) - 1
    , 'y'             = y
    , 'x'             = X[, -1]
    , 'scale_icept'   = scale_icept
    , 'scale_global'  = scale_global
    , 'nu_global'     = nu_global
    , 'nu_local'      = nu_local
    , 'slab_scale'    = slab_scale
    , 'slab_df'       = slab_df
  )
  
  return(data_list)
}

#' Sample from posterior distribution of Regularized Horseshoe model
#' 
#' This is a wrapper function to obtain posterior samples for Regularized Horseshoe
#' shrinkage.
#' 
#' @param formula two-sided formula explaining how the outcome relates to treatment, covariates, and subgroups.
#' @param data data set
#' @param scale_icept prior std for the intercept
#' @param p_0 number of parameters anticipated to be nonzero
#' @param nu_global degrees of freedom for the half-t priors for tau
#' @param nu_local degrees of freedom for the half-t priors for lambdas
#' @param slab_scale slab scale for the regularized horseshoe
#' @param slab_df slab degrees of freedom for the regularized horseshoe
#' @param ... other arguments to pass onto `cmdstanr::sample`
#' 
#' @return an object of class `draws_df` giving posterior draws. See the `posterior` package for more details.
#' @examples 
#' 
hs_reg.mcmc <- function(formula, data, 
                        scale_icept = 10, p_0 = 5,
                        nu_global = 1, nu_local = 1, 
                        slab_scale = 10, slab_df = 2, ...) {
  # Get standata to pass to model
  standata <- get.standata.hs_reg(formula, data, 
                                  scale_icept, p_0,
                                  nu_global, nu_local, 
                                  slab_scale, slab_df)
  
  # Obtain posterior draws
  fit_hs_reg <- hs_reg$sample(data = standata, ...)
  
  # Reformat into dataframe
  fit_hs_reg <- fit_hs_reg$draws(format = 'draws_df')
  
  # Rename regression coefficients
  index <- grep('^beta\\[', colnames(fit_hs_reg))
  colnames(fit_hs_reg)[index] <- colnames(standata$x)
  
  index_std <- grep('beta_std\\[', colnames(fit_hs_reg))
  colnames(fit_hs_reg)[index_std] <- paste0(colnames(standata$x), "_std")
  
  return(fit_hs_reg)
}

rm(file_hs_reg)