library(cmdstanr)
library(tidyverse)
library(posterior)

hghp.logistic <- cmdstanr::cmdstan_model("Stan/HGHP/hghp_logistic.stan")
hghp.normal   <- cmdstanr::cmdstan_model("Stan/HGHP/hghp_normal.stan")
hghp.pwe      <- cmdstanr::cmdstan_model("Stan/HGHP/hghp_pwe.stan")

#' Get stan data for horseshoe shrinkage method
#' 
#' Obtains stan data for the horseshoe shrinkage method based on user input. This
#' is an internal function.
#' 
#' @param formula two-sided formula explaining how the outcome relates to treatment, covariates, and subgroups. For survival models, right-hand side must be of type returned by the `Surv` function.
#' @param data data set
#' @param dist name of the distribution. Available distributions are `"logistic"`, `"normal"`, and `"pwe"` (piecewise exponential).
#' @param factors character vector giving the names of the factors present in `formula` to perform shrinkage
#' @param factor.force.in character vector giving the main effects to force into the model. If NULL, no factors / interactions will be shrunk
#' @param beta.force.sd standard deviation for covariates forced into the model (e.g., intercept)
#' @param sigma.mean mean hyperparameter for half-normal prior on standard deviation. Ignored if `dist != "normal"`
#' @param sigma.sd  standard deviation hyperparameter for half-normal prior on standard deviation. Ignored if `dist != "normal"`
#' @param hazard.mean mean parameter for half-normal prior on baseline hazards. Ignored if `dist != "pwe"`
#' @param hazard.sd standard deviation parameter for half-normal prior on baseline hazards. Ignored if `dist != "pwe"`
#' @param breaks a (J+1)-dimensional vector giving endpoints of the intervals for the piecewise exponential model. If NULL, defaults to J = 5 endpoints using quantiles of observed (non-censored) failure times.
#' 
#' @return a list giving data to be passed onto `rstan` or `cmdstanr` models.
#' @examples 
#' ##' A normal example
#' set.seed(123)
#' trt <- rbinom(20, 1, 0.5)   ##' treatment indicator
#' x <- rbinom(20, 1, 0.5)   ##' factor
#' y <- rnorm(20, 1 + trt - x + trt*x, sd = 2)  ##' outcome
#' data <- data.frame(y = y, x = x)
#' subgroup.hs.standata(y ~ x, data, dist = 'normal', factors = 'x')
#' 
#' ##' A piecewise-exponential example
#' set.seed(123)
#' trt <- rbinom(50, 1, 0.5)   ##' treatment indicator
#' x1  <- rbinom(50, 1, 0.5)   ##' first factor
#' x2  <- rbinom(50, 1, 0.5)   ##' second factor
#' C  <- rexp(50, 1)           ##' (unobserved) censoring variable
#' t  <- rexp(50, exp(-1 + trt + trt * x1))   ##' (partially unobserved) failure time
#' y  <- pmin(C, t)   ##' observed failure time
#' event <- as.numeric( t <= C )   ##' event indicator
#' data <- data.frame(y = y, event = event, trt = trt, x1 = x1, x2 = x2)
#' subgroup.hs.standata(Surv(y, event) ~ trt*x1*x2, data, dist = 'pwe', factors = c('x1', 'x2'))
subgroup.hs.standata <- function(
  formula, data, dist
  , factors
  , factor.force.in = NULL
  , beta.force.sd = 10
  , sigma.mean = 0, sigma.sd = 10
  , hazard.mean = 0, hazard.sd = 10
  , breaks = NULL
) {
  ##' Check if formula is two-sided
  if(length(formula) != 3 | class(formula) != 'formula')
    stop('formula must be a two-sided object of class "formula"')
  
  ##' Check if intercept is included if dist is normal or logistic
  if ( dist %in% c('normal', 'logistic')) {
    if ( attr(terms(formula, data = data), 'intercept') != 1 )
      stop('formula must include an intercept')
  }
  
  ##' Make sure variables listed as factors are factors
  factor.indx <- which(names(data) %in% factors)
  for ( i in factor.indx ) {
    data[, i] <- factor(data[, i])
  }
  
  ##' Construct formula and get design matrix
  X <- model.matrix(formula, data)
  
  ##' Remove intercept term if dist == "pwe"
  if ( (dist == "pwe") & ("(Intercept)" %in% colnames(X)) ) {
    message("Note: formula contained an intercept term. Removing the intercept term for the proportional hazards model")
    X <- X[, -1]
  }
  
  ##' Obtain response variable (and event indicator if applicable)
  formula.lhs <- as.character( formula[[2]] )
  if ( length(formula.lhs) > 1 ) {
    if ( !(dist %in% c('pwe') ) )
      stop("invalid formula")
    y        <- data[[formula.lhs[2]]]
    eventind <- data[[formula.lhs[3]]]
    if ( any( y<= 0 ) )
      stop('Times must be positive')
  } else {
    if ( !( dist %in% c('logistic', 'normal')) )
      stop('invalid dist or formula')
    y <- data[[formula.lhs[1]]]
    eventind <- rep(NA, length(y))
  }
  
  ##' For the piecewise exponential model, create the number of breaks (if necessary)
  ##' and give the index of the interval into which each observation was censored or
  ##' experienced an event (called intindx)
  if ( dist == 'pwe' ) {
    if ( is.null(breaks) ) {
      J         <- 5
      events    <- y[eventind==1]
      quantiles <- quantile(events, probs = (1:(J-1)) / J)
      breaks    <- c(0, quantiles, Inf)
    }
    J <- length(breaks)-1
    intindx <- as.numeric( cut(y, breaks) )
  }
  
  ##' Get information needed to index parameters properly
  ##'    nfactors = scalar giving number of factors
  ##'    levels = nfactors-dim vector giving name of the levels for each factor
  ##'    levels_m_ref = same as levels but with reference level removed
  ##'    nlevelparms = total number of levels with reference removed
  nfactors         <- length(factors)
  levels           <- lapply(factors, function(f) levels( as.factor( data[, f] ) ) )
  names(levels)    <- factors
  levels_m_ref     <- lapply(levels, function(l) l[-1])
  nlevelparms      <- sum( sapply( levels_m_ref, function(x) length(x) ) )
  
  ##' Construct pxG matrix giving to which level (if any) each column of X pertains
  ##'    lvls.mtx[j,g] = 0 if X[,j] is not associated with any factor, j = 1, ..., p; g = 1, ..., G
  ##'    lvls.mtx[j,g] = level of group g pertaining to x otherwise, j = 1, ..., p; g = 1, ..., G
  lvls.mtx <- matrix(0, nrow = ncol(X), ncol = nfactors)
  rownames(lvls.mtx) <- colnames(X)
  colnames(lvls.mtx) <- factors
  for ( g in 1:nfactors ) {
    fctr_g <- factors[g]
    lvls_g <- levels_m_ref[[fctr_g]]
    for ( l in 1:length(lvls_g) ) {
      name_gl <- paste0(fctr_g, lvls_g[l])
      indx_gl <- grepl(name_gl, rownames(lvls.mtx))
      lvls.mtx[indx_gl, g] <- l
    }
  }
  ##' Change lvls.mtx so the indices are unique, e.g.,
  ##'    1, 0      1, 0
  ##'    0, 1  --> 0, 2
  ##'    0, 2      0, 3
  if ( nfactors > 1 ) {
    cntr <- max(lvls.mtx[, 1])
    for ( g in 2:nfactors ) {
      lvls.mtx[, g] <- (lvls.mtx[, g] > 0) * (lvls.mtx[, g] + cntr)
      cntr <- max(lvls.mtx[, g])
    }
  }
  
  ##' Get nlevelparms-dim vector that gives to which group the flattened vector of levels correspond, e.g.,
  ##'   if G = 2, L[1] = 2, L[2] = 3, then
  ##'   lvlgrp.indx = c(1, 2, 2)
  lvlgrp.indx <- unlist( lapply( seq_along(levels_m_ref), function(i) rep(i, length(levels_m_ref[[i]]) ) ) )
  
  ##' Change lvls.mtx to reflect any forced-in main effects (i.e., set to 0)
  if(!(is.null(factor.force.in))) {
    for ( g in factor.force.in ) {
      xnames_g <- paste0( g, levels_m_ref[[g]] )
      lvls.mtx[xnames_g, g] <- 0
    }
  }
  
  ##' Get indices for forced-in covariates and variable selection covariates
  indx.force  <- which(rowSums(lvls.mtx) == 0)
  indx.select <- which(rowSums(lvls.mtx) > 0)
  p.force     <- length(indx.force)
  p.select    <- length(indx.select)
  
  ##' Return stan data
  lst <- list(
      'n'             = nrow(X)
    , 'p'             = ncol(X)
    , 'G'             = nfactors
    , 'psel'          = p.select
    , 'pforce'        = p.force
    , 'indx_select'   = as.array(indx.select)
    , 'indx_force'    = as.array(indx.force)
    , 'nlevelparms'   = nlevelparms
    , 'lvl_indx'      = lvls.mtx
    , 'lvlgrp_indx'   = as.array(lvlgrp.indx)
    , 'y'             = y
    , 'X'             = X
    , 'beta_force_sd' = beta.force.sd
  )
  if ( dist == 'normal' ) {
    lst[['sigma_mean']] <- sigma.mean
    lst[['sigma_sd']]   <- sigma.sd
  }
  if ( dist == 'pwe' ) {
    if ( length(hazard.mean) != J )
      hazard.mean <- rep(hazard.mean[1], J)
    if ( length(hazard.sd) != J )
      hazard.sd <- rep(hazard.sd[1], J)
    lst[['hazard_mean']] <- hazard.mean
    lst[['hazard_sd']]   <- hazard.sd
    lst[['eventind']]    <- eventind
    lst[['intindx']]     <- intindx
    lst[['J']]           <- J
    lst[['breaks']]      <- breaks
  }
  attr(lst, 'levels_m_ref') <- levels_m_ref
  return(lst)
}



#' Sample from posterior distribution of horseshoe model
#' 
#' This is a wrapper function to obtain posterior samples for grouped horseshoe
#' shrinkage. Shrinkage is performed at
#' \itemize{
#'    \item Overall shrinkage
#'    \item The group (factor)
#'    \item Each level of the group
#'    \item Individual-covariate shrinkage
#' }
#' 
#' @param formula two-sided formula explaining how the outcome relates to treatment, covariates, and subgroups. For survival models, right-hand side must be of type returned by the `Surv` function.
#' @param data data set
#' @param dist name of the distribution. Available distributions are `"logistic"`, `"normal"`, and `"pwe"` (piecewise exponential).
#' @param factors character vector giving the names of the factors present in `formula` to perform shrinkage
#' @param factor.force.in character vector giving the main effects to force into the model. If NULL, no factors / interactions will be shrunk
#' @param beta.force.sd standard deviation for covariates forced into the model (e.g., intercept)
#' @param sigma.mean mean hyperparameter for half-normal prior on standard deviation. Ignored if `dist != "normal"`
#' @param sigma.sd  standard deviation hyperparameter for half-normal prior on standard deviation. Ignored if `dist != "normal"`
#' @param hazard.mean mean parameter for half-normal prior on baseline hazards. Ignored if `dist != "pwe"`
#' @param hazard.sd standard deviation parameter for half-normal prior on baseline hazards. Ignored if `dist != "pwe"`
#' @param breaks a (J+1)-dimensional vector giving endpoints of the intervals for the piecewise exponential model. If NULL, defaults to J = 5 endpoints using quantiles of observed (non-censored) failure times. Ignored if `dist != "pwe"`
#' @param keep.all whether to keep all variables (including intermediate ones) in the output. By default, only the final model variables are returned.
#' @param ... other arguments to pass onto `cmdstanr::sample`
#' 
#' @return an object of class `draws_df` giving posterior draws. See the `posterior` package for more details.
#' @examples 
#' ##' A normal example
#' set.seed(123)
#' trt <- rbinom(20, 1, 0.5)   ##' treatment indicator
#' x <- rbinom(20, 1, 0.5)   ##' factor
#' y <- rnorm(20, 1 + trt - x + trt*x, sd = 2)  ##' outcome
#' data <- data.frame(y = y, x = x)
#' subgroup.hs.mcmc(y ~ x, data, dist = 'normal', factors = 'x', chains = 1, parallel_chains = 1)
#' 
#' ##' A piecewise-exponential example
#' set.seed(123)
#' trt <- rbinom(50, 1, 0.5)   ##' treatment indicator
#' x1  <- rbinom(50, 1, 0.5)   ##' first factor
#' x2  <- rbinom(50, 1, 0.5)   ##' second factor
#' C  <- rexp(50, 1)           ##' (unobserved) censoring variable
#' t  <- rexp(50, exp(-1 + trt + trt * x1))   ##' (partially unobserved) failure time
#' y  <- pmin(C, t)   ##' observed failure time
#' eventind <- as.numeric(t <= C)
#' data <- data.frame(y, eventind, trt, x1, x2)
#' subgroup.hs.mcmc(Surv(y, eventind) ~ trt*x1*x2, data, dist = 'pwe', factors = c('x1', 'x2'), chains = 1, parallel_chains = 1)
hghp.mcmc <- function(
  formula, data, dist
  , factors
  , factor.force.in = NULL
  , beta.force.sd = 10
  , sigma.mean = 0, sigma.sd = 10
  , hazard.mean = 0, hazard.sd = 10
  , breaks = NULL
  , keep.all = FALSE
  , ...
)  {
  names.force <- factor.force.in
  standat <- subgroup.hs.standata(  
    formula, data, dist, factors, factor.force.in = names.force, beta.force.sd = 10
    , sigma.mean = 0, sigma.sd = 10, hazard.mean = 0, hazard.sd = 10, breaks = NULL
  )
  if ( dist == 'logistic' ) {
    fit <- hghp.logistic$sample(data = standat, ...)
  } else if ( dist == 'normal' ) {
    fit <- hghp.normal$sample(data = standat, ...)
  } else if ( dist == 'pwe' ) {
    fit <- hghp.pwe$sample(data = standat, ...)
  } else {
    stop('dist must be one of "logistic", "normal", or "pwe"')
    return(NA) ## never reached
  }
  fit <- fit$draws(format = 'draws_df')
  attr(fit, 'standata') <- standat
  
  ## keep only relevant variables;
  if ( !keep.all ) {
    ## Drop any "raw" parameters
    fit <- fit[, !grepl('_raw', colnames(fit))]
    ## Keep relevant variables
    keeps.lst <- c('^lp__', '^beta\\[', '^sigma', '_std')
    indx.mtx  <- sapply(keeps.lst, function(x) grepl(x, colnames(fit) ) )
    fit       <- fit[, apply(indx.mtx, 1, any)] 
  }
  
  ## Rename regression coefficients
  indx.beta <- grepl('^beta\\[', colnames(fit))
  colnames(fit)[indx.beta] <- colnames(standat$X)
  
  ## Rename regression coefficients for standardized X
  indx <- grepl('beta_std\\[', colnames(fit))
  colnames(fit)[indx] <- paste0(colnames(standat$X), '_std')
  
  ## Rename lambdas so that they reflect covariate names
  indx <- grepl('lambda\\[', colnames(fit))
  colnames(fit)[indx] <- paste0('lambda[', colnames(standat$X[, standat$indx_select]), ']')
  
  ## Rename deltas so that they reflect factor / level names
  indx <- grepl('delta_grp\\[', colnames(fit))
  colnames(fit)[indx] <- paste0('delta_grp[', factors, ']')
  
  indx     <- grepl('delta_lvl\\[', colnames(fit))
  levels_m_ref <- attr(standat, 'levels_m_ref')
  lvlnames <- lapply( seq_along(levels_m_ref), function(i) paste0(names(levels_m_ref)[i], levels_m_ref[[i]] ) )
  lvlnames <- unlist(lvlnames)
  colnames(fit)[indx] <- paste0('delta_lvl[', lvlnames, ']')
  
  attr(fit, 'standata') <- standat
  fit
}



