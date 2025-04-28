# Load the functions for fitting the model
source("Data/libraries.R")
source("R-Wrappers/regularized_hs_wrappers.R")

# Load the data
df <- readRDS("Data/sample_df.Rds")

# Create formula; include treatment-covariate interactions
formula <- y ~ trt*(x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9)

# Fit model
fit_hs_reg <- hs_reg.mcmc(formula, data = df, chains = 1, parallel_chains = 1, 
                          iter_warmup = 4000, iter_sampling = 20000)

# Limit to columns corresponding to coefficients
fit_hs_reg_coef <- fit_hs_reg[, c(2, 1182:1210)]

# Obtain summary
res_hs_reg <- summarise_draws(fit_hs_reg_coef, "mean", "median", "sd", "mad",
                              ~quantile(.x, probs = c(0.025, 0.975), na.rm = TRUE),      
                              "rhat", "ess_bulk", "ess_tail")

res_hs_reg[1, 1] <- "(Intercept)"


