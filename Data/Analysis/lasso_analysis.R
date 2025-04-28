# Load the functions for fitting the model
source("Data/libraries.R")
source("R-Wrappers/lasso_wrappers.R")

# Load the data
df <- readRDS("Data/sample_df.Rds")

# Create formula; include treatment-covariate interactions
formula <- y ~ trt*(x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9)

# Fit model
fit_blasso <- blasso.mcmc(formula, data = df, chains = 1, parallel_chains = 1, 
                          iter_warmup = 2000, iter_sampling = 10000 
                          , index_force = 2
                          )

# Obtain summary
summary_blasso <- summarise_draws(fit_blasso, "mean", "median", "sd", "mad",
                                  ~quantile(.x, probs = c(0.025, 0.975), na.rm = TRUE),      
                                  "rhat", "ess_bulk", "ess_tail")

# Remove unwanted parameters
summary_blasso <- summary_blasso %>% 
  filter(!(variable %in% c("lp__", "lasso_scale_uniform", "lasso_scale")))

# Extract unstandardized estimates
index_std <- grep('_std', summary_blasso$variable)
res_blasso <- summary_blasso[-index_std,]








