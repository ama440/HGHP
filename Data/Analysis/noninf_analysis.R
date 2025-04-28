# Load the functions for fitting the model
source("Data/libraries.R")
source("R-Wrappers/noninf_wrappers.R")

# Load the data
df <- readRDS("Data/sample_df.Rds")

# Create formula; include treatment-covariate interactions
formula <- y ~ trt*(x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9)

# Fit model
fit_noninf <- noninf.mcmc(formula, data = df, 
                          chains = 1, parallel_chains = 1, 
                          iter_warmup = 2000, iter_sampling = 10000)

# Obtain summary
summary_noninf <- summarise_draws(fit_noninf, "mean", "median", "sd", "mad",
                                  ~quantile(.x, probs = c(0.025, 0.975), na.rm = TRUE),      
                                  "rhat", "ess_bulk", "ess_tail")

# Remove unwanted variables
res_noninf <- summary_noninf %>% filter(variable != "lp__")















