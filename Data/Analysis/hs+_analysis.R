# Load the functions for fitting the model
source("Data/libraries.R")

# Load the data
df <- readRDS("Data/sample_df.Rds")

# Create formula; include treatment-covariate interactions
formula <- y ~ trt*(x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9)

# Fit model
fit_hs_plus <- bayesreg(formula, data = df %>% mutate(y = as.factor(y)), 
                        model = "binomial", prior = "hs+", 
                        n.samples = 20000, burnin = 4000, thin = 1)

# Obtain draws and summary
draws_hs_plus <- cbind(as_draws_df(t(fit_hs_plus$beta0)), 
                       as_draws_df(t(fit_hs_plus$beta))) %>% 
  select(-c(.chain, .iteration, .draw))

res_hs_plus <- summarise_draws(draws_hs_plus, 
                               "mean", "median", "sd", "mad",
                               ~quantile(.x, probs = c(0.025, 0.975)),      
                               "rhat", "ess_bulk", "ess_tail")

res_hs_plus[1, 1] <- "(Intercept)"
