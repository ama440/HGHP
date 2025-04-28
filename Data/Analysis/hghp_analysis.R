# Load the functions for fitting the model
source("Data/libraries.R")
source("R-Wrappers/hghp_subgroup_wrappers.R")

# Load the data
df <- readRDS("Data/sample_df.Rds")

# Create formula; include treatment-covariate interactions
formula <- y ~ trt*(x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9)

# Store categorical variables
fctrs <- names(Filter(is.factor, df))

# Fit model; obtain 10,000 posterior draws after 2,000-draw burn-in
# Uses noninformative prior for treatment so that this term is forced into the model
fit_hghp <- hghp.mcmc(formula, data = df, dist = 'logistic', 
                      factors = fctrs, chains = 1, parallel_chains = 1, 
                      iter_warmup = 2000, iter_sampling = 10000,
                      factor.force.in = "trt")

# Summary of draws
summary_hghp <- summarise_draws(fit_hghp, "mean", "median", "sd", "mad",
                                ~quantile(.x, probs = c(0.025, 0.975)),
                                "rhat", "ess_bulk", "ess_tail")

# Extract unstandardized estimates and remove row corresponding to lp
index_std <- grep('_std', summary_hghp$variable)
res_hghp <- summary_hghp[-c(1,index_std),]


