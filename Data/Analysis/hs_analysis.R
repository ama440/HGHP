# Load the functions for fitting the model
source("Data/libraries.R")

# Load the data
df <- readRDS("Data/sample_df.Rds")

# Create formula; include treatment-covariate interactions
formula <- y ~ trt*(x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9)

# Fit model
fit_hs <- brm(formula, data = df, family = bernoulli(), 
              prior = prior("horseshoe(df = 1)", class = "b"),
              chains = 1, cores = 1, warmup = 4000, iter = 22000)

# Obtain summary table
summary_hs <- summary(fit_hs)
res_hs <- as.data.frame(summary_hs$fixed)
res_hs$variable <- rownames(res_hs)
rownames(res_hs) <- NULL
res_hs <- res_hs %>% select(variable, Estimate, `l-95% CI`, `u-95% CI`)
names(res_hs) <- c("variable", "mean", "2.5%", "97.5%")
res_hs[1,1] <- "(Intercept)"
