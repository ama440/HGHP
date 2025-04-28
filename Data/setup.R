# Load libraries
source("Data/libraries.R")

# Set sample size
n <- 1000

# Generate some covariates
trt <- as.factor(sample(c(0, 1), size = n, replace = TRUE))
x1 <- as.factor(sample(c(0, 1), size = n, replace = TRUE))
x2 <- as.factor(sample(c(0, 1), size = n, replace = TRUE))
x3 <- as.factor(sample(c(0, 1), size = n, replace = TRUE))
x4 <- as.factor(sample(c(0, 1), size = n, replace = TRUE))
x5 <- as.factor(sample(c(0, 1, 2), size = n, replace = TRUE))
x6 <- as.factor(sample(c(0, 1, 2), size = n, replace = TRUE))
x7 <- as.factor(sample(c(0, 1, 2), size = n, replace = TRUE))
x8 <- as.factor(sample(c(0, 1, 2), size = n, replace = TRUE))
x9 <- as.factor(sample(c(0, 1, 2), size = n, replace = TRUE))

# Combine into dataframe
df <- data.frame(trt, x1, x2, x3, x4, x5, x6, x7, x8, x9)

# Initialize y
df$y <- rep(0, n)

# Create formula and set coefficients based on sparse setting
formula <- y ~ trt*(x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9)
X <- model.matrix(formula, data = df)
beta <- rep(0, length(colnames(X)))
beta[which(colnames(X) == "trt1")] <- 0.4
beta[which(colnames(X) == "x61")] <- -0.2
beta[which(colnames(X) == "x62")] <- -0.75
beta[which(colnames(X) == "trt1:x61")] <- 0.3
beta[which(colnames(X) == "trt1:x61")] <- 0.5

# Generate binary response variable
log_odds <- X %*% beta                    # mean response
pr <- 1 / (1 + exp(-log_odds))            # pass through an inv-logit function
df$y <- rbinom(nrow(df), 1, pr)           # bernoulli response variable

# Save df
# saveRDS(df, "Data/sample_df.Rds")