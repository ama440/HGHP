This repository contains the code to fit a regression model using the Hierarchical Grouped Horseshoe Prior (HGHP). 
This method was published in [Statistics in Medicine](https://onlinelibrary.wiley.com/doi/10.1002/sim.70246) in 2025.

The Stan code for the HGHP is found in the folder `Stan/HGHP`. Since our application
used a binary outcome, we use the hghp_logistic.stan file, but the prior is also
implemented for normal and time-to-event responses. The R file containing the 
functions that obtain the necessary Stan data and fit the Stan model are found in 
the "R-Wrappers" folder.

Comparison methods include the standard Horseshoe, Horseshoe+, Regularized Horseshoe,
Bayesian LASSO, and noninformative prior. The standard Horseshoe and Horseshoe+ priors
are implemented in existing R packages, while the Stan and R code for the Regularized
Horseshoe, Bayesian LASSO, and noninformative prior are found in the "Stan/Comparison"
and `R-Wrappers` folders.

An example workflow is contained in the folder `Data`. The file `setup.R` creates
a dataset with categorical covariates and a binary response generated using a 
logistic regression model with sparse coefficient settings. An example dataset
is already included as `sample_df.Rds`. The model is fit in the file 
`Data/Analysis/hghp_analysis.R`. Code for fitting the comparison methods is also 
included in the `Data/Analysis` folder. `Data/Analysis/comparison.R` runs
the HGHP and comparison analysis files and creates a forestplot for visualizing
the different parameter estimates. The results of these fits are contained in
the `Data/Analysis/Results` folder.