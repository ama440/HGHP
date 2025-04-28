# Load BhGLM package
library(BhGLM)

# Source data
source("/proj/ibrahimlab/subgroup/Github/Data/setup.R")

# Create formula
fmla <- BINARY1 ~ ARMCD + 
  ARMCD*(AGEGR1 + SEX + RACE + ETHNIC + SPMONDDY + CMBT_Hyper + 
           CMBT_Diabe + CMBT_Cardio + CMBT_Cancer + CMBT_Obesity + CMBT_Lipide + 
           BSC_Plas + BSC_Steo + BSC_Remd + NIAIDBL)

# Define groups
grouping <- list("ARMCDA", 
                 c("AGEGR1>= 65 years", "ARMCDA:AGEGR1>= 65 years"),
                 c("SEXM", "ARMCDA:SEXM"),
                 c("RACEBLACK OR AFRICAN AMERICAN", "ARMCDA:RACEBLACK OR AFRICAN AMERICAN"),
                 c("RACEOTHER", "ARMCDA:RACEOTHER"),
                 c("ETHNICNOT HISPANIC OR LATINO", "ARMCDA:ETHNICNOT HISPANIC OR LATINO"),
                 c("SPMONDDY", "ARMCDA:SPMONDDY"),
                 c("CMBT_Hyper1", "ARMCDA:CMBT_Hyper1"),
                 c("CMBT_Diabe1", "ARMCDA:CMBT_Diabe1"),
                 c("CMBT_Cardio1", "ARMCDA:CMBT_Cardio1"),
                 c("CMBT_Cancer1", "ARMCDA:CMBT_Cancer1"),
                 c("CMBT_Obesity1", "ARMCDA:CMBT_Obesity1"),
                 c("CMBT_Lipide1", "ARMCDA:CMBT_Lipide1"),
                 c("BSC_Plas1", "ARMCDA:BSC_Plas1"),
                 c("BSC_Steo1", "ARMCDA:BSC_Steo1"),
                 c("BSC_Remd1", "ARMCDA:BSC_Remd1"),
                 c("NIAIDBL4", "ARMCDA:NIAIDBL4"))

# Fit model
fit_bhglm <- bglm(fmla, data = covid_df,
                  family = "binomial", prior = mde(), verbose = T, group = grouping)

# Obtain summary table
summary_bhglm <- cbind(summary.bh(fit_bhglm, digits = 10)[,1], 
                       log(summary.bh(fit_bhglm, digits = 10)[, c(4,5)])
)
colnames(summary_bhglm) <- c("mean", "2.5%", "97.5%")
summary_bhglm <- as.data.frame(summary_bhglm)
summary_bhglm$variable <- rownames(summary_bhglm)
rownames(summary_bhglm) <- NULL
summary_bhglm <- summary_bhglm %>% select(variable, mean, `2.5%`, `97.5%`)

# Save results
#saveRDS(summary_bhglm, "/proj/ibrahimlab/subgroup/Github/Data/Analysis/Results/res_bhglm.Rds")
