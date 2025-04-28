# Run all the analysis methods
source("Data/Analysis/hghp_analysis.R")
source("Data/Analysis/hs_analysis.R")
source("Data/Analysis/hs+_analysis.R")
source("Data/Analysis/hs_reg_analysis.R")
source("Data/Analysis/lasso_analysis.R")
source("Data/Analysis/noninf_analysis.R")

# Assign method names
res_hghp$method    <- "HGHP"
res_hs$method      <- "Horseshoe"
res_hs_plus$method <- "Horseshoe+"
res_hs_reg$method  <- "Regularized Horseshoe"
res_blasso$method  <- "Bayesian LASSO"
res_noninf$method  <- "Noninformative"

# Select only estimate and interval bounds from tables
res_hghp    <- res_hghp    %>% select(variable, mean, `2.5%`, `97.5%`, method)
res_blasso  <- res_blasso  %>% select(variable, mean, `2.5%`, `97.5%`, method)
res_noninf  <- res_noninf  %>% select(variable, mean, `2.5%`, `97.5%`, method)
res_hs      <- res_hs      %>% select(variable, mean, `2.5%`, `97.5%`, method)
res_hs_plus <- res_hs_plus %>% select(variable, mean, `2.5%`, `97.5%`, method)
res_hs_reg  <- res_hs_reg  %>% select(variable, mean, `2.5%`, `97.5%`, method)

# Form combined results table
res <- rbind(res_hghp, res_hs, res_hs_plus, res_hs_reg, res_blasso, res_noninf)
names(res) <- c("labeltext", "mean", "lower", "upper", "method")
res <- as.data.frame(res)

# Save combined results table
# saveRDS(res, "Data/Analysis/Results/res_all.Rds")

# res <- readRDS("/proj/ibrahimlab/subgroup/Github/Data/Analysis/Results/res_all.Rds")

# Forest plot of coefficients
scale <- 1
scaletext <- 1.4
pdf("Data/Analysis/Results/forest.pdf", width = 12*scale, height = 14*scale)
res |>
  group_by(method) |>
  forestplot(fn.ci_norm = c(fpDrawCircleCI, fpDrawCircleCI, fpDrawCircleCI, fpDrawCircleCI, fpDrawCircleCI, fpDrawCircleCI),
             # title="Comparison of Parameter Estimates",
             txt_gp=fpTxtGp(label=gpar(cex = scaletext*0.9, fontfamily = "Helvetica", fontface="plain"),
                            ticks=gpar(cex = scaletext*0.9),
                            title=gpar(cex = scaletext*0.9),
                            legend=gpar(cex = scaletext*0.9),
                            summary=gpar(cex = scaletext*1),
                            legend.title = gpar(cex = scaletext*1, fontface="plain"),
                            xlab=gpar(cex = scaletext*1)),
             legend_args = fpLegend(# pos = list(x = 1.3, y = 0.5),
                                    # gp = gpar(col = "#CCCCCC", fill = "#F9F9F9"),
                                    pos = "top",
                                    title = "Method"),
             boxsize = 0.12, # We set the box size to better visualize the type
             line.margin = 0.1, # We need to add this to avoid crowding
             clip = c(-2, 2),
             xlab = "Parameter Estimate",
             xticks = seq(from = -2, to = 2, by = 1),
             lwd.ci=1, ci.vertices=TRUE, ci.vertices.height = 0.00,
             col=fpColors(box=c("#1170aa", "#5aa469", "#7851a9", "#eeca3b", "#fc7d0b", "#a3acb9"), 
                          lines=c("#1170aa", "#5aa469", "#7851a9", "#eeca3b", "#fc7d0b", "#a3acb9"), 
                          zero = "gray50")
             , mar = unit(c(.1, .5, .1, .5), "cm")
             ) |>
  fp_add_header("Variable") |> 
  fp_set_zebra_style("#EFEFEF") |>
  fp_add_lines(h_18 = gpar(col = "black", lty = 2, lwd = 2))
dev.off()