#######################################################
# MODEL TUNING
#######################################################

bg_pts = bg_points.df %>%
  dplyr::select(
    lon = decimalLongitude , 
    lat = decimalLatitude
  )

nrow(bg_pts)

head(bg_pts)
head(sp_data_thin)

# We need a list of the feature class (fc) and regularisation multipliers (rm) to test
list_settings = list(
  fc = c("L","Q","H","LQH"), 
  rm = c(1,2,3,4,5)
)

# Run model tuning experiments 

# Set reproducible seed
set.seed(2023)

# Start the MaxEnt java application (after pushing system memory up)

options(java.parameters = "- Xmx2024m")

dismo::maxent()

# Run model tuning 

tuning_results = 
  ENMeval::ENMevaluate(
    occs = sp_data_thin,
    envs = reduced_preds_all19r2, 
    bg = bg_pts,
    tune.args = list_settings, 
    partitions = "block",
    algorithm = "maxent.jar",
    doClamp = FALSE
  )

str(tuning_results, max.level=2)

# Visualise results 

# Plot the model tuning results

ENMeval::evalplot.stats(
  e = tuning_results,              # Variable containing ENMevaluate results 
  stats = c(                       # Which metrics to plot?
    "auc.val",                     # - Make a plot for AUC
    "or.mtp",                      # - Make a plot for omission rate (minimum training presence)
    "or.10p"                       # - Make a plot for omission rate (10th percentile)
  ),   
  color = "fc",                    # Colours lines/bars by feature class (fc)
  x.var = "rm",                    # Variable to plot on x-axis
  error.bars = FALSE               # Don't show error bars 
)


# Select the optimal model settings 

# Extract the model tuning results to a data.frame 

res = ENMeval::eval.results(tuning_results)
head(res)

write.csv(res, "data/tuning_results.csv", row.names = FALSE)

# Select the model settings (RM and FC) that optimised AICc (delta AICc == 0)
best_delta_aicc = res %>% 
  dplyr::filter(delta.AICc == 0) ;best_delta_aicc
best_delta_aicc_output = best_delta_aicc %>% t(.)
write.table(best_delta_aicc_output, "data/bestAICc_LQH1_model.txt", quote = FALSE)

# select the model that optimised AUC (highest AUC value)
best_auc = res %>% 
  dplyr::filter(auc.val.avg == max(auc.val.avg)) ;best_auc
best_auc_test_output = best_auc %>% t(.)
write.table(best_auc_test_output, "data/bestAUC_H1_model.txt", quote = FALSE)

# CBI values don't work for some reason!
# select the model that optimised CBI (highest CBI value)
# best_cbi = res %>% 
#   dplyr::filter(cbi.val.avg == max(cbi.val.avg)) ;best_cbi
# best_cbi_output = best_cbi %>% t(.)
# write.table(best_cbi_output, "data/bestCBI_Q2_model.txt", quote = FALSE)

# select the model that optimised the 10% omission rate (lowest or.10p value)
best_or.10p.avg = res %>% 
  dplyr::filter(or.10p.avg == min(or.10p.avg)) ;best_or.10p.avg
best_or.10p.avg_output = best_or.10p.avg %>% t(.)
write.table(best_or.10p.avg_output, "data/bestOR10_H2_model.txt", quote = FALSE)

# default model output
default_mod_results = res %>% 
  dplyr::filter(tune.args == "fc.LQH_rm.1") ;default_mod_results
default_mod_results = default_mod_results %>% t(.)
write.table(default_mod_results, "data/best_model_default.txt", quote = FALSE)

best_models = dplyr::bind_rows(
  best_delta_aicc %>% dplyr::mutate(opt_criterion = "AICc (delta = 0)"),
  best_auc       %>% dplyr::mutate(opt_criterion = "AUC (max val AUC)"),
  best_or.10p.avg %>% dplyr::mutate(opt_criterion = "10% omission (min)")
)

best_summary = best_models %>%
  dplyr::select(
    opt_criterion,
    fc, rm, 
    auc.val.avg, auc.val.sd,
    or.10p.avg, or.10p.sd,
    AICc, delta.AICc
  )

best_summary = best_summary %>%
  dplyr::mutate(
    dplyr::across(where(is.numeric), ~ round(.x, 3))
  )

best_summary = best_summary %>%
  dplyr::mutate(
    AUC = paste0(
      format(round(auc.val.avg, 3), nsmall = 3),
      " ± ",
      format(round(auc.val.sd, 3), nsmall = 3)
    ),
    `10% omission rate` = paste0(
      format(round(or.10p.avg, 3), nsmall = 3),
      " ± ",
      format(round(or.10p.sd, 3), nsmall = 3)
    )
  ) %>%
  dplyr::select(opt_criterion, fc, rm, AUC, `10% omission rate`, AICc, delta.AICc)

best_summary 

write.csv(best_summary , "data/best_model_summary.csv", row.names = FALSE)

####################
# PCA
####################

pca_data = speciesWd.df %>%
  dplyr::select(!c(species, decimalLatitude, decimalLongitude))
head(pca_data)

pca_res = prcomp(pca_data)
summary(pca_res)

loadings = pca_res$x
loadings

PCA = factoextra::fviz_pca_var(pca_res,
                               col.var = "contrib", # Color by contributions to the PC
                               gradient.cols = c("blue", "gold", "red"),
                               repel = TRUE) + theme_classic(); PCA

# Contributions of variables to PC1
a = factoextra::fviz_contrib(pca_res, choice = "var", axes = 1)
# Contributions of variables to PC2
b = factoextra::fviz_contrib(pca_res, choice = "var", axes = 2)
pca_contribs = gridExtra::grid.arrange(a, b, ncol=2, top='Contribution of the variables to the first two PCs')
