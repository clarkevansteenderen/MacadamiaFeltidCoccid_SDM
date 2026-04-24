
predictors = terra::rast("data/predictors_current.tif")
PRES.POINTS = readRDS("data/PRESPOINTS.rds")
BACK.POINTS = readRDS("data/BACKPOINTS.rds")

coords.pres = PRES.POINTS %>%
  dplyr::select(decimalLongitude, decimalLatitude)

coords.back = BACK.POINTS %>%
  dplyr::select(decimalLongitude, decimalLatitude)

all.coords = rbind(coords.pres, coords.back)

pres = PRES.POINTS %>%
  dplyr::select(-c(species, decimalLatitude, decimalLongitude))

backs = BACK.POINTS %>%
  dplyr::select(-c(species, decimalLatitude, decimalLongitude))

# bind the presence and absence data together into one data frame
full.data = dplyr::bind_rows(pres, backs)

rownames(full.data) = NULL

# Create a vector containing 0 (indicating background points) and 1 (indicating presence points)
data.vector = c(
  replicate(nrow(pres), "1"),
  replicate(nrow(backs), "0")
) 

data.vector

length(data.vector)

myBiomodData = BIOMOD_FormatingData(
  # Ensure this is a simple 0/1 numeric vector
  resp.var = as.numeric(data.vector),
  # Ensure this is a SpatRaster (from terra)
  expl.var = predictors,
  # Force to a basic data.frame
  resp.xy = as.data.frame(all.coords), 
  resp.name = "Acanthococcus ironsidei",
  PA.nb.rep = 0,
  filter.raster = TRUE
)

myBiomodData

myBiomodOptions = bm_ModelingOptions(
  data.type = 'binary',
  models = c("ANN", "CTA", "FDA", "GAM", "GBM", "GLM", "MARS", "MAXENT", "MAXNET", "RF"),
  strategy = 'default',
  bm.format = myBiomodData
)

# specify hinge
myBiomodOptions@options$MAXNET.binary.maxnet.maxnet@args.values[['_allData_allRun']]$classes = "h"
# specify rm = 1
myBiomodOptions@options$MAXNET.binary.maxnet.maxnet@args.values[['_allData_allRun']]$regmult = 1.0

print(myBiomodOptions@options$MAXNET.binary.maxnet.maxnet@args.values[['_allData_allRun']])

myBiomodModelOut = BIOMOD_Modeling(
  bm.format = myBiomodData,
  modeling.id = "Acanthococcus ironsidei",
  models = c("ANN", "CTA", "FDA", "GAM", "GBM", "GLM", "MARS", "MAXENT", "MAXNET", "RF"),
  bm.options = myBiomodOptions,
  CV.strategy = 'random',
  CV.nb.rep = 10,
  CV.perc = 0.75,
  metric.eval = c('TSS', 'ROC', 'KAPPA'),
  var.import = 3
)

myBiomodModelOut

myEval = get_evaluations(myBiomodModelOut)
head(myEval)
length(myEval)

# 2. Calculate Mean and SD using dplyr
eval_summary = myEval %>%
  group_by(algo, metric.eval) %>%
  summarise(
    Mean = mean(validation, na.rm = TRUE),
    SD   = sd(validation, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Round to 2 decimals
  mutate(across(c(Mean, SD), ~ round(.x, 2))) %>%
  # Format as "mean ± sd"
  mutate(Display = paste0(Mean, " ± ", SD))

eval_final = eval_summary %>%
  select(algo, metric.eval, Display) %>%
  pivot_wider(names_from = metric.eval, values_from = Display) %>%
  arrange(desc(ROC))

print(eval_final)

writexl::write_xlsx(eval_final, "figures/ensemble_model_metrics.xlsx")

biomod_performance = bm_PlotEvalBoxplot(
  bm.out = myBiomodModelOut,
  group.by = c('algo', 'algo'), # Group by Algorithm
  dataset = 'validation',      # Look at Validation (Testing) scores
  scales = 'free'
)

biomod_performance$tab

biomod_performance_box = biomod_performance$plot
biomod_performance_box = biomod_performance_box +
  labs(x = "Model", y = "Validation", subtitle = expression(italic("Acanthococcus ironsidei"))) +
  scale_fill_tableau()  +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.1)
  )
biomod_performance_box

# =========================================================
# 1. OBSERVED DATA
# =========================================================
obs = get_formal_data(myBiomodModelOut, "resp.var")

obs_df = data.frame(
  points = seq_along(obs),
  obs = obs
)

# =========================================================
# 2. PREDICTIONS
# =========================================================
preds = get_predictions(myBiomodModelOut)

# =========================================================
# 3. ALIGN OBSERVATIONS (CRITICAL FIX)
# =========================================================
df = preds %>%
  left_join(obs_df, by = "points")

# =========================================================
# 4. KEEP TARGET MODELS
# =========================================================
algs_keep = c(
  "ANN", "CTA", "FDA", "GAM", "GBM",
  "GLM", "MARS", "MAXENT", "MAXNET", "RF"
)

df = df %>%
  filter(algo %in% algs_keep) %>%
  mutate(pred = pred / 1000)   # only if required by biomod scaling

# =========================================================
# 5. SAFE ROC CURVE COMPUTATION
# =========================================================
roc_df = df %>%
  group_by(algo, run) %>%
  group_split() %>%
  lapply(function(sub) {
    
    keep = !is.na(sub$pred) & !is.na(sub$obs)
    
    # must have both classes
    if (length(unique(sub$obs[keep])) < 2) return(NULL)
    
    roc_obj = tryCatch(
      pROC::roc(
        response  = sub$obs[keep],
        predictor = sub$pred[keep],
        quiet = TRUE
      ),
      error = function(e) return(NULL)
    )
    
    if (is.null(roc_obj)) return(NULL)
    
    data.frame(
      fpr  = 1 - roc_obj$specificities,
      tpr  = roc_obj$sensitivities,
      algo = unique(sub$algo),
      run  = unique(sub$run)
    )
  }) %>%
  bind_rows()

# =========================================================
# 6. CLEAN ROC OUTPUT
# =========================================================
roc_df = roc_df %>%
  filter(!is.na(fpr), !is.na(tpr)) %>%
  group_by(algo, run) %>%
  filter(length(fpr) == length(tpr)) %>%
  ungroup()

# =========================================================
# 7. INTERPOLATE ONTO COMMON GRID
# =========================================================
fpr_grid = seq(0, 1, length.out = 200)

roc_interp = roc_df %>%
  group_by(algo, run) %>%
  group_split() %>%
  lapply(function(sub) {
    
    sub = sub[order(sub$fpr), ]
    
    # skip invalid curves
    if (length(unique(sub$fpr)) < 3) return(NULL)
    
    data.frame(
      algo = unique(sub$algo),
      run  = unique(sub$run),
      fpr  = fpr_grid,
      tpr  = approx(
        x = sub$fpr,
        y = sub$tpr,
        xout = fpr_grid,
        ties = mean
      )$y
    )
  }) %>%
  bind_rows()

# =========================================================
# 8. MEAN ± SD ROC CURVES
# =========================================================
roc_summary = roc_interp %>%
  group_by(algo, fpr) %>%
  summarise(
    tpr_mean = mean(tpr, na.rm = TRUE),
    tpr_sd   = sd(tpr, na.rm = TRUE),
    .groups = "drop"
  )

# =========================================================
# 9. AUC SUMMARY TABLE (CV-BASED)
# =========================================================
auc_table = df %>%
  group_by(algo, run) %>%
  summarise(
    auc = tryCatch(
      as.numeric(
        pROC::auc(
          pROC::roc(obs, pred, quiet = TRUE)
        )
      ),
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  group_by(algo) %>%
  summarise(
    auc_mean = mean(auc, na.rm = TRUE),
    auc_sd   = sd(auc, na.rm = TRUE),
    .groups = "drop"
  )

roc.plot = ggplot(roc_summary, aes(x = fpr, y = tpr_mean, colour = algo, fill = algo)) +
  
geom_ribbon(
  aes(ymin = tpr_mean - tpr_sd,
      ymax = tpr_mean + tpr_sd),
  alpha = 0.15,
  colour = NA
) +
  
geom_line(linewidth = 1) +

geom_abline(slope = 1, intercept = 0,
            linetype = "dashed",
            colour = "grey50") +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
coord_equal() +
  labs(
    x = "False Positive Rate",
    y = "True Positive Rate",
    colour = "Algorithm",
    fill = "Algorithm"
  ) +
  ggthemes::scale_colour_tableau() +
  ggthemes::scale_fill_tableau() +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right"
  )

roc.plot

ggsave("figures/biomod_performance_box.png", plot = biomod_performance_box, 
       width = 7, height = 4, dpi = 450)

ggsave("figures/ROC_plot.png", plot = roc.plot, 
       width = 5, height = 4, dpi = 450)

# Get variable importance scores
myVarImp = get_variables_importance(myBiomodModelOut)

# Plot it
var_importance = bm_PlotVarImpBoxplot(
  bm.out = myBiomodModelOut,
  group.by = c('expl.var', 'algo', 'algo') 
)

head(var_importance$tab)

var_importance$tab %>%
  dplyr::group_by(expl.var) %>%
  dplyr::summarise(mean_var_imp = mean(var.imp, na.rm = TRUE)) %>%
  dplyr::arrange(desc(mean_var_imp))

varcontrib.table = var_importance$tab %>%
  group_by(algo, expl.var) %>%
  summarise(
    var_imp = sprintf("%.2f ± %.2f", mean(var.imp), sd(var.imp)),
    .groups = "drop"
  ) %>%
  mutate(expl.var = factor(expl.var, levels = c("bio_13", "bio_1", "bio_10"))) %>%
  select(expl.var, algo, var_imp) %>%
  pivot_wider(names_from = algo, values_from = var_imp) %>%
  arrange(expl.var)

varcontrib.table

write.csv(varcontrib.table, file = "figures/var_contribs.csv", row.names = FALSE)

var_importance.plot = var_importance$plot
var_importance.plot = var_importance.plot +
  labs(x = "Model", y = "Contribution", subtitle = expression(italic("Acanthococcus ironsidei"))) +
  scale_fill_tableau() +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2)
  )

var_importance.plot

ggsave("figures/var_contrib.png", plot = var_importance.plot, 
       width = 5, height = 6, dpi = 450)

# Plot the response curves
myMessyPlot = bm_PlotResponseCurves(
  bm.out = myBiomodModelOut,
  models.chosen = 'all',
  fixed.var = 'median',
  do.plot = FALSE          # =- Don't draw the spaghetti!
)

# 2. Extract the underlying data table
resp_data = myMessyPlot$tab
resp_data$algo = sapply(strsplit(as.character(resp_data$pred.name), "_"), tail, 1)
head(resp_data)

aggregate(pred.val ~ algo, resp_data, range)

resp.plot1 = ggplot(resp_data, aes(x = expl.val, y = pred.val, color = algo, fill = algo)) +
  # geom_smooth automatically calculates the average consensus line and a confidence band
  geom_smooth(method = "loess", span = 0.5, se = FALSE) + 
  #geom_line(alpha = 0.3) +
  #stat_summary(aes(group = 1), fun = mean, geom = "line", size = 1.2) +
  scale_color_tableau() +
  scale_fill_tableau() +
  facet_wrap(~ expl.name, scales = "free_x") +
  theme_classic() +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.1)
  ) +
  labs(#title = "Model Response Curves",
       #subtitle = expression(italic("Acanthococcus ironsidei")),
       y = "Probability of Occurrence",
       x = "Environmental Variable", color = "Model", fill = "Model") +
  theme(legend.position = "right")

resp.plot1

resp.plot2 = ggplot(resp_data, aes(x = expl.val, y = pred.val, color = algo, fill = algo)) +
  # geom_smooth automatically calculates the average consensus line and a confidence band
  #geom_smooth(method = "loess", span = 0.5, se = FALSE) + 
  geom_line(alpha = 0.3) +
  stat_summary(aes(group = 1), fun = mean, geom = "line", size = 1.2) +
  scale_color_tableau() +
  scale_fill_tableau() +
  facet_wrap(~ expl.name, scales = "free_x") +
  theme_classic() +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.1)
  ) +
  labs(#title = "Model Response Curves",
       #subtitle = expression(italic("Acanthococcus ironsidei")),
       y = "Probability of Occurrence",
       x = "Environmental Variable", color = "Model", fill = "Model") +
  theme(legend.position = "none")

resp.plot2

library(patchwork)

resp.plot = resp.plot1 + resp.plot2 +
  plot_layout(nrow = 2) +
  plot_annotation(tag_levels = "A")

resp.plot

ggsave("figures/resp_plot.png", plot = resp.plot, 
       width = 10, height = 6, dpi = 450)

##############################################################################
# Create the ensemble model
##############################################################################

myBiomodEM = BIOMOD_EnsembleModeling(
  bm.mod = myBiomodModelOut,
  models.chosen = 'all',
  #em.by = 'all',
  em.by = 'all',
  
  # COMBINATION METHOD:
  # We use Weighted Mean to trust the better models (MAXNET) more
  em.algo = c('EMwmean'),
  
  # FILTER -> keep ROC > a threshold
  metric.select = 'ROC',
  metric.select.thresh = 0.75,
  
  # Outputs
  var.import = 3,
  prob.mean = TRUE,
  prob.ci = TRUE
)

# Check if it worked
myBiomodEM

myKeptModels = get_kept_models(myBiomodEM)

# Print the list --> we have 71 models
print(myKeptModels)
length(myKeptModels)



####################################################################################
# PROJECT THIS ENSEMBLE MODEL ONTO AFRICA
####################################################################################

# just South Africa
south_africa_ext = rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
  dplyr::filter(name %in% c("South Africa", "Lesotho", "eSwatini"))

# whole African continent
africa_ext = rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
  dplyr::filter(continent == "Africa")

# southern Africa
africa_ext = rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
  dplyr::filter(name %in% c(
    "South Africa", "Lesotho", "eSwatini",
    "Namibia", "Botswana", "Malawi",
    "Zimbabwe", "Mozambique", "Zambia"
  ))

africa_vect = vect(africa_ext)

years = c("current", "2050", "2070", "2100")
mycols = rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(100))
diffcols = rev(colorRampPalette(brewer.pal(9, "PuOr"))(100))

predict_list = list()

for (yr in years) {
  
  print(paste("Processing year:", yr))
  
  # Read raster
  raster_file = paste0("data/predictors_", yr, ".tif")
  r = terra::rast(raster_file)
  
  # Crop to Africa extent
  r = terra::crop(r, terra::ext(africa_vect)) # using terra::ext instead of raster::extent
  r = terra::mask(r, africa_vect)             # optional: clips out the ocean
  
  # Set CRS and Names
  crs(r) = "EPSG:4326"
  names(r) = names(predictors)
  
  # -------------------------------------------------------------
  # THE BIOMOD2 MAGIC HAPPENS HERE
  # -------------------------------------------------------------
  
  # Project the individual models onto the new climate raster
  myProj = BIOMOD_Projection(
    bm.mod = myBiomodModelOut,       # Note: Use the base models here!
    proj.name = paste0("Proj_", yr), # Gives each year a unique folder
    new.env = r,
    models.chosen = 'all',
    build.clamping.mask = FALSE
  )
  
  # Combine them into the Ensemble forecast
  myEnsProj = BIOMOD_EnsembleForecasting(
    bm.em = myBiomodEM,              # Now use the Ensemble weights
    bm.proj = myProj,                # Feed it the projections from 3a
    models.chosen = 'all'
  )
  
  # Extract the final spatial raster to save in your list
  # By default, biomod2 spits out prob.mean or prob.wmean based on your EM setup
  final_raster = get_predictions(myEnsProj)
  
  # Extract ensemble layer only
  ens_name = grep("EMwmean", names(final_raster), value = TRUE)
  
  if (length(ens_name) == 0) {
    stop("No EMwmean layer found — check ensemble model output")
  }
  
  final_raster = final_raster[[ens_name]]
  
  # Save in list
  predict_list[[yr]] = final_raster
  
  # Plotting (Optional)
  # plot(final_raster, main = paste("Ensemble Suitability:", yr), col = mycols)
}

names(predict_list$current)
predict_list$current

current_prob = predict_list$current / 1000
current_prob_masked = terra::crop(current_prob, africa_vect)
current_prob_masked = terra::mask(current_prob, africa_vect)
current_prob_masked = terra::clamp(current_prob_masked, lower = 0, upper = 1)
current_prob_masked 

current_prob_masked_SA = terra::mask(current_prob_masked, south_africa_ext)
terra::plot(current_prob_masked_SA)

# calculate a threshold for each map
pts = terra::vect(coords.pres, geom=c("decimalLongitude", "decimalLatitude"), crs="EPSG:4326")

# extract values specifically from the raster you are plotting
vals_at_presences = terra::extract(current_prob_masked, pts)

# calculate the NEW threshold from these exact values
by(resp_data$pred.val, resp_data$algo, function(x) quantile(x, 0.1, na.rm = TRUE))
THRESH_current = quantile(vals_at_presences[,2], probs = 0.1, na.rm = TRUE)
THRESH_current

#############################################################
# Macadamia farm coordinates -> 527 unique locations in SA
#############################################################

macfarms = readxl::read_excel("gps/mac_farms.xlsx") %>%
  janitor::clean_names()

macfarms$latitude = as.numeric(macfarms$latitude)
macfarms$longitude = as.numeric(macfarms$longitude)

macfarms = macfarms %>%
  distinct(longitude, latitude, .keep_all = TRUE) %>%
  dplyr::filter(longitude > 20)

str(macfarms)

macfarms_sf = st_as_sf(macfarms, coords = c("longitude", "latitude"), crs = 4326)

#############################################################################
# CURRENT PREDICTIONS
#############################################################################

provincial_borders = ne_states(country = "South Africa", returnclass = "sf")

current_farm_suitability = terra::extract(current_prob_masked, macfarms_sf)

macfarms_sf$suitability_current = current_farm_suitability[[2]]

map.curr = PLOT.MYMAP(IN.DATA = current_prob_masked, YEAR = "Current", 
                      THRESH = THRESH_current, LEGEND.POS = "none") +
  geom_sf(data = provincial_borders, fill = NA, color = "black") +
  coord_sf(
    xlim = c(10, 42),
    ylim = c(-35, -7),
    expand = FALSE
  )

map.curr


hist(current_prob_masked, 
     main="Distribution of Suitability Scores", 
     xlab="Probability", 
     col="skyblue", 
     breaks=50)

current.plot = CONT.MAP(RASTER = current_prob_masked, subtitle = "Current",
                        legend.pos = "none") +
  geom_sf(data = provincial_borders, fill = NA, color = "black") +
  coord_sf(
    xlim = c(10, 42),
    ylim = c(-35, -7),
    expand = FALSE
  )
  # coord_sf(
  #   xlim = c(16, 33),
  #   ylim = c(-35, -22),
  #   expand = FALSE
  # ) 

current.plot

current.plot.withfarmcoords = current.plot + 
  geom_point(
    data = macfarms,
    aes(x = longitude, y = latitude),
    shape = 21,
    size = 3,
    stroke = 0.6,
    fill = NA,
    colour = "black"
  )

current.plot.withfarmcoords

#############################################################################
# 2050 PREDICTIONS
#############################################################################

future_2050_prob = predict_list$`2050` / 1000
future_2050_masked = mask(future_2050_prob, africa_vect)
future_2050_masked = clamp(future_2050_masked, lower = 0, upper = 1)
future_2050_masked

future_2050_farm_suitability = terra::extract(future_2050_masked, macfarms_sf)

macfarms_sf$suitability_2050 = future_2050_farm_suitability[[2]]

map.2050 = PLOT.MYMAP(IN.DATA = future_2050_masked, YEAR = "2050", 
                      THRESH = THRESH_current, LEGEND.POS = "none") +
  geom_sf(data = provincial_borders, fill = NA, color = "black") +
  coord_sf(
    xlim = c(10, 42),
    ylim = c(-35, -7),
    expand = FALSE
  )
  # coord_sf(
  #   xlim = c(16, 33),
  #   ylim = c(-35, -22),
  #   expand = FALSE
  # ) 

map.2050

future.2050.plot = CONT.MAP(RASTER = future_2050_masked, subtitle = "2050",
                            legend.pos = "none") +
  geom_sf(data = provincial_borders, fill = NA, color = "black") +
  coord_sf(
    xlim = c(10, 42),
    ylim = c(-35, -7),
    expand = FALSE
  )
  # coord_sf(
  #   xlim = c(16, 33),
  #   ylim = c(-35, -22),
  #   expand = FALSE
  # ) 

future.2050.plot

#############################################################################
# 2070 PREDICTIONS
#############################################################################

future_2070_prob = predict_list$`2070` / 1000
future_2070_masked = mask(future_2070_prob, africa_vect)
future_2070_masked = clamp(future_2070_masked, lower = 0, upper = 1)
future_2070_masked

future_2070_farm_suitability = terra::extract(future_2070_masked, macfarms_sf)

macfarms_sf$suitability_2070 = future_2070_farm_suitability[[2]]

map.2070 = PLOT.MYMAP(IN.DATA = future_2070_masked, YEAR = "2070", 
                      THRESH = THRESH_current, LEGEND.POS = "none") +
  geom_sf(data = provincial_borders, fill = NA, color = "black") +
  coord_sf(
    xlim = c(10, 42),
    ylim = c(-35, -7),
    expand = FALSE
  )
  # coord_sf(
  #   xlim = c(16, 33),
  #   ylim = c(-35, -22),
  #   expand = FALSE
  # ) 

map.2070

future.2070.plot = CONT.MAP(RASTER = future_2070_masked, subtitle = "2070",
                            legend.pos = "none") +
  geom_sf(data = provincial_borders, fill = NA, color = "black") +
  coord_sf(
    xlim = c(10, 42),
    ylim = c(-35, -7),
    expand = FALSE
  )
  # coord_sf(
  #   xlim = c(16, 33),
  #   ylim = c(-35, -22),
  #   expand = FALSE
  # ) 

future.2070.plot

#############################################################################
# 2100 PREDICTIONS
#############################################################################

future_2100_prob = predict_list$`2100` / 1000
future_2100_masked = mask(future_2100_prob, africa_vect)
future_2100_masked = clamp(future_2100_masked, lower = 0, upper = 1)
future_2100_masked

future_2100_farm_suitability = terra::extract(future_2100_masked, macfarms_sf)

macfarms_sf$suitability_2100 = future_2100_farm_suitability[[2]]

map.2100 = PLOT.MYMAP(IN.DATA = future_2100_masked, YEAR = "2100", 
                      THRESH = THRESH_current, LEGEND.POS = "right") +
  geom_sf(data = provincial_borders, fill = NA, color = "black") +
  coord_sf(
    xlim = c(10, 42),
    ylim = c(-35, -7),
    expand = FALSE
  )
  # coord_sf(
  #   xlim = c(16, 33),
  #   ylim = c(-35, -22),
  #   expand = FALSE
  # ) 

map.2100

future.2100.plot = CONT.MAP(RASTER = future_2100_masked, subtitle = "2100",
                            legend.pos = "right") +
  geom_sf(data = provincial_borders, fill = NA, color = "black") +
  coord_sf(
    xlim = c(10, 42),
    ylim = c(-35, -7),
    expand = FALSE
  )
  # coord_sf(
  #   xlim = c(16, 33),
  #   ylim = c(-35, -22),
  #   expand = FALSE
  # ) 

future.2100.plot

continuous.prediction.plot = 
  (current.plot + 
     future.2050.plot + 
     future.2070.plot + 
     future.2100.plot) +
  plot_layout(nrow = 2) +
  plot_annotation(tag_levels = "A")

continuous.prediction.plot

ggsave("figures/continuous_prediction_maps_southernAfrica_v2.png", plot = continuous.prediction.plot, 
       width = 10, height = 10, dpi = 450)
ggsave("figures/continuous_prediction_maps_southernAfrica.svg", plot = continuous.prediction.plot, 
       width = 10, height = 10, dpi = 450)

binary.prediction.plot = 
  (map.curr + 
     map.2050 + 
     map.2070 + 
     map.2100) +
  plot_layout(nrow = 2) +
  plot_annotation(tag_levels = "A")

binary.prediction.plot

ggsave("figures/binary_prediction_maps_southernAfrica_v2.png", plot = binary.prediction.plot, 
       width = 10, height = 10, dpi = 450)
ggsave("figures/binary_prediction_maps_southernAfrica.svg", plot = binary.prediction.plot, 
       width = 10, height = 10, dpi = 450)

current.plot.withfarmcoords

ggsave("figures/current_prediction_SA_withfarms.png", plot = current.plot.withfarmcoords, 
       width = 8, height = 8, dpi = 450)

# get the change between 2100 and current
macfarms_sf$delta_rel = (macfarms_sf$suitability_2100 -
                           macfarms_sf$suitability_current) /
  macfarms_sf$suitability_current

macfarms_sf$farm_id = seq_len(nrow(macfarms_sf))

macfarms_long = macfarms_sf %>%
  mutate(
    lon = st_coordinates(geometry)[, 1],
    lat = st_coordinates(geometry)[, 2]
  ) %>%
  st_drop_geometry() %>%
  pivot_longer(
    cols = starts_with("suitability"),
    names_to = "scenario",
    values_to = "suitability"
  ) %>%
  mutate(
    scenario = dplyr::recode(scenario,
                             suitability_current = "Current",
                             suitability_2050 = "2050",
                             suitability_2070 = "2070",
                             suitability_2100 = "2100"
    ),
    scenario = factor(scenario, levels = c("Current", "2050", "2070", "2100"))
  )


macfarms_long$exceeds_thresh = ifelse(macfarms_long$suitability > THRESH_current, 1, 0)
head(macfarms_long)


macfarms_suitability_2100 = macfarms_long %>%
  dplyr::filter(scenario == "2100")

macfarms_suitability_2100
str(macfarms_suitability_2100)

# convert points
farmpoints_sf = st_as_sf(macfarms_long, coords = c("lon", "lat"), crs = 4326)

# spatial join
farmpoints_sf = st_join(farmpoints_sf, provincial_borders ["name"])
macfarms_long$province = farmpoints_sf$name
head(macfarms_long)

future.plot.withfarmcoords = future.2100.plot + 
  ggnewscale::new_scale_fill() +
  geom_point(
    data = macfarms_suitability_2100,
    aes(x = lon, y = lat, fill = factor(exceeds_thresh)),
    shape = 21,
    size = 3,
    stroke = 1,
    colour = "black"
  ) +
  scale_fill_manual(
    values = c("0" = "white",  # unsuitable
               "1" = "darkgreen"), # suitable
    labels = c("0" = "Below threshold",
               "1" = "Above threshold"),
    name = "Suitability"
  ) +
  coord_sf(
    xlim = c(27, 33),
    ylim = c(-31.6, -22.6),
    expand = FALSE
  )

future.plot.withfarmcoords

ggsave("figures/suitable_farms_by_2100.png", plot = future.plot.withfarmcoords, 
       width = 8, height = 8, dpi = 450)
ggsave("figures/suitable_farms_by_2100.svg", plot = future.plot.withfarmcoords, 
       width = 8, height = 8, dpi = 450)

future.plot.withfarmcoords.binary = map.2100 + 
  ggnewscale::new_scale_fill() +
  geom_point(
    data = macfarms_suitability_2100,
    aes(x = lon, y = lat, fill = factor(exceeds_thresh)),
    shape = 21,
    size = 3,
    stroke = 1,
    colour = "black"
  ) +
  scale_fill_manual(
    values = c("0" = "white",  # unsuitable
               "1" = "darkgreen"), # suitable
    labels = c("0" = "Below threshold",
               "1" = "Above threshold"),
    name = "Suitability"
  ) +
  coord_sf(
    xlim = c(26, 33),
    ylim = c(-32, -22),
    expand = FALSE
  )

future.plot.withfarmcoords.binary

summary_thresh = macfarms_long %>%
  group_by(scenario, province) %>%
  summarise(
    n_farms = n(),
    n_suitable = sum(exceeds_thresh, na.rm = TRUE),
    prop_suitable = mean(exceeds_thresh, na.rm = TRUE)
  )

summary_thresh

write.csv(summary_thresh, file = "figures/suitable_farms.csv", row.names = F)


prop_suitable_farms.plot = macfarms_long %>%
  dplyr::filter(!is.na(exceeds_thresh)) %>%
  ggplot(aes(x = scenario, fill = factor(exceeds_thresh))) +
  scale_fill_manual(values = c("grey90", "forestgreen")) +
  geom_bar(position = "fill", colour = "black", alpha = 0.7) +
  #scale_y_continuous(labels = scales::percent) +
  labs(x = "Climate scenario", y = "Proportion of farms that are suitable") +
  labs(fill = "Suitable \n(> threshold)") +
  theme_classic() +
  theme(legend.position = "top") 

prop_suitable_farms.plot.prov = 
  dplyr::filter(macfarms_long, province != "Eastern Cape") %>%
  dplyr::filter(!is.na(exceeds_thresh)) %>%
  ggplot(aes(x = scenario, fill = factor(exceeds_thresh))) +
  scale_fill_manual(values = c("grey90", "forestgreen")) +
  geom_bar(position = "fill", colour = "black", alpha = 0.7) +
  facet_wrap(~ province) +
  labs(
    x = "Climate scenario",
    y = "Proportion of farms that are suitable",
    fill = "Suitable \n(> threshold)"
  ) +
  theme_classic() +
  theme(legend.position = "right")

prop_suitable_farms.plot.prov

ggsave("figures/prop_suitable_farms.prov.png", plot = prop_suitable_farms.plot.prov, 
       width = 8, height = 3, dpi = 450)

ggplot(macfarms_long, aes(x = scenario, y = suitability, fill = scenario)) +
  geom_boxplot() +
  theme_classic()

farm.change.plot = ggplot(macfarms_long, aes(
  x = scenario,
  y = suitability,
  group = farm_id
)) +
  geom_line(colour = "grey70", alpha = 0.4, linewidth = 0.6) +
  
  stat_summary(
    aes(group = 1),
    fun.data = mean_cl_normal,
    geom = "ribbon",
    alpha = 0.4,
    fill = "grey5"
  ) +
  
  stat_summary(
    aes(group = 1),
    fun = mean,
    geom = "line",
    linewidth = 1,
    colour = "#1B4F72"   # deep navy blue
  ) +
  
  geom_hline(
    yintercept = THRESH_current,
    linetype = "solid",
    colour = "#C0392B",   # muted red (less harsh than pure red)
    linewidth = 0.8
  ) +
  
  theme_classic() +
  
  labs(
    x = "Climate scenario",
    y = "Predicted climatic suitability"
  ) +
  
  scale_y_continuous(
    limits = c(0.2, 1),
    breaks = seq(0.2, 1, by = 0.1)
  )

farm.change.plot

ggsave("figures/farm_change.png", plot = farm.change.plot, 
       width = 8, height = 5, dpi = 450)

prov_summary <- macfarms_long %>%
  group_by(province) %>%
  summarise(
    n_farms = n_distinct(farm_id),
    .groups = "drop"
  )

prov_summary

unique(macfarms_long$province)

farm.change.plot.prov = ggplot(
  dplyr::filter(macfarms_long, province != "Eastern Cape"), 
aes(
  x = scenario,
  y = suitability,
  group = farm_id
)) +
  
  geom_line(colour = "grey70", alpha = 0.4, linewidth = 0.6) +
  
  stat_summary(
    aes(group = 1),
    fun.data = mean_cl_normal,
    geom = "ribbon",
    alpha = 0.4,
    fill = "grey5"
  ) +
  
  stat_summary(
    aes(group = 1),
    fun = mean,
    geom = "line",
    linewidth = 1,
    colour = "#1B4F72"
  ) +
  
  geom_hline(
    yintercept = THRESH_current,
    linetype = "solid",
    colour = "#C0392B",
    linewidth = 0.8
  ) +
  
  facet_wrap(~ province) +
  
  theme_classic() +
  
  labs(
    x = "Climate scenario",
    y = "Predicted climatic suitability"
  ) +
  
  scale_y_continuous(
    limits = c(0.2, 1),
    breaks = seq(0.2, 1, by = 0.1)
  )

farm.change.plot.prov

ggsave("figures/farm_change_prov.png", plot = farm.change.plot.prov, 
       width = 9, height = 3, dpi = 450)


suitable_farms_plots = prop_suitable_farms.plot +
  farm.change.plot +
  plot_annotation(tag_levels = "A")

ggsave("figures/farm_plots.png", plot = suitable_farms_plots, 
       width = 12, height = 4, dpi = 450)


# Look at changes over time

change.raster.2100 = change.plot(future_2100_masked, current_prob_masked,
                                 MIN = -0.2, MAX = 0.6, STEP = 0.2,
                                 subtitle = "Change: 2100 - Current") +
  geom_sf(data = provincial_borders, fill = NA, color = "black") +
  coord_sf(
    xlim = c(10, 42),
    ylim = c(-35, -7),
    expand = FALSE
  )

change.raster.2070 = change.plot(future_2070_masked, current_prob_masked,
                                 MIN = -0.2, MAX = 0.6, STEP = 0.2,
                                 subtitle = "Change: 2070 - Current", legend.pos = "none") +
  geom_sf(data = provincial_borders, fill = NA, color = "black") +
  coord_sf(
    xlim = c(10, 42),
    ylim = c(-35, -7),
    expand = FALSE
  )


change.raster.2050 = change.plot(future_2050_masked, current_prob_masked,
                                 MIN = -0.2, MAX = 0.6, STEP = 0.2,
                                 subtitle = "Change: 2050 - Current", legend.pos = "none") +
  geom_sf(data = provincial_borders, fill = NA, color = "black") +
  coord_sf(
    xlim = c(10, 42),
    ylim = c(-35, -7),
    expand = FALSE
  )


change.raster.2100
change.raster.2070
change.raster.2050

change.plots = 
  change.raster.2050 +
  change.raster.2070 +
  change.raster.2100 +
  plot_annotation(tag_levels = "A") +
  plot_layout(nrow = 2) 

change.plots

ggsave("figures/change_plots_southernAfrica_v2.png", plot = change.plots, 
       width = 10, height = 10, dpi = 450)
ggsave("figures/change_plots_southernAfrica_v2.svg", plot = change.plots, 
       width = 10, height = 10, dpi = 450)
