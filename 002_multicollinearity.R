#######################################################
# MULTICOLLINEARITY
#######################################################

# We start with all 19 predictors

all_19_all = c(
  "wc2.1_2.5m_bio_1",
  "wc2.1_2.5m_bio_2",
  "wc2.1_2.5m_bio_3",
  "wc2.1_2.5m_bio_4",
  "wc2.1_2.5m_bio_5",
  "wc2.1_2.5m_bio_6",
  "wc2.1_2.5m_bio_7",
  "wc2.1_2.5m_bio_8",
  "wc2.1_2.5m_bio_9",
  "wc2.1_2.5m_bio_10",
  "wc2.1_2.5m_bio_11",
  "wc2.1_2.5m_bio_12",
  "wc2.1_2.5m_bio_13",
  "wc2.1_2.5m_bio_14",
  "wc2.1_2.5m_bio_15",
  "wc2.1_2.5m_bio_16",
  "wc2.1_2.5m_bio_17",
  "wc2.1_2.5m_bio_18",
  "wc2.1_2.5m_bio_19"
)

# Extract climate values at focal taxon GPS points
# Focal taxon points are stored in `sp_data_thin`
head(sp_data_thin)

###############################################################
# Start with all 19 predictors
###############################################################

# Extract climate data at these points 
clim_sp = terra::extract(
  x = pred_clim_current,         # SpatRast containing all 19 WORLDCLIM layers
  y = sp_data_thin               # SpatVect or data.frame containing GPS of study taxon (lon, lat)
)

clim_sp = dplyr::select(clim_sp, !ID)
head(clim_sp)

######################################################
# All 19 reduced by R-squared
# Using the correlation wheel approach
######################################################
options(encoding = "latin1")
source("correlation_wheel.R")

corwh = corrWheel(Data = clim_sp, Threshold = 0.7)# Specify the R2 value we want as our lowest limit here

###########################
# make a correlation matrix
###########################

clim_sp_corr = clim_sp
new_colnames = gsub(".*_(\\d+)$", "\\1", names(clim_sp_corr))
names(clim_sp_corr) = new_colnames

cor_mat = cor(clim_sp_corr, method = "spearman")

# get rid of the correlations of 1 (comparing the same variable to itself, to remove diagonals)
# cor_mat[cor_mat == 1 ] = 0
# 
# corrplot(cor_mat, type = "upper", order = "hclust", 
#          tl.col = "black", tl.srt = 0)

# Keep only correlations < 0.7
cor_mat[cor_mat < 0.7 ] = 0

# Make an Igraph object from this matrix:
network = graph_from_adjacency_matrix( cor_mat, weighted=T, mode="undirected", diag=F)

# Basic chart
plot(network,
     vertex.color = "yellow")

all_19_reduced_r2 = c(
  "bio_1",
  "bio_10",
  #"12", 12 and 18 do not seem to have any effect, based on the dismo plots
  "bio_13"
  #"bio_18"
)

######################################################
# All 19 reduced by VIF
######################################################

# Identify collinear variables that should be excluded (VIF > 5)
# VIF of 5 might be best
usdm::vifstep(clim_sp, th = 5)

######################################################
# All 19 reduced by R-squared, then VIF
######################################################

reduced_preds_all19r2 = terra::subset(x = pred_clim_current, 
                                      subset = all_19_reduced_r2)

terra::writeRaster(reduced_preds_all19r2, "data/predictors_current.tif")

clim_sp_reduced = terra::extract(
  x = reduced_preds_all19r2,         
  y = sp_data_thin              
)

clim_sp_reduced = dplyr::select(clim_sp_reduced, !ID)

# we'll use all_19_reduced_r2
# bio 1, 10, and 13

# write rasters for future climates

# 2050
terra::subset(x = pred_clim_2050, subset = all_19_reduced_r2) %>%
              terra::writeRaster("data/predictors_2050.tif")

# 2070
terra::subset(x = pred_clim_2070, subset = all_19_reduced_r2) %>%
  terra::writeRaster("data/predictors_2070.tif")

# 2100
terra::subset(x = pred_clim_2100, subset = all_19_reduced_r2) %>%
  terra::writeRaster("data/predictors_2100.tif")
