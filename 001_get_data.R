#######################################################
# DOWNLOAD CURRENT CLIMATE RASTER LAYERS
#######################################################

source("setup.R")

climatedata.path = "C:/Users/Clarke v Steenderen/Documents/climate_layers_masterfile"

# Download the WORLDCLIM raster layers for current time period to your PC

# wc_current = geodata::worldclim_global(
#     var = "bio",
#     res = 2.5,      # Minute degree resolution of raster layers
#     path = here::here("./data/environmental_layers/current/"),
#     version = "2.1"
# )

# Load the WORLDCLIM rasters layers we already have downloaded 
# - We don't need to run the download code above each new R session 

pred_clim_current = terra::rast( list.files(
  here::here(paste0(climatedata.path, "/current/wc2.1_2.5m/")) ,
  full.names = TRUE,
  pattern = '.tif'
))

terra::crs(pred_clim_current) = "epsg:4326"
terra::crs(pred_clim_current, describe = T)

names(pred_clim_current) = paste0("bio_", sub(".*_", "", names(pred_clim_current)))

#############################################################################
# FUTURE CLIMATES --> RUN THIS JUST ONCE TO DONWLOAD
#############################################################################

# UNCOMMENT TO RE-DOWNLOAD OR DOWNLOAD FOR THE FIRST TIME

# # Define the models you want to use
# models_to_use = c("MIROC6", "BCC-CSM2-MR", "CanESM5")
# 
# # Function to get an averaged ensemble for a specific time window
# get_future_ensemble = function(time_window, ssp_code = "245") {
# 
#   message(paste("Processing Ensemble for:", time_window))
# 
#   model_stacks = list()
# 
#   for (m in models_to_use) {
#     # Download data into a temporary directory to avoid folder clutter
#     dat = cmip6_world(
#       model = m,
#       ssp = ssp_code,
#       time = time_window,
#       res = 2.5,
#       var = "bioc",
#       path = tempdir()
#     )
#     model_stacks[[m]] = dat
#   }
#   # Average the models (Mean of the 3 stacks)
#   ensemble = (model_stacks[[1]] + model_stacks[[2]] + model_stacks[[3]]) / 3
# 
#   # Standardize projection and names
#   terra::crs(ensemble) = "epsg:4326"
#   names(ensemble) = paste0("bio_", 1:19)
# 
#   return(ensemble)
# }
# 
# # generate the ensembles
# #pred_clim_2030 = get_future_ensemble("2021-2040")
# pred_clim_2050 = get_future_ensemble("2041-2060")
# pred_clim_2070 = get_future_ensemble("2061-2080")
# pred_clim_2100 = get_future_ensemble("2081-2100")
# 
# # save them as TIFs so you don't have to download again
# terra::writeRaster(pred_clim_2050, 
#                    here::here(paste0(climatedata.path, "/future/ensemble_2050_ssp245.tif")),
#                    overwrite = TRUE)
# terra::writeRaster(pred_clim_2070, 
#                    here::here(paste0(climatedata.path, "/future/ensemble_2070_ssp245.tif")),
#                    overwrite = TRUE)
# terra::writeRaster(pred_clim_2100, 
#                    here::here(paste0(climatedata.path, "/future/ensemble_2100_ssp245.tif")),
#                    overwrite = TRUE)

read_future_env = function(file_path) {
  # Load the raster
  r = terra::rast(file_path)
  
  # Ensure names match your training data (e.g., bio_1, bio_2...)
  # This assumes you have 19 bioclim variables in order
  names(r) = paste0("bio_", 1:19)
  
  # Set the CRS just in case it wasn't preserved
  terra::crs(r) = "epsg:4326"
  
  return(r)
}

#pred_clim_2030 = read_future_env(here::here("climate_layers/future/ensemble_2030_ssp245.tif"))
pred_clim_2050 = read_future_env(here::here(paste0(climatedata.path, "/future/ensemble_2050_ssp245.tif")))
pred_clim_2070 = read_future_env(here::here(paste0(climatedata.path, "/future/ensemble_2070_ssp245.tif")))
pred_clim_2100 = read_future_env(here::here(paste0(climatedata.path, "/future/ensemble_2100_ssp245.tif")))

##############################################################################
# KOPPEN-GEIGER DATA
##############################################################################

# read in the koppen-geiger (KG) data
kg_layer = sf::st_read(dsn = paste0(climatedata.path, "/koppen_geiger"))

# Reproject KG layer
geo_proj = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
kg_layer = sf::st_transform(kg_layer, geo_proj)

#######################################################
# DOWNLOAD SPECIES GPS DATA FROM GBIF 
#######################################################

# Below, we will download GPS data from GBIF for Diaphorina citri,
# - We can download records from GBIF, or import GPS records from a .csv 
#   file that we have stored on our PC somewhere 
# - Pick the option that works for you 

# Option #1: Download species occurrences (GPS) from GBIF
 
species_name = "Acanthococcus ironsidei"

# Get keys for the list of species
results_list = rgbif::name_backbone_checklist(name_data = species_name)
species_keys = results_list$speciesKey

# download from GBIF
gbif_download = rgbif::occ_download(
  # exclude any records with geospatial issues
  rgbif::pred("hasGeospatialIssue", FALSE),
  # keep only records with available GPS coordinates
  rgbif::pred("hasCoordinate", TRUE),
  # remove absent records
  rgbif::pred("occurrenceStatus","PRESENT"), 
  # automatically looks for unique species keys (no duplication by default)
  rgbif::pred_in("speciesKey", species_keys), 
  format = "SIMPLE_CSV",
  user = "clarke.vansteenderen",
  pwd = "roxie2@!",
  email = "vsteenderen@gmail.com"
)

rgbif::occ_download_wait(gbif_download, 
                         quiet = FALSE, 
                         status_ping = 3) 

result = rgbif::occ_download_get(key = gbif_download, 
                                 overwrite = TRUE, 
                                 path = paste0("gps/"))

# but no GPS coordinates :/
gbif_data = read_table("gps/0000054-260409193756587/0000054-260409193756587.csv")

# import a csv file containing Rosali's GPS data 
sp_gps = readr::read_csv("gps/MFC_positive_records.csv") %>%
  dplyr::select(
    species,
    lat = decimalLatitude,
    lon = decimalLongitude
  )

head(sp_gps)

# Let's just keep the columns of interest
sp_data = sp_gps %>%
  dplyr::select(
    species,
    lon,
    lat
  )

#################################
# Remove duplicate GPS data and 
# remove NAs
#################################

sp_data = sp_data %>%
  dplyr::distinct(lon, lat, .keep_all= TRUE) %>%
tidyr::drop_na(lon, lat) 

# Convert one of our environmental predictors into a raster layer 
r = raster::raster(pred_clim_current[[1]])

# Extract longitude and latitude into a dataframe
xy = sp_data %>%
  dplyr::select(
    lon,
    lat
  ) %>%
  # Coerce into a data.frame object (can't be a tibble!!!)
  as.data.frame(.)
head(xy)

# Retain only 1 GPS record per grid cell
set.seed(2012)
sp_data_thin = dismo::gridSample(
  xy = xy,     # Data.frame containing lon/lat columns only 
  r = r ,      # Environmental raster (must be a raster, not spatRast object)
  n = 1        # Number of records to keep per cell
)

nrow(sp_data)
nrow(sp_data_thin)

# split the data into training and testing sets
# `%not_in%`= Negate(`%in%`)
# testData=dplyr::slice_sample(sp_data_thin, prop=.2)
# trainData=dplyr::filter(sp_data_thin, lon %not_in% testData$lon & lat %not_in% testData$lat)

#################################
# Get world map 
#################################

world_map = rnaturalearth::ne_countries(
  scale = "medium", 
  returnclass = "sf"
) 

# Plot GPS points on world map to check our locality data is correct 
global_distr = ggplot() +
  # Add raster layer of world map 
  geom_sf(data = world_map, alpha = 0.5) +
  # Add GPS points 
  geom_point(
    data = sp_data_thin, 
    size = 1.2,
    shape = 21,
    stroke = 0.2,
    color =  "red",
    fill = "red",
    aes(
      x = lon, 
      y = lat
    )
  ) +
  # Set world map CRS 
  coord_sf(
    crs = 4326,
    ylim = c(-60, 85),
    expand = FALSE
  ) + 
  xlab("Longitude") + 
  ylab("Latitude")

global_distr

ggsave("figures/global_distribution_map.png", global_distr, dpi = 450,
       height = 4, width = 10)

# Plot GPS points on world map to check our locality data is correct 
# training_testing_map = ggplot() +
#   # Add raster layer of world map 
#   geom_sf(data = world_map, alpha = 0.5) +
#   # Add GPS points 
#   geom_point(
#     data = trainData, 
#     size = 0.6,
#     color =  "royalblue",
#     aes(
#       x = lon, 
#       y = lat
#     )
#   ) +
#   # Add GPS points 
#   geom_point(
#     data = testData, 
#     size = 0.6,
#     color =  "red",
#     aes(
#       x = lon, 
#       y = lat
#     )
#   ) +
#   # Set world map CRS 
#   coord_sf(
#     crs = 4326,
#     expand = FALSE
#   ) + 
#   xlab("Longitude") + 
#   ylab("Latitude")
# 
# training_testing_map
# 
# ggsave("training_testing_map.png", training_testing_map, dpi = 450)
