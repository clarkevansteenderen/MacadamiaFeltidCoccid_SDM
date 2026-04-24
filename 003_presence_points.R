##################
# PRESENCE POINTS
##################

reduced_preds_all19r2
sp_data_thin

speciesEnv = base::data.frame(
    raster::extract(reduced_preds_all19r2, cbind(sp_data_thin$lon, sp_data_thin$lat) )
)
  
speciesWd = cbind(sp_data_thin, speciesEnv)

head(speciesWd)
  
# Process data 
speciesWd = as.data.frame(speciesWd) %>%
# Keep only the lat/lon columns and key environmental variables
dplyr::select(decimalLatitude = lat,
              decimalLongitude = lon,
              all_19_reduced_r2) %>%
dplyr::mutate(species = "Acanthococcus ironsidei") %>%
dplyr::select(species, everything()) %>%
# Drop rows with NA
tidyr::drop_na()

speciesWd.df = speciesWd
  
# Reproject CRS 
coordinates(speciesWd) <-  ~ decimalLongitude + decimalLatitude
crs_wgs84 <- CRS(SRS_string = "EPSG:4326")
slot(speciesWd, "proj4string") <- crs_wgs84
  
PRES.POINTS = speciesWd.df

saveRDS(object = PRES.POINTS, file = "data/PRESPOINTS.rds")
