#######################################################
# LOAD REQUIRED PACKAGES AND SET UP PREFERENCES
#######################################################

if (!require("pacman"))
  install.packages("pacman")

pacman::p_load(
  tidyverse,
  dismo,
  raster,
  here,
  corrplot,
  Hmisc,
  patchwork,
  ecospat,
  # kuenm,
  gridSVG,
  gridExtra,
  grid,
  ENMeval,
  spThin,
  viridis,
  viridisLite,
  mapdata,
  maptools,
  scales,
  geosphere,
  rgdal,
  ggtext,
  rJava,
  rgeos,
  sp,
  sf,
  ggspatial,
  ecospat,
  rnaturalearth,
  rnaturalearthdata,
  megaSDM,
  caret, 
  terra,
  geodata,
  usdm,
  tidyterra,
  DHARMa,
  car,
  ggeffects,
  devtools,
  caret,
  factoextra,
  ggfortify,
  tcltk2,
  igraph,
  blockCV,
  biomod2,
  curl,
  ggthemes
)

if(!require("InformationValue")){
  devtools::install_github("selva86/InformationValue")
}

if(!require("megaSDM")){
  devtools::install_github("brshipley/megaSDM", build_vignettes = TRUE)
}

if(!require("ggConvexHull")){
  devtools::install_github("cmartin/ggConvexHull")
}

# Change ggplot theme
theme_set(
  theme_classic() +
    theme(
      panel.border = element_rect(colour = "black",
                                  fill = NA),
      axis.text = element_text(colour = "black"),
      axis.title.x = element_text(margin = unit(c(2, 0, 0, 0),
                                                "mm")),
      axis.title.y = element_text(margin = unit(c(0, 4, 0, 0),
                                                "mm")),
      legend.position = "none"
    )
)

# Set the theme for the maps
theme_opts <- list(
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_rect(fill = 'white', colour = NA),
    plot.background = element_rect(),
    axis.line = element_blank(),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    axis.title.x = element_text(colour = "black"),
    axis.title.y = element_text(colour = "black"),
    plot.title = element_text(colour = "black"),
    panel.border = element_rect(fill = NA),
    legend.key = element_blank()
  )
)


###########################################################################################
# Function for plotting a thresholded map
###########################################################################################

# PLOT.MYMAP = function(IN.DATA, YEAR, MAXP = 1, THRESH, LEGEND.POS = "right"){
# 
#   bin2 = THRESH + (MAXP - THRESH) * 0.33
#   bin3 = THRESH + (MAXP - THRESH) * 0.66
# 
#   mymap = ggplot() +
#     # 1. Plot the raster (automatically handles terra objects)
#     tidyterra::geom_spatraster(data = IN.DATA) +
# 
#     # 2. Add the country borders over the top
#     tidyterra::geom_spatvector(data = africa_vect, fill = NA, color = "black", linewidth = 0.3) +
# 
#     scale_fill_stepsn(
#       # unsuitable -> low suitability -> moderate suitability -> high suitability
#       #colors = c("grey90", "#fee08b", "#fdae61", "#d73027"), # Grey to Red
#       colors = c(
#         "grey92",   # unsuitable
#         "#ffffcc",  # low (pale yellow)
#         "#fd8d3c",  # medium (strong orange)
#         "#d73027"   # high (red)
#       ),
#       #colors = viridisLite::inferno(4),
#       breaks = c(0, THRESH, bin2, bin3),
#       limits = c(0, MAXP),
#       na.value = "white",
#       name = "Suitability",
#       labels = c("Unsuitable", "Low", "Medium", "High"),
#       # This ensures the labels are centered in the color blocks
#       guide = guide_colorsteps(
#         barheight = unit(6, "cm"),
#         ticks = FALSE,
#         even.steps = TRUE, # This makes the 4 boxes equal height!
#         show.limits = FALSE
#       )
#     ) +
#     theme_classic() +
#     labs(
#       title = paste0(YEAR),
#       x = "Longitude",
#       y = "Latitude"
#     ) +
#     theme(
#       legend.position = LEGEND.POS,
#       panel.background = element_rect(fill = "white", color = NA)
#     )
#   return(mymap)
# }

PLOT.MYMAP <- function(IN.DATA, YEAR, THRESH, LEGEND.POS = "right") {
  
  MAXP <- terra::global(IN.DATA, "max", na.rm = TRUE)[1,1]
  
  bin2 <- THRESH + (MAXP - THRESH) * 0.33
  bin3 <- THRESH + (MAXP - THRESH) * 0.66
  
  IN.CLASS <- terra::classify(
    IN.DATA,
    rcl = matrix(c(
      -Inf, THRESH, 1,
      THRESH, bin2, 2,
      bin2, bin3, 3,
      bin3, Inf, 4
    ), ncol = 3, byrow = TRUE)
  )
  
  ggplot() +
    tidyterra::geom_spatraster(
      data = IN.CLASS,
      aes(fill = factor(after_stat(value)))
    ) +
    tidyterra::geom_spatvector(
      data = africa_vect,
      fill = NA,
      color = "black",
      linewidth = 0.3
    ) +
    scale_fill_manual(
      values = c(
        "1" = "#f0f0f0",
        "2" = "#ffff80",
        "3" = "#fd8d3c",
        "4" = "#800026"
      ),
      labels = c(
        "Unsuitable",
        "Low",
        "Medium",
        "High"
      ),
      name = "Suitability",
      na.value = "white",
      na.translate = FALSE,
      drop = FALSE
    ) +
    coord_sf(expand = FALSE) +
    theme_classic() +
    labs(
      title = YEAR,
      x = "Longitude",
      y = "Latitude"
    ) +
    theme(
      legend.position = LEGEND.POS,
      panel.background = element_rect(fill = "white", color = NA)
    ) +
    # Add scale bar to bottom-right of map
    ggspatial::annotation_scale(
      location = "br",          # 'bl' = bottom left
      style = "ticks",
      width_hint = 0.2
    ) +
    # Add north arrow
    ggspatial::annotation_north_arrow(
      location = "br",
      which_north = "true",
      pad_x = unit(0.175, "in"),
      pad_y = unit(0.3, "in"),
      style = ggspatial::north_arrow_fancy_orienteering
    )
}

###########################################################################################
# Function for plotting a continuous prediction map
###########################################################################################

CONT.MAP = function(RASTER, subtitle = "", legend.pos = "right"){
  
  current.plot = ggplot() +
    # Plot MaxEnt prediction raster
    tidyterra::geom_spatraster(
      data = RASTER,
      maxcell = 5e+7         # maxcell = Inf
    ) +
    # Control raster colour and legend
    tidyterra::scale_fill_whitebox_c(
      palette = "muted",
      breaks = seq(0, 1, 0.2),
      limits = c(0, 1)
    ) +
    # Plot Africa boundary
    geom_sf(data = africa_ext, fill = NA, color = "black", size = 0.1) +
    # Control axis and legend labels 
    labs(
      subtitle = subtitle,
      x = "Longitude",
      y = "Latitude",
      fill = "P(suitability)"
    ) +
    # geom_point(
    #   data = coords.pres,
    #   size = 1,
    #   aes(x = decimalLongitude, y = decimalLatitude)
    # )  +
    # Crops map to just the geographic extent of Africa
    # coord_sf(
    #   xlim = c(16, 33),
    #   ylim = c(-35, -22),
    #   expand = FALSE
    # ) +
    # Create title for the legend
    theme(legend.position = legend.pos) +
    # Add scale bar to bottom-right of map
    ggspatial::annotation_scale(
      location = "br",          # 'bl' = bottom left
      style = "ticks",
      width_hint = 0.2
    ) +
    # Add north arrow
    ggspatial::annotation_north_arrow(
      location = "br",
      which_north = "true",
      pad_x = unit(0.175, "in"),
      pad_y = unit(0.3, "in"),
      style = ggspatial::north_arrow_fancy_orienteering
    ) +
    # Change appearance of the legend
    guides(
      fill = guide_colorbar(ticks = FALSE)
    )
  
}


############################################################################################
# Plot to visualise suitability change over time, e.g. between 2100 and the current climate
############################################################################################

change.plot = function(FUTURE.RASTER, CURRENT.RASTER, MIN, MAX, STEP, subtitle,
                       legend.pos = "right"){
  
  change.raster = FUTURE.RASTER - CURRENT.RASTER
  
  plot = ggplot() +
    # Plot MaxEnt prediction raster
    tidyterra::geom_spatraster(
      data = change.raster,
      maxcell = 5e+7         # maxcell = Inf
    ) +
    # Control raster colour and legend
    #tidyterra::scale_fill_whitebox_c(
    tidyterra::scale_fill_grass_c(
      #palette = "gn_yl",
      palette = "magma",
      breaks = seq(MIN, MAX, STEP),
      #limits = c(0, 1)
    ) +
    # Plot Africa boundary
    geom_sf(data = africa_ext, fill = NA, color = "black", size = 0.1) +
    # Control axis and legend labels 
    labs(
      subtitle = subtitle,
      x = "Longitude",
      y = "Latitude",
      fill = "Change P(suitability)"
    ) +
    # Crops map to just the geographic extent of Africa
    # coord_sf(
    #   xlim = c(16, 33),
    #   ylim = c(-35, -22),
    #   expand = FALSE
    # ) +
    # Create title for the legend
    theme(legend.position = legend.pos) +
    # Add scale bar to bottom-right of map
    ggspatial::annotation_scale(
      location = "br",          # 'bl' = bottom left
      style = "ticks",
      width_hint = 0.2
    ) +
    # Add north arrow
    ggspatial::annotation_north_arrow(
      location = "br",
      which_north = "true",
      pad_x = unit(0.175, "in"),
      pad_y = unit(0.3, "in"),
      style = ggspatial::north_arrow_fancy_orienteering
    ) +
    # Change appearance of the legend
    guides(
      fill = guide_colorbar(ticks = FALSE)
    )
  
  return(plot)
  
}
