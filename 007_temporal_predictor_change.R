
PRES.POINTS = readRDS("data/PRESPOINTS.rds")

bio_vals_current = dplyr::select(PRES.POINTS, c(bio_1, bio_10, bio_13)) %>%
  dplyr::mutate(id = dplyr::row_number())

preds_2050 = terra::rast("data/predictors_2050.tif")
preds_2070 = terra::rast("data/predictors_2070.tif")
preds_2100 = terra::rast("data/predictors_2100.tif")

bio_vals_2050 <- terra::extract(preds_2050, pts) %>%
  dplyr::rename(id = ID) %>%
  dplyr::rename_with(~ paste0(.x, "_2050"), -id)

bio_vals_2070 <- terra::extract(preds_2070, pts) %>%
  dplyr::rename(id = ID) %>%
  dplyr::rename_with(~ paste0(.x, "_2070"), -id)

bio_vals_2100 <- terra::extract(preds_2100, pts) %>%
  dplyr::rename(id = ID) %>%
  dplyr::rename_with(~ paste0(.x, "_2100"), -id)

full_bio <- bio_vals_current %>%
  dplyr::left_join(bio_vals_2050, by = "id") %>%
  dplyr::left_join(bio_vals_2070, by = "id") %>%
  dplyr::left_join(bio_vals_2100, by = "id")

head(full_bio)

full_long <- full_bio %>%
  dplyr::mutate(dplyr::across(everything(), ~ .)) %>%  # ensures tibble behaviour
  tidyr::pivot_longer(
    cols = -id,
    names_to = c("variable", "year"),
    names_pattern = "(bio_\\d+)(?:_(\\d{4}))?",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    year = dplyr::na_if(year, ""),          # convert "" → NA (just in case)
    year = dplyr::coalesce(year, "Current") # NA → "current"
  )

full_long <- full_long %>%
  dplyr::mutate(
    year = factor(year, levels = c("Current", "2050", "2070", "2100"))
  )


summary_df <- full_long %>%
  group_by(variable, year) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    sd_value   = sd(value, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(summary_df, file ="figures/predictor_changes.csv", row.names = F)

pred_changes = ggplot(full_long, aes(x = year, y = value, group = id)) +
  geom_line(alpha = 0.2) +
  # SD ribbon (uncertainty envelope)
  geom_ribbon(
    data = summary_df,
    aes(
      x = year,
      ymin = mean_value - sd_value,
      ymax = mean_value + sd_value,
      group = 1
    ),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  # mean trajectory
  geom_line(
    data = mean_df,
    aes(x = year, y = mean_value, group = 1),
    colour = "darkred",
    linewidth = 1.2
  ) +
  labs(
    x = "Climate scenario",
    y = "Predictor value"
  ) +
  facet_wrap(~ variable, scales = "free_y")

ggsave("figures/predictor_change.png", plot = pred_changes, 
       width = 9, height = 3, dpi = 450)
