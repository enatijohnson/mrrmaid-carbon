################################################################################
# 03_ponderosa_evergreen.R
#
# Ponderosa pine (SPCD 122) allometric model + height-dependent TPA scaling.
# Requires: 01_data_preparation.R objects in environment.
#
# Replicates the cottonwood approach:
#   - Quadratic allometric model (biomass_kg ~ HT_m + I(HT_m^2))
#   - Split-weight blended TPA with monotonicity constraint
#   - Parametric smooth for GEE export
#
# Valley bottom filtering: PHYSCLCD + SLOPE < 15% + spatial constraint
################################################################################

# ==============================================================================
# 1. VALLEY BOTTOM FILTERING
# ==============================================================================

study_area_valley <- st_read("studyareazip.shp")

# Toggle for topographic position filter
USE_TOPO_FILTER <- FALSE

# Step 1: Condition table — physiographic class + slope
valley_conditions <- cond_df %>%
  filter(
    SLOPE < 15,
    PHYSCLCD %in% c(21, 22, 23, 24, 25, 29, 31, 32, 34, 39)
  )

# Step 2: Plot table — topographic position (optional)
if (USE_TOPO_FILTER) {
  valley_plots <- plot_df %>%
    filter(TOPO_POSITION_PNW %in% c(5, 6, 7, 8, 9) |
             PLT_CN %in% valley_conditions$PLT_CN)
} else {
  valley_plots <- plot_df %>%
    filter(PLT_CN %in% valley_conditions$PLT_CN)
}

# Step 3: Spatial constraint
plot_coords <- valley_plots %>%
  filter(!is.na(LON), !is.na(LAT)) %>%
  select(PLT_CN, LON, LAT) %>%
  st_as_sf(coords = c("LON", "LAT"), crs = 4269)

plots_in_valley <- st_join(plot_coords, study_area_valley, join = st_within) %>%
  filter(!is.na(names(study_area_valley)[1]))

valley_plots_final <- valley_plots %>%
  filter(PLT_CN %in% plots_in_valley$PLT_CN)

cat("Valley bottom filtering:\n")
cat("  Conditions (physiographic + slope):", nrow(valley_conditions), "\n")
cat("  Plots after spatial constraint:", nrow(valley_plots_final), "\n\n")

# ==============================================================================
# 2. EXTRACT PONDEROSA PINE DATA
# ==============================================================================

pine_spcd <- c(122)

tpa_pine <- tpa(fia_clipped, byPlot = TRUE, treeType = "live",
                treeDomain = SPCD %in% pine_spcd) %>%
  filter(TPA > 0, TPA < 500,
         PLT_CN %in% valley_plots_final$PLT_CN)

pine_trees <- tree_df %>%
  filter(SPCD %in% pine_spcd,
         PLT_CN %in% valley_plots_final$PLT_CN,
         !is.na(DRYBIO_AG), !is.na(HT), HT > 0,
         HT_m < 40, biomass_kg > 0)

# Unfiltered comparison
pine_trees_all <- tree_df %>%
  filter(SPCD %in% pine_spcd, !is.na(DRYBIO_AG), !is.na(HT),
         HT_m < 40, biomass_kg > 0)

cat("Pine data summary:\n")
cat("  All ponderosa in sagebrush biome:", nrow(pine_trees_all), "\n")
cat("  Valley bottom ponderosa:", nrow(pine_trees), "\n")
cat("  Reduction:", round((1 - nrow(pine_trees)/nrow(pine_trees_all)) * 100, 1), "%\n")
cat("  Valley bottom plots with pine:", nrow(tpa_pine), "\n")
cat("  Height range:", round(min(pine_trees$HT_m), 1), "-",
    round(max(pine_trees$HT_m), 1), "m\n\n")

if (nrow(pine_trees) < 30) {
  warning("Small pine sample (n = ", nrow(pine_trees),
          "). Consider relaxing valley bottom filters.")
}

# ==============================================================================
# 3. FIT ALLOMETRIC MODELS
# ==============================================================================

library(brms)
library(bayesplot)

pine_tree_data <- pine_trees %>%
  select(PLT_CN, TREE, HT, HT_m, biomass_kg)

# Bayesian Gamma GLM — linear (reference)
cat("Fitting Bayesian Gamma GLM (linear) for pine...\n")
pine_model_linear <- brm(
  biomass_kg ~ HT_m,
  data = pine_tree_data,
  family = Gamma(link = "log"),
  chains = 4, iter = 2000, warmup = 1000,
  cores = 8, seed = 123, silent = 2, refresh = 0
)
pine_model_linear
plot(pine_model_linear)

# Bayesian Gamma GLM — quadratic (primary)
cat("Fitting Bayesian Gamma GLM (quadratic) for pine...\n")
pine_model_quad <- brm(
  biomass_kg ~ HT_m + I(HT_m^2),
  data = pine_tree_data,
  family = Gamma(link = "log"),
  chains = 4, iter = 2000, warmup = 1000,
  cores = 8, seed = 123, silent = 2, refresh = 0
)
pine_model_quad
plot(pine_model_quad)

# Model comparison via LOO
cat("\n=== LOO Comparison ===\n")
loo_comp_pine <- loo_compare(loo(pine_model_linear), loo(pine_model_quad))
print(loo_comp_pine)

# Extract quadratic coefficients
pine_coef_quad <- fixef(pine_model_quad)
pine_b0 <- pine_coef_quad[1, "Estimate"]
pine_b1 <- pine_coef_quad[2, "Estimate"]
pine_b2 <- pine_coef_quad[3, "Estimate"]

cat("\n=== Quadratic Model (Primary) ===\n")
cat("Family: Gamma(link = 'log')\n")
cat("biomass_kg ~ HT_m + I(HT_m^2)\n")
print(pine_coef_quad)
cat("\nBayesian R²:\n")
print(bayes_R2(pine_model_quad))

cat("\n=== Linear Model (Reference) ===\n")
print(fixef(pine_model_linear))
print(bayes_R2(pine_model_linear))

# Posterior predictive check
pp_pine <- pp_check(pine_model_quad) + theme_minimal() +
  labs(title = "Posterior Predictive Check — Ponderosa Pine (Quadratic)")
print(pp_pine)

# Prediction function: exp(b0 + b1*h + b2*h^2) via Gamma log-link
predict_pine_biomass_kg <- function(height_m) {
  exp(pine_b0 + pine_b1 * height_m + pine_b2 * height_m^2)
}

# Sanity check
cat("\nPrediction check:\n")
for (h in c(5, 10, 15, 20, 25, 30)) {
  cat(sprintf("  %2dm -> %.1f kg/tree\n", h, predict_pine_biomass_kg(h)))
}


performance::mae(pine_model_quad)
plot_predictions(pine_model_quad, condition = "HT_m")
# ==============================================================================
# 4. TPA SCALING (REPLICATING COTTONWOOD APPROACH)
# ==============================================================================

# Bin edges
# Pine tends to be taller than cottonwood
PINE_BIN_EDGES <- c(0, 5, 8, 11, 14, 17, 20, 24, 28, 34, 45)
PINE_W_LOW  <- 0.5
PINE_W_HIGH <- 0.7
PINE_HT_THRESHOLD <- 20

pine_tpa_result <- run_tpa_scaling(
  tree_data = pine_tree_data,
  tpa_data = tpa_pine,
  predict_fn = predict_pine_biomass_kg,
  bin_edges = PINE_BIN_EDGES,
  w_low = PINE_W_LOW,
  w_high = PINE_W_HIGH,
  ht_threshold = PINE_HT_THRESHOLD,
  species_label = "Ponderosa Pine"
)
# Quick diagnostic
cat("PLT_CN in pine_tree_data:", class(pine_tree_data$PLT_CN), "\n")
cat("PLT_CN in tpa_pine:", class(tpa_pine$PLT_CN), "\n")
cat("Overlap:", length(intersect(pine_tree_data$PLT_CN, tpa_pine$PLT_CN)), "\n")
cat("Pine tree PLT_CNs:", n_distinct(pine_tree_data$PLT_CN), "\n")
cat("TPA pine PLT_CNs:", n_distinct(tpa_pine$PLT_CN), "\n")
# ==============================================================================
# 5. VISUALIZATIONS
# ==============================================================================

# --- 5A: Individual tree allometry ---
pred_range_pine <- tibble(HT_m = seq(1, 40, by = 0.5)) %>%
  mutate(
    Quadratic = exp(pine_b0 + pine_b1 * HT_m + pine_b2 * HT_m^2),
    Linear = exp(fixef(pine_model_linear)[1, "Estimate"] +
                   fixef(pine_model_linear)[2, "Estimate"] * HT_m)
  ) %>%
  pivot_longer(-HT_m, names_to = "Model", values_to = "biomass_kg")

p_pine_allometry <- ggplot() +
  geom_point(data = pine_tree_data, aes(x = HT_m, y = biomass_kg),
             alpha = 0.15, size = 1, color = "gray50") +
  geom_line(data = pred_range_pine, aes(x = HT_m, y = biomass_kg, color = Model),
            linewidth = 1.1) +
  scale_color_manual(values = c("Quadratic" = "darkblue", "Linear" = "steelblue")) +
  labs(title = "Ponderosa Pine: Individual Tree Allometry",
       subtitle = "Valley bottom plots only | Aboveground biomass (DRYBIO_AG)",
       x = "Tree Height (m)", y = "Aboveground Biomass (kg)",
       color = "Model") +
  theme_minimal() +
  theme(legend.position = c(0.2, 0.8))

print(p_pine_allometry)
ggsave("fig_pine_allometry.png", p_pine_allometry, width = 9, height = 6, dpi = 200)

# --- 5B: TPA scaling ---
pine_constrained <- pine_tpa_result$constrained

p_pine_tpa <- ggplot() +
  geom_point(data = pine_tpa_result$merged, aes(x = HT_m, y = TPA),
             alpha = 0.15, size = 1, color = "gray50")

if (!is.null(pine_tpa_result$param_fit)) {
  smooth_seq_pine <- tibble(mean_HT = seq(min(pine_constrained$mean_HT),
                                          max(pine_constrained$mean_HT), by = 0.5))
  smooth_seq_pine$TPA <- predict(pine_tpa_result$param_fit, newdata = smooth_seq_pine)
  
  p_pine_tpa <- p_pine_tpa +
    geom_line(data = smooth_seq_pine, aes(x = mean_HT, y = TPA),
              color = "darkorange", linewidth = 1.2)
}

p_pine_tpa <- p_pine_tpa +
  geom_line(data = pine_constrained, aes(x = mean_HT, y = TPA_constrained),
            color = "cyan3", linewidth = 0.8, linetype = "dashed") +
  geom_point(data = pine_constrained, aes(x = mean_HT, y = TPA_constrained),
             color = "cyan3", size = 3, shape = 17) +
  labs(title = "Ponderosa Pine: Height-Dependent TPA Scaling",
       subtitle = sprintf("Split-weight: w=%.1f below %dm, w=%.1f above %dm",
                          PINE_W_LOW, PINE_HT_THRESHOLD, PINE_W_HIGH, PINE_HT_THRESHOLD),
       x = "Canopy Height (m)", y = "Trees per Acre (TPA)") +
  theme_minimal()

print(p_pine_tpa)
ggsave("fig_pine_tpa.png", p_pine_tpa, width = 10, height = 6, dpi = 200)

# --- 5C: Final biomass curve ---
height_seq_pine <- seq(1, 40, by = 0.5)

if (!is.null(pine_tpa_result$param_fit)) {
  pine_final_tpa <- pmax(predict(pine_tpa_result$param_fit,
                                 newdata = data.frame(mean_HT = height_seq_pine)), 1)
} else {
  loess_pine <- loess(TPA_constrained ~ mean_HT, data = pine_constrained, span = 0.75)
  pine_final_tpa <- pmax(predict(loess_pine, newdata = data.frame(mean_HT = height_seq_pine)), 1)
}

pine_biomass_curve <- tibble(
  height_m = height_seq_pine,
  tree_kg = predict_pine_biomass_kg(height_seq_pine),
  TPA = pine_final_tpa,
  biomass_Mg_ha = calc_biomass_Mg_ha(tree_kg, TPA),
  biomass_Mg_pixel = biomass_Mg_ha * 0.01
)

p_pine_biomass <- ggplot(pine_biomass_curve, aes(x = height_m, y = biomass_Mg_ha)) +
  geom_line(color = "darkblue", linewidth = 1.3) +
  labs(title = "Ponderosa Pine: Final Biomass Prediction Curve",
       subtitle = "Quadratic allometry × smoothed split-weight TPA",
       x = "Canopy Height (m)", y = "Aboveground Biomass (Mg/ha)") +
  theme_minimal()

print(p_pine_biomass)
ggsave("fig_pine_biomass_curve.png", p_pine_biomass, width = 9, height = 6, dpi = 200)

# ==============================================================================
# 6. DIAGNOSTIC TABLE
# ==============================================================================

cat("\n=== Ponderosa Pine Final Diagnostics ===\n\n")

pine_diag <- tibble(height_m = c(3, 5, 8, 10, 15, 20, 25, 30, 35)) %>%
  mutate(
    tree_kg = round(predict_pine_biomass_kg(height_m), 1),
    TPA = round(if (!is.null(pine_tpa_result$param_fit)) {
      pmax(predict(pine_tpa_result$param_fit, newdata = data.frame(mean_HT = height_m)), 1)
    } else { NA }, 1),
    Mg_ha = round(calc_biomass_Mg_ha(tree_kg, TPA), 2),
    Mg_pixel = round(Mg_ha * 0.01, 4)
  )

print(pine_diag, n = Inf, width = Inf)

# ==============================================================================
# 7. STORE RESULTS FOR GEE EXPORT
# ==============================================================================

pine_quad_coefs <- fixef(pine_model_quad)

ponderosa_export <- list(
  class = "evergreen",
  species = "Pinus ponderosa",
  allometric_model = "Bayesian Gamma GLM (quadratic, log-link)",
  allometric_equation = "biomass_kg = exp(b0 + b1*h + b2*h^2)",
  b0 = pine_b0,
  b1 = pine_b1,
  b2 = pine_b2,
  coef_table = pine_quad_coefs,
  bayes_r2 = bayes_R2(pine_model_quad),
  n_trees = nrow(pine_tree_data),
  tpa_model = "parametric_exponential_decay",
  tpa_fit = pine_tpa_result$param_fit,
  tpa_coefs = if (!is.null(pine_tpa_result$param_fit)) coef(pine_tpa_result$param_fit) else NULL,
  w_low = PINE_W_LOW,
  w_high = PINE_W_HIGH,
  ht_threshold = PINE_HT_THRESHOLD,
  max_height = 35,
  biomass_curve = pine_biomass_curve,
  brms_model = pine_model_quad,
  valley_filter = list(
    slope_max = 15,
    physclcd = c(21, 22, 23, 24, 25, 29, 31, 32, 34, 39),
    topo_filter = USE_TOPO_FILTER
  )
)