################################################################################
# 02_cottonwood_deciduous.R
#
# Cottonwood (SPCD 747) allometric model + height-dependent TPA scaling.
# Requires: 01_data_preparation.R objects in environment.
#
# Allometric model: quadratic (biomass_kg ~ HT_m + I(HT_m^2))
# TPA scaling: split-weight blended + monotonicity constrained + parametric smooth
################################################################################

# ==============================================================================
# 1. EXTRACT COTTONWOOD DATA AND COMPUTE DOMINANCE
# ==============================================================================

cottonwood_spcd <- c(747)

# Cottonwood TPA by plot
tpa_cottonwood <- tpa(fia_clipped, byPlot = TRUE, treeType = "live",
                      treeDomain = SPCD %in% cottonwood_spcd) %>%
  filter(TPA > 0, TPA < 500)

selected_plots <- unique(tpa_cottonwood$PLT_CN)

# Extract cottonwood trees
cottonwood_trees <- tree_df %>%
  filter(SPCD %in% cottonwood_spcd, PLT_CN %in% selected_plots,
         !is.na(DRYBIO_AG), !is.na(HT), HT > 0, biomass_kg > 0) %>%
  mutate(HT_m = HT * 0.3048)

# Cottonwood dominance (>= 50% of plot aboveground biomass)
cw_biomass_by_plot <- cottonwood_trees %>%
  group_by(PLT_CN) %>%
  summarise(cw_biomass = sum(biomass_kg, na.rm = TRUE), .groups = "drop")

total_biomass_by_plot <- tree_df %>%
  filter(PLT_CN %in% selected_plots, !is.na(DRYBIO_AG), biomass_kg > 0) %>%
  group_by(PLT_CN) %>%
  summarise(total_biomass = sum(biomass_kg, na.rm = TRUE), .groups = "drop")

dominance <- cw_biomass_by_plot %>%
  left_join(total_biomass_by_plot, by = "PLT_CN") %>%
  mutate(cw_prop = cw_biomass / total_biomass,
         is_dominant = cw_prop >= 0.5)

cat("Plots with cottonwood:", length(selected_plots), "\n")
cat("Cottonwood-dominant plots:", sum(dominance$is_dominant), "\n")
cat("Total cottonwood trees:", nrow(cottonwood_trees), "\n\n")

# ==============================================================================
# 2. FIT ALLOMETRIC MODELS
# ==============================================================================

library(brms)
library(bayesplot)
library(performance)
library(marginaleffects)

# Prepare tree data for modeling
cottonwood_tree_data <- cottonwood_trees %>%
  select(PLT_CN, TREE, HT, HT_m, biomass_kg)

# Bayesian Gamma GLM — linear (reference)
cat("Fitting Bayesian Gamma GLM (linear)...\n")
cw_model_linear <- brm(
  biomass_kg ~ HT_m,
  data = cottonwood_tree_data,
  family = Gamma(link = "log"),
  chains = 4, iter = 2000, warmup = 1000,
  cores = 8, seed = 123, silent = 2, refresh = 0
)
plot(cw_model_linear)
# Bayesian Gamma GLM — quadratic (primary)
cat("Fitting Bayesian Gamma GLM (quadratic)...\n")
cw_model_quad <- brm(
  biomass_kg ~ HT_m + I(HT_m^2),
  data = cottonwood_tree_data,
  family = Gamma(link = "log"),
  chains = 4, iter = 2000, warmup = 1000,
  cores = 8, seed = 123, silent = 2, refresh = 0
)
plot(cw_model_quad)
# Model comparison via LOO
cat("\n=== LOO Comparison ===\n")
loo_comp <- loo_compare(loo(cw_model_linear), loo(cw_model_quad))
print(loo_comp)

# Extract quadratic coefficients
cw_coef_quad <- fixef(cw_model_quad)
cw_b0 <- cw_coef_quad[1, "Estimate"]  # Intercept
cw_b1 <- cw_coef_quad[2, "Estimate"]  # HT_m
cw_b2 <- cw_coef_quad[3, "Estimate"]  # I(HT_m^2)

cat("\n=== Quadratic Model (Primary) ===\n")
cat("Family: Gamma(link = 'log')\n")
cat("biomass_kg ~ HT_m + I(HT_m^2)\n")
print(cw_coef_quad)
cat("\nBayesian R²:\n")
print(bayes_R2(cw_model_quad))
summary(cw_model_quad)

cat("\n=== Linear Model (Reference) ===\n")
print(fixef(cw_model_linear))
print(bayes_R2(cw_model_linear))

# Posterior predictive check
pp <- pp_check(cw_model_quad) + theme_minimal() +
  labs(title = "Posterior Predictive Check — Cottonwood (Quadratic)")
print(pp)

# Prediction function: exp(b0 + b1*h + b2*h^2) via Gamma log-link
predict_cw_biomass_kg <- function(height_m) {
  exp(cw_b0 + cw_b1 * height_m + cw_b2 * height_m^2)
}

performance::mae(cw_model_quad)
plot_predictions(cw_model_quad, condition = "HT_m")
# ==============================================================================
# 3. TPA SCALING
# ==============================================================================

CW_BIN_EDGES <- c(0, 4, 7, 10, 13, 16, 19, 22, 26, 32, 45)
CW_W_LOW  <- 0.3   # blend weight below threshold
CW_W_HIGH <- 0.7   # blend weight above threshold
CW_HT_THRESHOLD <- 20  # meters

cw_tpa_result <- run_tpa_scaling(
  tree_data = cottonwood_tree_data,
  tpa_data = tpa_cottonwood,
  predict_fn = predict_cw_biomass_kg,
  bin_edges = CW_BIN_EDGES,
  w_low = CW_W_LOW,
  w_high = CW_W_HIGH,
  ht_threshold = CW_HT_THRESHOLD,
  species_label = "Cottonwood"
)

# ==============================================================================
# 4. VISUALIZATIONS
# ==============================================================================

# --- 4A: Individual tree allometry (quadratic vs linear Gamma) ---
pred_range <- tibble(HT_m = seq(1, 35, by = 0.5)) %>%
  mutate(
    Quadratic = exp(cw_b0 + cw_b1 * HT_m + cw_b2 * HT_m^2),
    Linear = exp(fixef(cw_model_linear)[1, "Estimate"] +
                   fixef(cw_model_linear)[2, "Estimate"] * HT_m)
  ) %>%
  pivot_longer(-HT_m, names_to = "Model", values_to = "biomass_kg")

p_allometry <- ggplot() +
  geom_point(data = cottonwood_tree_data, aes(x = HT_m, y = biomass_kg),
             alpha = 0.15, size = 1, color = "gray50") +
  geom_line(data = pred_range, aes(x = HT_m, y = biomass_kg, color = Model),
            linewidth = 1.1) +
  scale_color_manual(values = c("Quadratic" = "darkgreen", "Linear" = "steelblue")) +
  labs(title = "Cottonwood: Individual Tree Allometry",
       subtitle = "Aboveground biomass (DRYBIO_AG)",
       x = "Tree Height (m)", y = "Aboveground Biomass (kg)",
       color = "Model") +
  theme_minimal() +
  theme(legend.position = c(0.2, 0.8))

print(p_allometry)
ggsave("fig_cottonwood_allometry.png", p_allometry, width = 9, height = 6, dpi = 200)

# --- 4B: TPA scaling with GAM overlay ---
cw_gam <- gam(TPA ~ s(HT_m), data = cw_tpa_result$merged, family = Gamma(link = "log"))
gam_seq <- data.frame(HT_m = seq(0.5, 40, by = 0.5))
gam_p <- predict(cw_gam, newdata = gam_seq, type = "response", se.fit = TRUE)
gam_seq$fit <- gam_p$fit
gam_seq$lower <- gam_p$fit - 2 * gam_p$se.fit
gam_seq$upper <- gam_p$fit + 2 * gam_p$se.fit

constrained_pts <- cw_tpa_result$constrained

p_tpa <- ggplot() +
  geom_point(data = cw_tpa_result$merged, aes(x = HT_m, y = TPA),
             alpha = 0.15, size = 1, color = "gray50") +
  geom_ribbon(data = gam_seq, aes(x = HT_m, ymin = lower, ymax = upper),
              fill = "blue", alpha = 0.1) +
  geom_line(data = gam_seq, aes(x = HT_m, y = fit),
            color = "blue", linewidth = 0.6, alpha = 0.5)

# Add parametric smooth if available
if (!is.null(cw_tpa_result$param_fit)) {
  smooth_seq <- tibble(mean_HT = seq(min(constrained_pts$mean_HT),
                                     max(constrained_pts$mean_HT), by = 0.5))
  smooth_seq$TPA <- predict(cw_tpa_result$param_fit, newdata = smooth_seq)
  
  p_tpa <- p_tpa +
    geom_line(data = smooth_seq, aes(x = mean_HT, y = TPA),
              color = "darkorange", linewidth = 1.2)
}

p_tpa <- p_tpa +
  geom_line(data = constrained_pts, aes(x = mean_HT, y = TPA_constrained),
            color = "cyan3", linewidth = 0.8, linetype = "dashed") +
  geom_point(data = constrained_pts, aes(x = mean_HT, y = TPA_constrained),
             color = "cyan3", size = 3, shape = 17) +
  labs(title = "Cottonwood: Height-Dependent TPA Scaling",
       subtitle = sprintf("Split-weight: w=%.1f below %dm, w=%.1f above %dm; orange = parametric smooth",
                          CW_W_LOW, CW_HT_THRESHOLD, CW_W_HIGH, CW_HT_THRESHOLD),
       x = "Canopy Height (m)", y = "Trees per Acre (TPA)") +
  theme_minimal()

print(p_tpa)
ggsave("fig_cottonwood_tpa.png", p_tpa, width = 10, height = 6, dpi = 200)

# --- 4C: Final biomass curve ---
height_seq <- seq(1, 35, by = 0.5)

if (!is.null(cw_tpa_result$param_fit)) {
  cw_final_tpa <- pmax(predict(cw_tpa_result$param_fit,
                               newdata = data.frame(mean_HT = height_seq)), 1)
} else {
  loess_fit <- loess(TPA_constrained ~ mean_HT, data = constrained_pts, span = 0.75)
  cw_final_tpa <- pmax(predict(loess_fit, newdata = data.frame(mean_HT = height_seq)), 1)
}

cw_biomass_curve <- tibble(
  height_m = height_seq,
  tree_kg = predict_cw_biomass_kg(height_seq),
  TPA = cw_final_tpa,
  biomass_Mg_ha = calc_biomass_Mg_ha(tree_kg, TPA),
  biomass_Mg_pixel = biomass_Mg_ha * 0.01
)

p_biomass <- ggplot(cw_biomass_curve, aes(x = height_m, y = biomass_Mg_ha)) +
  geom_line(color = "darkgreen", linewidth = 1.3) +
  labs(title = "Cottonwood: Final Biomass Prediction Curve",
       subtitle = "Quadratic allometry × smoothed split-weight TPA",
       x = "Canopy Height (m)", y = "Aboveground Biomass (Mg/ha)") +
  theme_minimal()

print(p_biomass)
ggsave("fig_cottonwood_biomass_curve.png", p_biomass, width = 9, height = 6, dpi = 200)

# ==============================================================================
# 5. DIAGNOSTIC TABLE
# ==============================================================================

cw_diag <- tibble(height_m = c(3, 5, 8, 10, 15, 20, 25, 30)) %>%
  mutate(
    tree_kg = round(predict_cw_biomass_kg(height_m), 1),
    TPA = round(if (!is.null(cw_tpa_result$param_fit)) {
      pmax(predict(cw_tpa_result$param_fit, newdata = data.frame(mean_HT = height_m)), 1)
    } else { NA }, 1),
    Mg_ha = round(calc_biomass_Mg_ha(tree_kg, TPA), 2),
    Mg_pixel = round(Mg_ha * 0.01, 4)
  )

print(cw_diag, n = Inf, width = Inf)

# ==============================================================================
# 6. STORE RESULTS FOR GEE EXPORT
# ==============================================================================

cw_quad_coefs <- fixef(cw_model_quad)

cottonwood_export <- list(
  class = "deciduous",
  species = "Populus trichocarpa (black cottonwood)",
  allometric_model = "Bayesian Gamma GLM (quadratic, log-link)",
  allometric_equation = "biomass_kg = exp(b0 + b1*h + b2*h^2)",
  b0 = cw_b0,
  b1 = cw_b1,
  b2 = cw_b2,
  coef_table = cw_quad_coefs,
  bayes_r2 = bayes_R2(cw_model_quad),
  n_trees = nrow(cottonwood_tree_data),
  tpa_model = "parametric_exponential_decay",
  tpa_fit = cw_tpa_result$param_fit,
  tpa_coefs = if (!is.null(cw_tpa_result$param_fit)) coef(cw_tpa_result$param_fit) else NULL,
  w_low = CW_W_LOW,
  w_high = CW_W_HIGH,
  ht_threshold = CW_HT_THRESHOLD,
  max_height = 35,
  biomass_curve = cw_biomass_curve,
  brms_model = cw_model_quad
)
