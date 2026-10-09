################################################################################
# 05_gee_export.R
#
# Compile all vegetation class models, generate GEE-ready exports, and create
# publication-quality comparison figures.
#
# Requires: all *_export objects from scripts 02, 03, 04 in environment.
################################################################################

library(jsonlite)

cat("\n########## GEE EXPORT AND FINAL FIGURES ##########\n\n")

# ==============================================================================
# 1. PRINT FINAL EQUATIONS FOR ALL CLASSES
# ==============================================================================

cat("================================================================\n")
cat("  FINAL BIOMASS EQUATIONS FOR GEE IMPLEMENTATION\n")
cat("================================================================\n\n")

# --- Cottonwood ---
cat("1. DECIDUOUS (Cottonwood) — Bayesian Gamma GLM (quadratic, log-link) + TPA scaling\n")
cat(sprintf("   tree_kg = exp(%.6f + %.6f * h + %.8f * h²)\n",
            cottonwood_export$b0, cottonwood_export$b1, cottonwood_export$b2))
if (!is.null(cottonwood_export$tpa_coefs)) {
  cat(sprintf("   TPA(h) = %.2f + (%.2f - %.2f) * exp(-%.4f * h)\n",
              cottonwood_export$tpa_coefs["tpa_floor"],
              cottonwood_export$tpa_coefs["tpa_max"],
              cottonwood_export$tpa_coefs["tpa_floor"],
              cottonwood_export$tpa_coefs["k"]))
}
cat(sprintf("   Bayesian R² = %.4f | n = %d trees\n",
            cottonwood_export$bayes_r2[1, "Estimate"], cottonwood_export$n_trees))
cat(sprintf("   TPA split-weight: w=%.1f below %dm, w=%.1f above %dm\n\n",
            cottonwood_export$w_low, cottonwood_export$ht_threshold,
            cottonwood_export$w_high, cottonwood_export$ht_threshold))

# --- Ponderosa ---
cat("2. EVERGREEN (Ponderosa Pine) — Bayesian Gamma GLM (quadratic, log-link) + TPA scaling\n")
cat(sprintf("   tree_kg = exp(%.6f + %.6f * h + %.8f * h²)\n",
            ponderosa_export$b0, ponderosa_export$b1, ponderosa_export$b2))
if (!is.null(ponderosa_export$tpa_coefs)) {
  cat(sprintf("   TPA(h) = %.2f + (%.2f - %.2f) * exp(-%.4f * h)\n",
              ponderosa_export$tpa_coefs["tpa_floor"],
              ponderosa_export$tpa_coefs["tpa_max"],
              ponderosa_export$tpa_coefs["tpa_floor"],
              ponderosa_export$tpa_coefs["k"]))
}
cat(sprintf("   Bayesian R² = %.4f | n = %d trees\n",
            ponderosa_export$bayes_r2[1, "Estimate"], ponderosa_export$n_trees))
cat(sprintf("   Valley bottom filtered: slope < 15%%, PHYSCLCD valley classes\n\n"))

# --- Willow ---
cat("3. SHRUB (Willow) — Chronosequence\n")
cat("   age = 2 + 17 * ((h - 0.5) / 9.5)\n")
cat("   Biomass(Mg/ha) = 13.40 × ln(age) - 7.020\n")
cat(sprintf("   Height range: %.1f - %.1f m\n\n",
            willow_export$min_height, willow_export$max_height))

# --- Herbaceous ---
cat("4. HERBACEOUS — Constant\n")
cat(sprintf("   Biomass = %.1f Mg/ha (%.3f Mg/pixel)\n",
            herbaceous_export$biomass_Mg_ha, herbaceous_export$biomass_Mg_pixel))
cat(sprintf("   Source: %s\n\n", herbaceous_export$citation))

# ==============================================================================
# 2. GENERATE GEE LOOKUP TABLES
# ==============================================================================

cat("=== Generating GEE lookup tables ===\n\n")

# Deciduous: use stored biomass curve
write.csv(
  cottonwood_export$biomass_curve %>%
    select(height_m, biomass_Mg_ha, biomass_Mg_pixel),
  "gee_lookup_deciduous.csv", row.names = FALSE
)

# Evergreen: use stored biomass curve
write.csv(
  ponderosa_export$biomass_curve %>%
    select(height_m, biomass_Mg_ha, biomass_Mg_pixel),
  "gee_lookup_evergreen.csv", row.names = FALSE
)

# Shrub: generate from chronosequence function
willow_lookup <- tibble(
  height_m = seq(0, 15, by = 0.1)
) %>%
  mutate(
    biomass_Mg_ha = sapply(height_m, willow_export$biomass_fn),
    biomass_Mg_pixel = biomass_Mg_ha * 0.01
  )
write.csv(willow_lookup, "gee_lookup_shrub.csv", row.names = FALSE)

cat("CSV lookup tables saved:\n")
cat("  gee_lookup_deciduous.csv\n")
cat("  gee_lookup_evergreen.csv\n")
cat("  gee_lookup_shrub.csv\n\n")

# ==============================================================================
# 3. EXPORT GEE JAVASCRIPT MODULE
# ==============================================================================

# Generate a JavaScript module with all coefficients for direct GEE use.

gee_js <- sprintf('
// Auto-generated GEE biomass coefficients
// Source: FIA-based allometric models with height-dependent TPA scaling

exports.DECIDUOUS = {
  // Quadratic allometric model: tree_kg = max(0, b0 + b1*h + b2*h^2)
  b0: %.6f,
  b1: %.6f,
  b2: %.6f,
  // TPA scaling: TPA(h) = floor + (max - floor) * exp(-k * h)
  tpa_floor: %.4f,
  tpa_max: %.4f,
  tpa_k: %.6f,
  max_height: %d,
  species: "Populus trichocarpa"
};

exports.EVERGREEN = {
  b0: %.6f,
  b1: %.6f,
  b2: %.6f,
  tpa_floor: %.4f,
  tpa_max: %.4f,
  tpa_k: %.6f,
  max_height: %d,
  species: "Pinus ponderosa"
};

exports.SHRUB = {
  // Chronosequence: age = 2 + 17 * ((h - 0.5) / 9.5)
  //                 biomass_Mg_ha = 13.40 * ln(age) - 7.020
  chrono_a: 13.40,
  chrono_b: -7.020,
  min_height: %.1f,
  max_height: %.1f,
  min_age: %d,
  max_age: %d,
  species: "Salix spp."
};

exports.HERBACEOUS = {
  biomass_Mg_ha: %.1f,
  biomass_Mg_pixel: %.4f,
  species: "Herbaceous vegetation"
};

// Helper: predict deciduous or evergreen biomass (Mg/ha) from height
// Note: Gamma GLM with log-link, so tree_kg = exp(b0 + b1*h + b2*h^2)
exports.predictTreeBiomass = function(height_m, classParams) {
  var tree_kg = Math.exp(classParams.b0 + classParams.b1 * height_m +
                          classParams.b2 * height_m * height_m);
  var tpa = classParams.tpa_floor +
            (classParams.tpa_max - classParams.tpa_floor) *
            Math.exp(-classParams.tpa_k * height_m);
  tpa = Math.max(tpa, 1);
  return (tree_kg * tpa) / 1000 / 0.404686;
};

// Helper: predict willow biomass (Mg/ha) from height
exports.predictWillowBiomass = function(height_m) {
  var p = exports.SHRUB;
  var h = Math.max(p.min_height, Math.min(p.max_height, height_m));
  var age = p.min_age + (p.max_age - p.min_age) *
            ((h - p.min_height) / (p.max_height - p.min_height));
  var biomass = p.chrono_a * Math.log(age) + p.chrono_b;
  return Math.max(0, biomass);
};
',
                  # Deciduous
                  cottonwood_export$b0, cottonwood_export$b1, cottonwood_export$b2,
                  if (!is.null(cottonwood_export$tpa_coefs)) cottonwood_export$tpa_coefs["tpa_floor"] else 60,
                  if (!is.null(cottonwood_export$tpa_coefs)) cottonwood_export$tpa_coefs["tpa_max"] else 160,
                  if (!is.null(cottonwood_export$tpa_coefs)) cottonwood_export$tpa_coefs["k"] else 0.03,
                  cottonwood_export$max_height,
                  # Evergreen
                  ponderosa_export$b0, ponderosa_export$b1, ponderosa_export$b2,
                  if (!is.null(ponderosa_export$tpa_coefs)) ponderosa_export$tpa_coefs["tpa_floor"] else 40,
                  if (!is.null(ponderosa_export$tpa_coefs)) ponderosa_export$tpa_coefs["tpa_max"] else 120,
                  if (!is.null(ponderosa_export$tpa_coefs)) ponderosa_export$tpa_coefs["k"] else 0.03,
                  ponderosa_export$max_height,
                  # Shrub
                  willow_export$min_height, willow_export$max_height,
                  willow_export$min_age, willow_export$max_age,
                  # Herbaceous
                  herbaceous_export$biomass_Mg_ha, herbaceous_export$biomass_Mg_pixel
)

writeLines(gee_js, "gee_biomass_coefficients.js")
cat("GEE JavaScript module saved: gee_biomass_coefficients.js\n\n")

# ==============================================================================
# 4. PUBLICATION FIGURES: ALL-CLASS COMPARISON
# ==============================================================================

# --- 4A: All classes biomass curves on one plot ---
height_all <- seq(0.5, 35, by = 0.5)

# Cottonwood
cw_curve <- cottonwood_export$biomass_curve %>%
  select(height_m, biomass_Mg_ha) %>%
  filter(height_m <= 35) %>%
  mutate(class = "Deciduous (cottonwood)")

# Pine
pine_curve <- ponderosa_export$biomass_curve %>%
  select(height_m, biomass_Mg_ha) %>%
  filter(height_m <= 35) %>%
  mutate(class = "Evergreen (ponderosa)")

# Willow
willow_curve <- tibble(
  height_m = seq(0.5, 10, by = 0.1),
  biomass_Mg_ha = sapply(seq(0.5, 10, by = 0.1), willow_export$biomass_fn),
  class = "Shrub (willow)"
)

# Herbaceous
herb_curve <- tibble(
  height_m = seq(0, 5, by = 0.1),
  biomass_Mg_ha = herbaceous_export$biomass_Mg_ha,
  class = "Herbaceous"
)

all_curves <- bind_rows(cw_curve, pine_curve, willow_curve, herb_curve)

p_all_classes <- ggplot(all_curves, aes(x = height_m, y = biomass_Mg_ha, color = class)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c(
    "Deciduous (cottonwood)" = "#2E7D32",
    "Evergreen (ponderosa)"  = "#1565C0",
    "Shrub (willow)"         = "#7B1FA2",
    "Herbaceous"             = "#F57F17"
  )) +
  labs(title = "Aboveground Biomass by Vegetation Class",
       subtitle = "Height-dependent estimates for Intermountain West riparian zones",
       x = "Canopy Height (m)",
       y = "Aboveground Biomass (Mg/ha)",
       color = "Vegetation Class") +
  theme_minimal(base_size = 13) +
  theme(legend.position = c(0.25, 0.8),
        legend.background = element_rect(fill = "white", color = "gray80"))

print(p_all_classes)
ggsave("fig_all_classes_biomass.png", p_all_classes, width = 10, height = 7, dpi = 300)

# --- 4B: Per-pixel version ---
all_curves_pixel <- all_curves %>%
  mutate(biomass_Mg_pixel = biomass_Mg_ha * 0.01)

p_all_pixel <- ggplot(all_curves_pixel, aes(x = height_m, y = biomass_Mg_pixel, color = class)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c(
    "Deciduous (cottonwood)" = "#2E7D32",
    "Evergreen (ponderosa)"  = "#1565C0",
    "Shrub (willow)"         = "#7B1FA2",
    "Herbaceous"             = "#F57F17"
  )) +
  labs(title = "Per-Pixel Aboveground Biomass by Vegetation Class",
       subtitle = "10m × 10m pixel estimates",
       x = "Canopy Height (m)",
       y = "Biomass (Mg per pixel)",
       color = "Vegetation Class") +
  theme_minimal(base_size = 13) +
  theme(legend.position = c(0.25, 0.8),
        legend.background = element_rect(fill = "white", color = "gray80"))

print(p_all_pixel)
ggsave("fig_all_classes_per_pixel.png", p_all_pixel, width = 10, height = 7, dpi = 300)

# --- 4C: TPA comparison between tree classes ---
cw_tpa_curve <- cottonwood_export$biomass_curve %>%
  select(height_m, TPA) %>%
  mutate(class = "Cottonwood")

pine_tpa_curve <- ponderosa_export$biomass_curve %>%
  select(height_m, TPA) %>%
  mutate(class = "Ponderosa")

tpa_comparison <- bind_rows(cw_tpa_curve, pine_tpa_curve)

p_tpa_compare <- ggplot(tpa_comparison, aes(x = height_m, y = TPA, color = class)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("Cottonwood" = "#2E7D32", "Ponderosa" = "#1565C0")) +
  labs(title = "TPA Scaling Functions: Cottonwood vs Ponderosa",
       x = "Canopy Height (m)", y = "Trees per Acre (TPA)",
       color = "Species") +
  theme_minimal(base_size = 13) +
  theme(legend.position = c(0.8, 0.8))

print(p_tpa_compare)
ggsave("fig_tpa_comparison.png", p_tpa_compare, width = 9, height = 6, dpi = 300)

# ==============================================================================
# 5. COMPREHENSIVE COEFFICIENT TABLE
# ==============================================================================

coef_table <- tibble(
  Class = c("Deciduous", "Evergreen", "Shrub", "Herbaceous"),
  Species = c("Populus trichocarpa", "Pinus ponderosa", "Salix spp.", "Herbaceous"),
  Model = c("Quadratic + TPA", "Quadratic + TPA", "Chronosequence", "Constant"),
  b0 = c(cottonwood_export$b0, ponderosa_export$b0, NA, NA),
  b1 = c(cottonwood_export$b1, ponderosa_export$b1, NA, NA),
  b2 = c(cottonwood_export$b2, ponderosa_export$b2, NA, NA),
  TPA_floor = c(
    if (!is.null(cottonwood_export$tpa_coefs)) cottonwood_export$tpa_coefs["tpa_floor"] else NA,
    if (!is.null(ponderosa_export$tpa_coefs)) ponderosa_export$tpa_coefs["tpa_floor"] else NA,
    NA, NA
  ),
  TPA_max = c(
    if (!is.null(cottonwood_export$tpa_coefs)) cottonwood_export$tpa_coefs["tpa_max"] else NA,
    if (!is.null(ponderosa_export$tpa_coefs)) ponderosa_export$tpa_coefs["tpa_max"] else NA,
    NA, NA
  ),
  TPA_k = c(
    if (!is.null(cottonwood_export$tpa_coefs)) cottonwood_export$tpa_coefs["k"] else NA,
    if (!is.null(ponderosa_export$tpa_coefs)) ponderosa_export$tpa_coefs["k"] else NA,
    NA, NA
  ),
  n_trees = c(cottonwood_export$n_trees, ponderosa_export$n_trees, NA, NA),
  Bayes_R2 = c(cottonwood_export$bayes_r2[1, "Estimate"],
               ponderosa_export$bayes_r2[1, "Estimate"], NA, NA),
  Height_range = c("1-35m", "1-40m", "0.5-10m", "N/A"),
  Fixed_Mg_ha = c(NA, NA, NA, herbaceous_export$biomass_Mg_ha)
)

cat("\n=== Complete Coefficient Table ===\n")
print(coef_table, width = Inf)

write.csv(coef_table, "biomass_model_coefficients.csv", row.names = FALSE)
cat("\nCoefficient table saved: biomass_model_coefficients.csv\n")