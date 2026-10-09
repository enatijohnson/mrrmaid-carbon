################################################################################
# 04_willow_shrub_herbaceous.R
#
# Willow (shrub class): height -> age -> biomass chronosequence
# Herbaceous: constant 1.5 Mg/ha
# Requires: 01_data_preparation.R objects in environment.
################################################################################

# ==============================================================================
# 1. WILLOW HEIGHT-AGE-BIOMASS CHRONOSEQUENCE
# ==============================================================================
#
# Age-biomass equation: Biomass(Mg/ha) = 13.40 × ln(age) - 7.020
# Height-to-age: linear mapping from [0.5m, 10m] -> [2yr, 19yr]
# No artificial biomass floor; equation gives ~2.3 Mg/ha at age 2 (min).

# --- Parameters ---
MIN_WILLOW_HEIGHT <- 0.5
MAX_WILLOW_HEIGHT <- 10.0
MIN_AGE <- 2
MAX_AGE <- 19

# --- Functions ---
height_to_age <- function(height_m) {
  height_m <- pmax(MIN_WILLOW_HEIGHT, pmin(MAX_WILLOW_HEIGHT, height_m))
  MIN_AGE + (MAX_AGE - MIN_AGE) *
    ((height_m - MIN_WILLOW_HEIGHT) / (MAX_WILLOW_HEIGHT - MIN_WILLOW_HEIGHT))
}

age_to_biomass_Mg_ha <- function(age) {
  pmax(0, 13.40 * log(age) - 7.020)
}

willow_height_to_biomass <- function(height_m) {
  age_to_biomass_Mg_ha(height_to_age(height_m))
}

# --- Validation against empirical data ---
cat("Validation against published empirical data:\n")
cat("  Age  5yr: equation =", round(13.40 * log(5) - 7.020, 1),
    "Mg/ha | observed = 14.1 Mg/ha\n")
cat("  Age 12yr: equation =", round(13.40 * log(12) - 7.020, 1),
    "Mg/ha | observed = 28.8 Mg/ha\n")
cat("  Age 19yr: equation =", round(13.40 * log(19) - 7.020, 1),
    "Mg/ha | observed = 31.6 Mg/ha\n\n")

# --- Diagnostic table ---
willow_diag <- tibble(
  height_m = seq(0.5, 10, by = 0.5)
) %>%
  mutate(
    age_yr = round(height_to_age(height_m), 1),
    biomass_Mg_ha = round(willow_height_to_biomass(height_m), 2),
    biomass_Mg_pixel = round(biomass_Mg_ha * 0.01, 4)
  )

cat("Willow height-to-biomass lookup:\n")
print(willow_diag, n = Inf)

# --- Visualization ---
willow_plot_data <- tibble(
  height_m = seq(0, 12, by = 0.1),
  biomass_Mg_ha = sapply(seq(0, 12, by = 0.1), willow_height_to_biomass)
)

p_willow <- ggplot(willow_plot_data, aes(x = height_m, y = biomass_Mg_ha)) +
  geom_line(color = "purple", linewidth = 1.3) +
  geom_vline(xintercept = 4.12, linetype = "dashed", color = "blue", alpha = 0.6) +
  geom_vline(xintercept = c(0.5, 10.0), linetype = "dotted", color = "red", alpha = 0.5) +
  annotate("text", x = 4.3, y = 5, label = "Mean height\n(4.12m)",
           color = "blue", size = 3, hjust = 0) +
  labs(title = "Willow: Height-Biomass Relationship",
       subtitle = "Chronosequence: height → age → Biomass(Mg/ha) = 13.40 × ln(age) − 7.020",
       x = "Canopy Height (m)", y = "Aboveground Biomass (Mg/ha)") +
  theme_minimal()

print(p_willow)
ggsave("fig_willow_biomass.png", p_willow, width = 9, height = 6, dpi = 200)

# --- Store for export ---
willow_export <- list(
  class = "shrub",
  species = "Salix spp. (willow)",
  method = "height_age_biomass_chronosequence",
  equation_age_biomass = "Biomass_Mg_ha = 13.40 * ln(age) - 7.020",
  equation_height_age = sprintf("age = %d + (%d - %d) * ((h - %.1f) / (%.1f - %.1f))",
                                 MIN_AGE, MAX_AGE, MIN_AGE,
                                 MIN_WILLOW_HEIGHT, MAX_WILLOW_HEIGHT, MIN_WILLOW_HEIGHT),
  min_height = MIN_WILLOW_HEIGHT,
  max_height = MAX_WILLOW_HEIGHT,
  min_age = MIN_AGE,
  max_age = MAX_AGE,
  chronosequence_coefs = list(a = 13.40, b = -7.020),
  biomass_fn = willow_height_to_biomass
)


# ==============================================================================
# 2. HERBACEOUS CLASS
# ==============================================================================


# Constant value from Dwire et al. (2004) — riparian meadows, NE Oregon
HERBACEOUS_BIOMASS_MG_HA    <- 1.5
HERBACEOUS_BIOMASS_MG_PIXEL <- HERBACEOUS_BIOMASS_MG_HA * 0.01

cat("Herbaceous biomass: constant 1.5 Mg/ha (0.015 Mg per 10m pixel)\n")
cat("Source: Dwire et al. (2004) — montane riparian meadows\n\n")

herbaceous_export <- list(
  class = "herbaceous",
  species = "Herbaceous vegetation",
  method = "constant",
  biomass_Mg_ha = HERBACEOUS_BIOMASS_MG_HA,
  biomass_Mg_pixel = HERBACEOUS_BIOMASS_MG_PIXEL,
  citation = "Dwire et al. (2004)"
)

