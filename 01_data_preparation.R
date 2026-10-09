################################################################################
# 01_data_preparation.R
#
# FIA data download, clipping, table extraction, and shared helper functions.
################################################################################

library(rFIA)
library(dplyr)
library(tidyr)
library(sf)
library(ggplot2)
library(mgcv)
library(tidyverse)

# ==============================================================================
# 1. LOAD FIA DATA
# ==============================================================================

setwd("~/R/FIA_analysis")

options(timeout = 3600)

fia_data <- getFIA(
  states = c("ID", "OR", "MT", "WA", "WY", "CO", "UT", "NV", "CA"),
  dir = "fia_data/",
  load = TRUE
)

study_area <- st_read("US_Sagebrush_Biome_2019.shp")
fia_clipped <- clipFIA(fia_data, mask = study_area)

# ==============================================================================
# 2. EXTRACT TABLES AND PREPARE TREE DATA
# ==============================================================================

tree_df <- fia_clipped$TREE
cond_df <- fia_clipped$COND
plot_df <- fia_clipped$PLOT

# Aboveground biomass only (DRYBIO_AG), convert units
tree_df <- tree_df %>%
  mutate(
    biomass_kg = DRYBIO_AG * 0.453592,   # lb -> kg
    HT_m       = HT * 0.3048            # ft -> m
  )

# ==============================================================================
# 3. SHARED HELPER FUNCTIONS
# ==============================================================================

# Biomass conversion: (biomass_kg * TPA) / 1000 / 0.404686 -> Mg/ha
calc_biomass_Mg_ha <- function(tree_biomass_kg, tpa) {
  (tree_biomass_kg * tpa) / 1000 / 0.404686
}

# Monotonicity constraint: walk up height axis, adjust TPA upward where
# the biomass product would otherwise decrease.
constrain_monotonic_tpa <- function(bin_data, tpa_column, label,
                                    predict_fn, biomass_fn) {
  result <- bin_data %>%
    arrange(mean_HT) %>%
    mutate(
      TPA_raw = .data[[tpa_column]],
      TPA_constrained = TPA_raw,
      tree_biomass_kg = predict_fn(mean_HT),
      was_adjusted = FALSE
    )
  
  result$biomass_Mg_ha <- biomass_fn(result$tree_biomass_kg, result$TPA_constrained)
  
  for (i in 2:nrow(result)) {
    if (result$biomass_Mg_ha[i] < result$biomass_Mg_ha[i - 1]) {
      required_biomass <- result$biomass_Mg_ha[i - 1]
      required_TPA <- (required_biomass * 1000 * 0.404686) / result$tree_biomass_kg[i]
      result$TPA_constrained[i] <- required_TPA * 1.01
      result$biomass_Mg_ha[i] <- biomass_fn(result$tree_biomass_kg[i],
                                             result$TPA_constrained[i])
      result$was_adjusted[i] <- TRUE
    }
  }
  
  result$approach <- label
  return(result)
}

# Full TPA scaling pipeline for a tree class:
#   1. Build per-tree data (tree heights + parent plot TPA)
#   2. Bin by height, compute local medians
#   3. Blend with global anchor (split-weight)
#   4. Apply monotonicity constraint
#   5. Smooth with parametric fit
#   6. Return final TPA function + coefficients
run_tpa_scaling <- function(tree_data, tpa_data, predict_fn,
                            bin_edges, w_low = 0.5, w_high = 0.7,
                            ht_threshold = 20, species_label = "Species") {
  
  # Per-tree merge: each tree gets its parent plot's TPA
  merged <- tree_data %>%
    left_join(tpa_data %>% select(PLT_CN, TPA), by = "PLT_CN")
  
  n_before <- nrow(merged)
  n_na <- sum(is.na(merged$TPA))
  
  if (n_na > 0) {
    cat(sprintf("\nNote: %d of %d trees (%.1f%%) had no matching TPA and were dropped.\n",
                n_na, n_before, 100 * n_na / n_before))
    cat("  (Plots present in tree table but absent from tpa() output.)\n")
    merged <- merged %>% filter(!is.na(TPA))
  }
  
  global_median_TPA <- median(merged$TPA, na.rm = TRUE)
  global_mean_TPA   <- mean(merged$TPA, na.rm = TRUE)
  
  cat(sprintf("\n=== TPA Scaling: %s ===\n", species_label))
  cat("Trees:", nrow(merged), " | Plots:", n_distinct(merged$PLT_CN), "\n")
  cat("Global median TPA (per-tree):", round(global_median_TPA, 1), "\n")
  cat("Global mean TPA (per-tree):", round(global_mean_TPA, 1), "\n")
  
  # Bin by height
  binned <- merged %>%
    mutate(height_bin = cut(HT_m, breaks = bin_edges, include.lowest = TRUE))
  
  bin_stats <- binned %>%
    group_by(height_bin) %>%
    summarise(
      n_trees = n(),
      n_plots = n_distinct(PLT_CN),
      mean_HT = mean(HT_m),
      local_median_TPA = median(TPA),
      .groups = "drop"
    ) %>%
    arrange(mean_HT)
  
  cat("\nBin statistics:\n")
  print(bin_stats, n = Inf, width = Inf)
  
  # Split-weight blending
  bin_stats <- bin_stats %>%
    mutate(
      blended_TPA = ifelse(
        mean_HT < ht_threshold,
        w_low * global_median_TPA + (1 - w_low) * local_median_TPA,
        w_high * global_median_TPA + (1 - w_high) * local_median_TPA
      )
    )
  
  # Monotonicity constraint
  constrained <- constrain_monotonic_tpa(
    bin_stats, "blended_TPA",
    sprintf("Constrained split (w=%.1f/<%dm, w=%.1f/>=%dm)",
            w_low, ht_threshold, w_high, ht_threshold),
    predict_fn, calc_biomass_Mg_ha
  )
  
  cat("\nConstrained TPA values:\n")
  constrained %>%
    select(height_bin, mean_HT, n_trees, TPA_raw, TPA_constrained,
           was_adjusted, biomass_Mg_ha) %>%
    mutate(across(where(is.numeric), ~round(., 2))) %>%
    print(n = Inf, width = Inf)
  
  # Parametric smooth: TPA(h) = floor + (max - floor) * exp(-k * h)
  smooth_d <- constrained %>%
    select(mean_HT, TPA_constrained) %>%
    filter(!is.na(TPA_constrained))
  
  param_fit <- tryCatch({
    nls(
      TPA_constrained ~ tpa_floor + (tpa_max - tpa_floor) * exp(-k * mean_HT),
      data = smooth_d,
      start = list(tpa_floor = min(smooth_d$TPA_constrained),
                   tpa_max = max(smooth_d$TPA_constrained),
                   k = 0.05),
      lower = c(1, 10, 0.001),
      algorithm = "port"
    )
  }, error = function(e) {
    cat("Parametric fit failed:", e$message, "\n")
    NULL
  })
  
  if (!is.null(param_fit)) {
    cat("\nSmoothed TPA equation:\n")
    cat(sprintf("  TPA(h) = %.2f + (%.2f - %.2f) * exp(-%.4f * h)\n",
                coef(param_fit)["tpa_floor"], coef(param_fit)["tpa_max"],
                coef(param_fit)["tpa_floor"], coef(param_fit)["k"]))
  }
  
  return(list(
    merged = merged,
    bin_stats = bin_stats,
    constrained = constrained,
    param_fit = param_fit,
    global_median_TPA = global_median_TPA,
    w_low = w_low,
    w_high = w_high,
    ht_threshold = ht_threshold
  ))
}
