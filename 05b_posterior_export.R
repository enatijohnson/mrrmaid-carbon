################################################################################
# 05b_posterior_export.R
# Run BEFORE 05_gee_export
#
# Extracts posterior covariance matrices from the brms allometric models
# for delta-method uncertainty propagation in GEE.
#
# WHY:
#   The current biomass prediction interval only captures height uncertainty
#   (from the conformal CHM prediction interval). It treats the allometric
#   coefficients (b0, b1, b2) as fixed point estimates. In reality, these
#   coefficients have their own uncertainty from the Bayesian model fit —
#   the posterior distribution over (b0, b1, b2) captures how confident
#   we are in the height-to-biomass relationship itself.
#
#   By exporting the posterior covariance matrix, GEE can compute the
#   allometric uncertainty per pixel using the delta method, without
#   needing to run brms or sample from the posterior in GEE.
#
# DELTA METHOD (for Gamma GLM with log link):
#   The model is: log(mu) = b0 + b1*h + b2*h^2
#   For a given height h, the design vector is: x = [1, h, h^2]
#   The variance of the log-scale prediction is: Var[log(mu)] = x' Sigma x
#   where Sigma is the 3x3 posterior covariance matrix of (b0, b1, b2).
#
#   On the response scale (biomass in kg):
#     CV_allometric = sqrt(exp(x' Sigma x) - 1)
#     SD_allometric = mu * CV_allometric
#   where mu = exp(b0 + b1*h + b2*h^2) is the point prediction.
#
#   This is an approximation that works well when the posterior is
#   approximately normal. We validate by comparing to Monte Carlo draws.
#
# OUTPUTS:
#   posterior_covariance_deciduous.csv   — 3×3 covariance matrix
#   posterior_covariance_evergreen.csv   — 3×3 covariance matrix
#   posterior_summary.csv                — mean, sd, quantiles per coefficient
#   gee_posterior_module.js              — GEE JavaScript module with matrices
#   fig_posterior_correlation.png        — correlation heatmap
#   fig_delta_vs_monte_carlo.png        — validation comparison
#
# REQUIRES:
#   cw_model_quad and pine_model_quad brms objects in the environment
#   (from scripts 02_cottonwood_deciduous.R and 03_ponderosa_evergreen.R)
################################################################################

library(tidyverse)
library(brms)
library(posterior)

cat("\n########## POSTERIOR COVARIANCE EXTRACTION ##########\n\n")

# ==============================================================================
# 1. EXTRACT POSTERIOR DRAWS AND COVARIANCE MATRICES
# ==============================================================================

extract_posterior_info <- function(model, label) {
  cat("--- ", label, " ---\n")
  
  # Posterior summary
  coefs <- fixef(model)
  cat("Coefficient summary:\n")
  print(coefs)
  
  # Full posterior draws (fixed effects only)
  draws <- as_draws_matrix(model, variable = "^b_", regex = TRUE)
  
  # The column names from brms are b_Intercept, b_HT_m, b_IHT_mE2
  # Rename for clarity.
  draws_df <- as.data.frame(draws) %>%
    select(starts_with("b_"))
  
  # Identify the three coefficient columns
  coef_cols <- colnames(draws_df)
  cat("Posterior draw columns:", paste(coef_cols, collapse = ", "), "\n")
  
  # Covariance matrix of the posterior draws
  cov_matrix <- cov(draws_df)
  cat("Posterior covariance matrix:\n")
  print(round(cov_matrix, 8))
  
  # Correlation matrix
  cor_matrix <- cor(draws_df)
  cat("\nPosterior correlation matrix:\n")
  print(round(cor_matrix, 4))
  cat("\n")
  
  return(list(
    label = label,
    coefs = coefs,
    draws = draws_df,
    cov_matrix = cov_matrix,
    cor_matrix = cor_matrix,
    coef_names = coef_cols
  ))
}

cw_post  <- extract_posterior_info(cw_model_quad, "Deciduous (Cottonwood)")
pine_post <- extract_posterior_info(pine_model_quad, "Evergreen (Ponderosa)")


# ==============================================================================
# 2. EXPORT COVARIANCE MATRICES AS CSV
# ==============================================================================
# Each CSV is a 3×3 matrix with row/column names matching the coefficient order:
# Intercept, HT_m, HT_m^2. GEE can read these as ee.FeatureCollection or
# parse them client-side in the JavaScript module.

export_cov_csv <- function(post_info, filename) {
  cov_df <- as.data.frame(post_info$cov_matrix)
  
  # Standardize column names for GEE consumption
  colnames(cov_df) <- c("Intercept", "HT_m", "HT_m2")
  cov_df$parameter <- c("Intercept", "HT_m", "HT_m2")
  cov_df <- cov_df %>% select(parameter, everything())
  
  write_csv(cov_df, filename)
  cat("Saved:", filename, "\n")
}

export_cov_csv(cw_post, "posterior_covariance_deciduous.csv")
export_cov_csv(pine_post, "posterior_covariance_evergreen.csv")

# Combined summary table
summary_table <- bind_rows(
  as.data.frame(cw_post$coefs) %>%
    rownames_to_column("parameter") %>%
    mutate(class = "deciduous"),
  as.data.frame(pine_post$coefs) %>%
    rownames_to_column("parameter") %>%
    mutate(class = "evergreen")
)
write_csv(summary_table, "posterior_summary.csv")
cat("Saved: posterior_summary.csv\n")


# ==============================================================================
# 3. GENERATE GEE JAVASCRIPT MODULE
# ==============================================================================
# Embeds the covariance matrices directly in JavaScript so GEE can compute
# delta-method uncertainty without loading external assets.

format_matrix_js <- function(mat, indent) {
  rows <- apply(mat, 1, function(row) {
    paste0(indent, "  [", paste(formatC(row, format = "e", digits = 8), collapse = ", "), "]")
  })
  paste0(indent, "[\n", paste(rows, collapse = ",\n"), "\n", indent, "]")
}

gee_js <- sprintf('/**
 * gee_posterior_uncertainty.js
 * ----------------------------
 * Posterior covariance matrices from brms allometric models for
 * delta-method uncertainty propagation.
 *
 * For a pixel with height h:
 *   x = [1, h, h^2]
 *   var_log_mu = x[0]*x[0]*S[0][0] + x[0]*x[1]*S[0][1]*2 + x[0]*x[2]*S[0][2]*2
 *              + x[1]*x[1]*S[1][1] + x[1]*x[2]*S[1][2]*2
 *              + x[2]*x[2]*S[2][2]
 *   cv_allometric = sqrt(exp(var_log_mu) - 1)
 *   sd_allometric_Mg_ha = biomass_Mg_ha * cv_allometric
 *
 * Generated by 05b_posterior_export.R.
 */
 
// Posterior covariance matrix: Cov(Intercept, HT_m, HT_m^2)
// Deciduous (cottonwood) — Bayesian Gamma GLM (log link, quadratic)
exports.COV_DECIDUOUS = %s;
 
// Evergreen (ponderosa) — Bayesian Gamma GLM (log link, quadratic)
exports.COV_EVERGREEN = %s;
 
// Posterior means (same as point estimates in the allometric models)
exports.MEANS_DECIDUOUS = [%.8f, %.8f, %.8f];
exports.MEANS_EVERGREEN = [%.8f, %.8f, %.8f];
 
/**
 * Compute allometric CV for a height image using the delta method.
 * Returns a single-band image of coefficient of variation (dimensionless).
 *
 * @param {ee.Image} heightImage - Height in meters.
 * @param {Array} covMatrix - 3x3 covariance matrix [[S00,S01,S02],[S10,S11,S12],[S20,S21,S22]].
 * @returns {ee.Image} CV image (multiply by biomass_Mg_ha to get SD in Mg/ha).
 */
exports.computeAllometricCV = function(heightImage, covMatrix) {
  var h = heightImage;
  var h2 = h.multiply(h);
 
  // Compute x^T Sigma x using the symmetric matrix expansion.
  // For x = [1, h, h^2] and symmetric S:
  //   x^T S x = S00 + 2*S01*h + 2*S02*h^2 + S11*h^2 + 2*S12*h^3 + S22*h^4
  var S = covMatrix;
  var varLogMu = ee.Image.constant(S[0][0])
    .add(h.multiply(2 * S[0][1]))
    .add(h2.multiply(2 * S[0][2]))
    .add(h2.multiply(S[1][1]))
    .add(h2.multiply(h).multiply(2 * S[1][2]))
    .add(h2.multiply(h2).multiply(S[2][2]));
 
  // CV = sqrt(exp(var_log_mu) - 1)
  // For small var_log_mu this is approximately sqrt(var_log_mu)
  var cv = varLogMu.exp().subtract(1).max(0).sqrt();
 
  return cv.rename("allometric_cv");
};
',
                  format_matrix_js(cw_post$cov_matrix, ""),
                  format_matrix_js(pine_post$cov_matrix, ""),
                  cw_post$coefs[1, "Estimate"], cw_post$coefs[2, "Estimate"], cw_post$coefs[3, "Estimate"],
                  pine_post$coefs[1, "Estimate"], pine_post$coefs[2, "Estimate"], pine_post$coefs[3, "Estimate"]
)

writeLines(gee_js, "gee_posterior_uncertainty.js")
cat("Saved: gee_posterior_uncertainty.js\n\n")

# ==============================================================================
# 4. VALIDATE DELTA METHOD VS MONTE CARLO
# ==============================================================================
# Compare the delta-method approximation to full Monte Carlo sampling
# from the posterior at a range of heights. If the delta method is accurate,
# the two CV curves should closely overlap.

validate_delta_vs_mc <- function(post_info, label, n_mc = 4000) {
  heights <- seq(1, 35, by = 0.5)
  draws <- as.matrix(post_info$draws)
  cov_mat <- post_info$cov_matrix
  
  results <- map_dfr(heights, function(h) {
    x <- c(1, h, h^2)
    
    # Delta method CV
    var_log_mu <- as.numeric(t(x) %*% cov_mat %*% x)
    cv_delta <- sqrt(exp(var_log_mu) - 1)
    
    # Monte Carlo: sample biomass_kg from each posterior draw
    log_mu_draws <- draws %*% x   # n_draws × 1
    mu_draws <- exp(log_mu_draws)
    cv_mc <- sd(mu_draws) / mean(mu_draws)
    
    # Also get the mean prediction for context
    mean_coefs <- post_info$coefs[, "Estimate"]
    mu_point <- exp(sum(mean_coefs * x))
    
    tibble(
      height_m = h,
      cv_delta = cv_delta,
      cv_mc = cv_mc,
      sd_delta_kg = mu_point * cv_delta,
      sd_mc_kg = sd(mu_draws),
      mean_kg = mu_point
    )
  })
  
  results$class <- label
  return(results)
}

cw_validation  <- validate_delta_vs_mc(cw_post, "Deciduous")
pine_validation <- validate_delta_vs_mc(pine_post, "Evergreen")
all_validation <- bind_rows(cw_validation, pine_validation)

# Print summary
cat("=== Delta Method vs Monte Carlo Validation ===\n")
cat("Max absolute CV difference (Deciduous):",
    max(abs(cw_validation$cv_delta - cw_validation$cv_mc)), "\n")
cat("Max absolute CV difference (Evergreen):",
    max(abs(pine_validation$cv_delta - pine_validation$cv_mc)), "\n\n")

# --- Plot: CV comparison ---
p_cv <- ggplot(all_validation, aes(x = height_m)) +
  geom_line(aes(y = cv_delta, color = "Delta method"), linewidth = 1.2) +
  geom_line(aes(y = cv_mc, color = "Monte Carlo"), linewidth = 1.2, linetype = "dashed") +
  facet_wrap(~class, scales = "free_y") +
  scale_color_manual(values = c("Delta method" = "#1a73e8", "Monte Carlo" = "#e8711a")) +
  labs(
    title = "Allometric Uncertainty: Delta Method vs Monte Carlo",
    subtitle = "CV of biomass prediction from posterior parameter uncertainty alone",
    x = "Canopy Height (m)",
    y = "Coefficient of Variation",
    color = "Method"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

print(p_cv)
ggsave("fig_delta_vs_monte_carlo.png", p_cv, width = 10, height = 5, dpi = 300)


# --- Plot: SD in kg at each height ---
p_sd <- ggplot(all_validation, aes(x = height_m)) +
  geom_ribbon(aes(ymin = mean_kg - sd_mc_kg, ymax = mean_kg + sd_mc_kg),
              fill = "#e8711a", alpha = 0.2) +
  geom_ribbon(aes(ymin = mean_kg - sd_delta_kg, ymax = mean_kg + sd_delta_kg),
              fill = "#1a73e8", alpha = 0.2) +
  geom_line(aes(y = mean_kg), linewidth = 1) +
  facet_wrap(~class, scales = "free_y") +
  labs(
    title = "Per-Tree Biomass with Allometric Uncertainty",
    subtitle = "Blue = delta method ±1 SD | Orange = Monte Carlo ±1 SD",
    x = "Canopy Height (m)",
    y = "Tree Biomass (kg)"
  ) +
  theme_minimal(base_size = 13)

print(p_sd)
ggsave("fig_allometric_uncertainty_bands.png", p_sd, width = 10, height = 5, dpi = 300)


# --- Plot: Posterior correlation heatmap ---
plot_cor_heatmap <- function(cor_mat, label) {
  cor_df <- as.data.frame(cor_mat) %>%
    mutate(param1 = rownames(cor_mat)) %>%
    pivot_longer(-param1, names_to = "param2", values_to = "correlation")
  
  # Clean parameter names
  cor_df <- cor_df %>%
    mutate(
      param1 = str_replace_all(param1, c("b_Intercept" = "b0", "b_HT_m" = "b1", ".*HT_mE2.*" = "b2")),
      param2 = str_replace_all(param2, c("b_Intercept" = "b0", "b_HT_m" = "b1", ".*HT_mE2.*" = "b2"))
    )
  
  ggplot(cor_df, aes(x = param1, y = param2, fill = correlation)) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = round(correlation, 3)), size = 5, fontface = "bold") +
    scale_fill_gradient2(low = "#d73027", mid = "white", high = "#4575b4",
                         midpoint = 0, limits = c(-1, 1)) +
    labs(title = paste("Posterior Correlation —", label),
         x = "", y = "") +
    theme_minimal(base_size = 13) +
    theme(panel.grid = element_blank())
}

p_cor_cw <- plot_cor_heatmap(cw_post$cor_matrix, "Deciduous")
p_cor_pine <- plot_cor_heatmap(pine_post$cor_matrix, "Evergreen")

library(patchwork)
p_cor_combined <- p_cor_cw + p_cor_pine
print(p_cor_combined)
ggsave("fig_posterior_correlation.png", p_cor_combined, width = 12, height = 5, dpi = 300)