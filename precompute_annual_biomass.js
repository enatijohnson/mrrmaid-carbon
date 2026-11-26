/**
 * precompute_annual_biomass.js
 * ----------------------------
 * Pre-computes annual biomass rasters (2017-2025) by combining:
 *   - Predicted CHM ensemble (v2, Int16 × 100, from mrrmaid-carbon)
 *   - Vegetation classification with soft probabilities (v3, from mrrmaid-carbon)
 *   - FIA-derived allometric models with height-dependent TPA scaling
 *
 * For each year, exports a 9-band Int16 × 10 image (divide by 10 for Mg/ha):
 *   biomass_mean               — hard class × mean CHM (primary estimate)
 *   biomass_lower              — hard class × PI_Lower CHM
 *   biomass_upper              — hard class × PI_Upper CHM
 *   biomass_weighted           — probability-weighted mixture, mean CHM
 *   biomass_weighted_lower     — probability-weighted mixture, PI_Lower CHM
 *   biomass_weighted_upper     — probability-weighted mixture, PI_Upper CHM
 *   prediction_interval_width  — band 6 minus band 5
 *   allometric_sd              - Bayesian posterior SD on biomass (delta method)
 *   allometric_cv              - coefficient of variation x 1000 (dimensionless)
 * 
 * UNCERTAINTY SOURCES:
 * The biomass output captures three uncertainty sources:
 *    1. CHM prediction interval: height uncertainty from conformal calibration.
 *        Captured in biomass_lower/upper (hard) and weighted_lower/upper (soft).
 *    2. Classification uncertainty: captured by the probability-weighted bands.
 *        Pixels with ambiguous class probabilities get blended estimates.
 *    3. Allometric parameter uncertainty: captured in allometric_sd/cv via the
 *        delta method on the brms posterior covariance matrix. This is new
 *        in v3 and quantifies how uncertain we are in the height-to-biomass
 *        relationship itself, independent of how uncertain we are about height.
 * 
 *   The three sources are conceptually independent and can be combined in
 *   downstream analysis. The script does not currently combine them into a
 *   single total uncertainty band. That requires deciding how to propagate
 *   them through the nonlinear allometric function. They are reported
 *   separately so the magnitude of each can be evaluated.
 *
 * Storage: Int16 × 10 (0.1 Mg/ha precision). Covers 0–3276.7 Mg/ha, well
 * above the maximum expected riparian biomass (~300 Mg/ha). The asset_loader
 * module handles rescaling on load — consumers get Float32 in Mg/ha.
 *
 * The hard-class estimates use only the argmax vegetation class for each
 * pixel. The probability-weighted estimates account for classification
 * uncertainty by computing biomass under all four classes and weighting
 * by the ensemble class probabilities.
 *
 * The prediction interval (PI) on biomass comes from propagating the CHM
 * prediction interval through the allometric equations. This captures
 * height uncertainty but not allometric parameter uncertainty (which would
 * require the Bayesian posterior covariance from brms — planned for a
 * future update).
 *
 * Yearly masks are not applied here. The biomass rasters cover the full
 * CHM extent. The app and temporal change scripts apply year-specific
 * masks at load time via asset_loader.applyYearlyMask().
 *
 * CRS: EPSG:4326 at 10 m. Exports to mrrmaid-carbon project.
 */
 
// =============================================================================
// MODULE IMPORTS
// =============================================================================
 
var Loader = require('users/eriknatijohnson/MRRMaidCarbon:finalstack/asset_loader');
var Posterior = require('users/eriknatijohnson/MRRMaidCarbon:finalstack/gee_posterior_uncertainty');
 
// =============================================================================
// CONFIGURATION
// =============================================================================
 
var STUDY_AREA = ee.FeatureCollection(Loader.PATHS.studyArea);
var GEOMETRY = STUDY_AREA.geometry();
 
var SCALE = 10;
var CRS = 'EPSG:4326';
 
var YEARS = Loader.YEARS;
var B = Loader.BAND_NAMES;
 
// =============================================================================
// ALLOMETRIC MODEL COEFFICIENTS (from R brms analysis)
// =============================================================================
// Bayesian Gamma GLM with log-link and quadratic height term:
//   tree_kg = exp(b0 + b1 * h + b2 * h²)
//
// Height-dependent TPA scaling (exponential decay with monotonicity constraint):
//   TPA(h) = floor + (max - floor) * exp(-k * h)
//
// Biomass conversion:
//   Mg/ha = (tree_kg × TPA) / 1000 / 0.404686
//   where 0.404686 converts acres⁻¹ to ha⁻¹
 
var MODELS = {
  evergreen: {   // Class 0 — Pinus ponderosa
    b0: 1.378819,
    b1: 0.347377,
    b2: -0.004520,
    tpa_floor: 103.8374,
    tpa_max: 219.4595,
    tpa_k: 0.165506
  },
  deciduous: {   // Class 1 — Populus trichocarpa
    b0: 1.327521,
    b1: 0.399839,
    b2: -0.006079,
    tpa_floor: 44.9257,
    tpa_max: 228.6310,
    tpa_k: 0.083137
  },
  herbaceous: {  // Class 2 — constant (Dwire et al. 2004)
    biomass_Mg_ha: 1.5
  },
  shrub: {       // Class 3 — Salix spp. (chronosequence)
    chrono_a: 13.40,
    chrono_b: -7.020,
    min_height: 0.5,
    max_height: 10.0,
    min_age: 2,
    max_age: 19
  }
};
 
 
// =============================================================================
// BIOMASS PREDICTION FUNCTIONS
// =============================================================================
// All functions operate on ee.Image objects (server-side pixel algebra).
// Heights are in meters (already rescaled by asset_loader.loadCHM).
 
/**
 * Predict tree biomass (Mg/ha) for evergreen or deciduous class.
 * Gamma GLM log-link: tree_kg = exp(b0 + b1*h + b2*h²)
 * Then scaled by height-dependent TPA and unit conversion.
 */
function predictTreeBiomass(heightImage, params) {
  var tree_kg = heightImage.multiply(params.b1)
      .add(heightImage.pow(2).multiply(params.b2))
      .add(params.b0)
      .exp();
 
  var tpa = heightImage.multiply(-params.tpa_k)
      .exp()
      .multiply(params.tpa_max - params.tpa_floor)
      .add(params.tpa_floor)
      .max(1);
 
  return tree_kg.multiply(tpa).divide(1000).divide(0.404686);
}
 
/**
 * Predict willow biomass (Mg/ha) via chronosequence.
 * height → age (linear interpolation) → biomass (log-linear).
 * Height clamped to [0.5, 10] m; biomass floored at 0.
 */
function predictWillowBiomass(heightImage) {
  var p = MODELS.shrub;
  var h = heightImage.max(p.min_height).min(p.max_height);
 
  var age = h.subtract(p.min_height)
      .divide(p.max_height - p.min_height)
      .multiply(p.max_age - p.min_age)
      .add(p.min_age);
 
  return age.log().multiply(p.chrono_a).add(p.chrono_b).max(0);
}
 
/**
 * Hard-class biomass: assigns each pixel the biomass for its argmax class.
 */
function predictBiomassHardClass(heightImage, vegClass) {
  var evergreen  = predictTreeBiomass(heightImage, MODELS.evergreen);
  var deciduous  = predictTreeBiomass(heightImage, MODELS.deciduous);
  var shrub      = predictWillowBiomass(heightImage);
  var herbaceous = ee.Image.constant(MODELS.herbaceous.biomass_Mg_ha);
 
  return ee.Image(0)
      .where(vegClass.eq(0), evergreen)
      .where(vegClass.eq(1), deciduous)
      .where(vegClass.eq(2), herbaceous)
      .where(vegClass.eq(3), shrub);
}
 
/**
 * Probability-weighted biomass: computes biomass under all four classes
 * and takes the probability-weighted sum. This accounts for classification
 * uncertainty. Pixels where the classifier is unsure between deciduous
 * and shrub (for example) get a blended estimate rather than committing
 * to one class.
 */
function predictBiomassWeighted(heightImage, vegClassImage) {
  var evergreen  = predictTreeBiomass(heightImage, MODELS.evergreen);
  var deciduous  = predictTreeBiomass(heightImage, MODELS.deciduous);
  var shrub      = predictWillowBiomass(heightImage);
  var herbaceous = ee.Image.constant(MODELS.herbaceous.biomass_Mg_ha);
 
  var pEv = vegClassImage.select(B.VEG.probEvergreen);
  var pDe = vegClassImage.select(B.VEG.probDeciduous);
  var pHe = vegClassImage.select(B.VEG.probHerbaceous);
  var pSh = vegClassImage.select(B.VEG.probShrub);
 
  return evergreen.multiply(pEv)
      .add(deciduous.multiply(pDe))
      .add(herbaceous.multiply(pHe))
      .add(shrub.multiply(pSh));
}
 
 
// =============================================================================
// LOAD VEGETATION CLASSIFICATION
// =============================================================================
 
var vegClassImage = Loader.loadVegClassification();
var hardClass = vegClassImage.select(B.VEG.hardClass);
 
 
// =============================================================================
// PROCESS AND EXPORT EACH YEAR
// =============================================================================
 
YEARS.forEach(function(year) {
 
  // Load CHM via asset_loader (handles Int16 ÷ 100 rescaling)
  var chm = Loader.loadCHM(year);
 
  var heightMean  = chm.select(B.CHM.mean);
  var heightLower = chm.select(B.CHM.piLower);
  var heightUpper = chm.select(B.CHM.piUpper);
 
  // --- Hard-class biomass ---
  var biomassMean  = predictBiomassHardClass(heightMean, hardClass)
      .rename(B.BIOMASS.mean);
  var biomassLower = predictBiomassHardClass(heightLower, hardClass)
      .rename(B.BIOMASS.lower);
  var biomassUpper = predictBiomassHardClass(heightUpper, hardClass)
      .rename(B.BIOMASS.upper);
 
  // --- Probability-weighted biomass ---
  var biomassWeighted = predictBiomassWeighted(heightMean, vegClassImage)
      .rename(B.BIOMASS.weighted);
  var biomassWeightedLower = predictBiomassWeighted(heightLower, vegClassImage)
      .rename(B.BIOMASS.weightedLower);
  var biomassWeightedUpper = predictBiomassWeighted(heightUpper, vegClassImage)
      .rename(B.BIOMASS.weightedUpper);
 
  // --- Prediction interval width ---
  var piWidth = biomassWeightedUpper.subtract(biomassWeightedLower)
      .rename(B.BIOMASS.piWidth);
 
   // --- Allometric parameter uncertainty (delta method on brms posterior) ---
  // Compute the height-dependent CV from each tree class's posterior
  // covariance matrix. The CV represents the uncertainty in the biomass
  // prediction arising from uncertainty in the allometric coefficients
  // themselves (b0, b1, b2), independent of height uncertainty.
  //
  // For shrub and herbaceous classes, allometric CV is set to 0 because:
  //   - Shrub uses a deterministic chronosequence (no posterior available)
  //   - Herbaceous is a literature constant
  // These could be added later from published parameter uncertainties.
  //
  // The allometric SD is the probability-weighted combination of:
  //   evergreen: cv_evergreen × evergreen_biomass × p_evergreen
  //   deciduous: cv_deciduous × deciduous_biomass × p_deciduous
  // (shrub and herbaceous contribute 0)
  //
  // This is the marginal allometric uncertainty per pixel; combining it
  // with the CHM PI uncertainty would require assumptions about independence
  // and propagation through the nonlinear allometric function. We export
  // it as a separate band for transparency.
 
  var cvEvergreen = Posterior.computeAllometricCV(heightMean, Posterior.COV_EVERGREEN);
  var cvDeciduous = Posterior.computeAllometricCV(heightMean, Posterior.COV_DECIDUOUS);
 
  // Per-class biomass at mean height (needed for CV → SD conversion)
  var bioEvergreen = predictTreeBiomass(heightMean, MODELS.evergreen);
  var bioDeciduous = predictTreeBiomass(heightMean, MODELS.deciduous);
 
  // Per-class allometric SD (Mg/ha)
  var sdEvergreen = bioEvergreen.multiply(cvEvergreen);
  var sdDeciduous = bioDeciduous.multiply(cvDeciduous);
 
  // Probability-weighted allometric SD across the two tree classes.
  // Variances add (assuming independence), so weighted SD comes from
  // the square root of the probability-weighted variance sum.
  var pEv = vegClassImage.select(B.VEG.probEvergreen);
  var pDe = vegClassImage.select(B.VEG.probDeciduous);
 
  var weightedVariance = sdEvergreen.pow(2).multiply(pEv)
      .add(sdDeciduous.pow(2).multiply(pDe));
  var allometricSD = weightedVariance.sqrt().rename('allometric_sd');
 
  // CV stored as × 1000 (3 decimal precision in Int16)
  var allometricCV = allometricSD.divide(biomassWeighted.max(0.01))
      .multiply(1000)
      .rename('allometric_cv');
 
  // --- Assemble and export as Int16 × 10 ---
  // Multiply by BIOMASS_SCALE before casting to Int16. This gives 0.1 Mg/ha
  // precision with an Int16 range of -3276.8 to 3276.7 Mg/ha — more than
  // sufficient for riparian biomass (typically 0-300 Mg/ha).
  // asset_loader.loadBiomass() divides by BIOMASS_SCALE on load.
  
  // the allometric cv band uses its own x 1000 scaling because CV is dimensionless
  // and small. It is assembled separately and combined with the rescaled biomass bands
  var BIOMASS_SCALE = Loader.BIOMASS_SCALE;
 
  var biomassBands = biomassMean
      .addBands(biomassLower)
      .addBands(biomassUpper)
      .addBands(biomassWeighted)
      .addBands(biomassWeightedLower)
      .addBands(biomassWeightedUpper)
      .addBands(piWidth)
      .addBands(allometricSD)
      .multiply(BIOMASS_SCALE)
      .toInt16();
 
   // allometric_cv is already scaled × 1000, just cast
  var cvBand = allometricCV.toInt16();
 
  var output = biomassBands.addBands(cvBand).clip(GEOMETRY);
 
  Export.image.toAsset({
    image: output,
    description: 'Biomass_' + year,
    assetId: Loader.PATHS.biomassFolder + '/' + Loader.PATHS.biomassPrefix + year,
    region: GEOMETRY,
    scale: SCALE,
    crs: CRS,
    maxPixels: 1e13
  });
 
});
 
