/**
 * asset_loader.js
 * ----------------
 * Shared utility module for loading pre-computed raster assets with
 * consistent rescaling, masking, and band naming.
 *
 * WHY THIS EXISTS:
 *   The v2 CHM assets are stored as Int16 × 100 to halve storage costs.
 *   Without a shared loader, every script that touches these assets
 *   (biomass precompute, temporal change, the app) would need its own
 *   division-by-100 step. This module centralizes that logic.
 *
 *   It also centralizes yearly mask loading and application, asset path
 *   construction, and band name constants.
 *
 * WHAT IT PROVIDES:
 *   loadCHM(year)          → CHM image in meters (Float32), all 7 bands
 *   loadBiomass(year)      → Annual biomass image (Float32), all bands
 *   loadVegClassification() → Vegetation classification image
 *   loadYearlyMask(year)   → Cleaned binary mask (1 = valid)
 *   applyYearlyMask(image, year) → Mask an image to valid riparian pixels
 *   getUnionMask()         → Union of all yearly masks
 *   BAND_NAMES             → Constants for band name references
 *   PATHS                  → Asset path constants
 *
 * CRS:
 *   This module is CRS-agnostic. It loads assets as-is and does not
 *   reproject. The consumer controls CRS at sample/export time.
 *
 * USAGE:
 *   var Loader = require('users/eriknatijohnson/MRRMaidCarbon:finalstack/asset_loader');
 *   var chm2020 = Loader.loadCHM(2020);
 *   var height = chm2020.select(Loader.BAND_NAMES.CHM.mean);  // in meters
 */
 
// =============================================================================
// ASSET PATHS
// =============================================================================
// Primary project for raw inputs (study area, training points, masks)
var SOURCE = 'projects/ee-eriknatijohnson/assets';
 
// Output project for all v2 products
var CARBON = 'projects/mrrmaid-carbon/assets';
 
var PATHS = {
  // --- Source project (read-only inputs) ---
  studyArea:     SOURCE + '/studyarea',
  trainingPoints: SOURCE + '/multiclass_points',
  yearlyMasks:   SOURCE + '/YEARLYMASKS',
 
  // --- Carbon project (v2 products) ---
  chmFolder:     CARBON + '/Pred_CHM_v2',
  chmPrefix:     'Predicted_CHM_v2_',
  vegClass:      CARBON + '/Veg_Classification_v3',
  globalValidMask: CARBON + '/global_valid_mask',
  biomassFolder: CARBON + '/Annual_Biomass',
  biomassPrefix: 'Biomass_',
  uncertaintyFolder: CARBON + '/Combined_Uncertainty',
  uncertaintyPrefix: 'CombinedUncertainty_',
  yoyChangePrefix: 'YoY_Change_',
  changeFolder:  CARBON + '/Biomass_Change'
};
 
exports.PATHS = PATHS;
 
// =============================================================================
// SCALE FACTORS
// =============================================================================
// CHM outputs are Int16 × 100. Divide by this to get meters.
// Biomass outputs are Int16 × 10. Divide by this to get Mg/ha.
// sd_log_total and rel_uncertainty are stored Int16 × 1000 (dimensionless,
// small). biomass_lower_total / biomass_upper_total use BIOMASS_SCALE (×10).
 
var CHM_SCALE = 100;
exports.CHM_SCALE = CHM_SCALE;
 
var BIOMASS_SCALE = 10;
exports.BIOMASS_SCALE = BIOMASS_SCALE;

var UNCERTAINTY_SCALE = 1000;
exports.UNCERTAINTY_SCALE = UNCERTAINTY_SCALE;
 
 
// =============================================================================
// BAND NAME CONSTANTS
// =============================================================================
// Centralized band name references. 
 
var BAND_NAMES = {
  CHM: {
    mean:    'Predicted_CHM',
    std:     'Ensemble_StdDev',
    piLower: 'PI_Lower',
    piUpper: 'PI_Upper',
    piWidth: 'PI_Width',
    rf:      'RF_Prediction',
    gbt:     'GBT_Prediction'
  },
  BIOMASS: {
    mean:           'biomass_mean',
    lower:          'biomass_lower',
    upper:          'biomass_upper',
    weighted:       'biomass_weighted',
    weightedLower:  'biomass_weighted_lower',
    weightedUpper:  'biomass_weighted_upper',
    piWidth:        'prediction_interval_width',
    allometricSD:   'allometric_sd',
    allometricCV:   'allometric_cv'
  },
  VEG: {
    hardClass:      'class',
    probEvergreen:  'prob_evergreen',
    probDeciduous:  'prob_deciduous',
    probHerbaceous: 'prob_herbaceous',
    probShrub:      'prob_shrub',
    setSize:        'set_size',
    inSetEvergreen: 'in_set_evergreen',
    inSetDeciduous: 'in_set_deciduous',
    inSetHerbaceous:'in_set_herbaceous',
    inSetShrub:     'in_set_shrub',
    maxProb:        'max_prob',
    entropyNorm:    'entropy_norm'
  },
  CHANGE: {
    absChange:   'abs_change_Mg_ha',
    relChange:   'rel_change',
    trendSlope:  'trend_Mg_ha_per_year',
    trendR2:     'trend_r2'
  },
  UNCERTAINTY: {
    sdLogTotal:     'sd_log_total',
    lowerTotal:     'biomass_lower_total',
    upperTotal:     'biomass_upper_total',
    relUncertainty: 'rel_uncertainty'
   },
};
 
exports.BAND_NAMES = BAND_NAMES;
 
 
// =============================================================================
// CLASS METADATA
// =============================================================================
 
exports.CLASS_NAMES = {0: 'Evergreen', 1: 'Deciduous', 2: 'Herbaceous', 3: 'Shrub'};
exports.CLASS_PALETTE = ['#1b7837', '#fdb863', '#d9d9d9', '#762a83'];
 
// Years covered by the pipeline
exports.YEARS = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025];
 
 
// =============================================================================
// CHM LOADING
// =============================================================================
 
/**
 * Loads a v2 CHM asset for the given year and rescales from Int16 × 100
 * to Float32 meters. Returns all 7 bands.
 *
 * @param {number} year - Calendar year (2017-2025).
 * @returns {ee.Image} 7-band Float32 image in meters:
 *   Predicted_CHM, Ensemble_StdDev, PI_Lower, PI_Upper, PI_Width,
 *   RF_Prediction, GBT_Prediction.
 */
exports.loadCHM = function(year) {
  var assetId = PATHS.chmFolder + '/' + PATHS.chmPrefix + year;
  var raw = ee.Image(assetId);
 
  // All 7 bands have the same scale factor (Int16 × 100)
  return raw.toFloat().divide(CHM_SCALE);
};
 
 
// =============================================================================
// BIOMASS LOADING
// =============================================================================
 
/**
 * Loads a pre-computed annual biomass asset and rescales from Int16 × 10
 * to Float32 Mg/ha.
 * 
 * allometric_cv band is stored as Int16 x 1000 because CV is dimensionless 
 * and small, needing finer precision than the biomass bands.
 *
 * @param {number} year - Calendar year (2017-2025).
 * @returns {ee.Image} Multi-band biomass image (Mg/ha, Float32):
 *   biomass_mean, biomass_lower, biomass_upper,
 *   biomass_weighted, biomass_weighted_lower, biomass_weighted_upper,
 *   prediction_interval_width.
 */
exports.loadBiomass = function(year) {
  var assetId = PATHS.biomassFolder + '/' + PATHS.biomassPrefix + year;
  var raw = ee.Image(assetId);
 
  // Biomass bands: divide by BIOMASS_SCALE (10) to get Mg/ha
  var biomassBands = raw.select([
    'biomass_mean', 'biomass_lower', 'biomass_upper',
    'biomass_weighted', 'biomass_weighted_lower', 'biomass_weighted_upper',
    'prediction_interval_width', 'allometric_sd'
  ]).toFloat().divide(BIOMASS_SCALE);
 
  // CV band: divide by 1000 to get dimensionless ratio
  var cvBand = raw.select('allometric_cv').toFloat().divide(1000);
 
  return biomassBands.addBands(cvBand);
};


// =============================================================================
// COMBINED UNCERTAINTY LOADING
// =============================================================================
/**
 * Loads the precomputed combined-uncertainty asset for a year and rescales
 * to Float32.
 *
 *   biomass_lower_total, biomass_upper_total → Mg/ha  (÷ BIOMASS_SCALE)
 *   sd_log_total          → dimensionless log-scale SD (÷ UNCERTAINTY_SCALE)
 *   rel_uncertainty       → dimensionless half-width/mean (÷ UNCERTAINTY_SCALE)
 *
 * The combined interval reflects CHM (conformal) height uncertainty and
 * allometric parameter uncertainty. It excludes TPA-scaling uncertainty, and
 * in shrub/herbaceous-dominated pixels reflects height uncertainty only.
 *
 * @param {number} year - Calendar year (2017-2025).
 * @returns {ee.Image} 4-band Float32 image.
 */
exports.loadCombinedUncertainty = function(year) {
  var assetId = exports.PATHS.uncertaintyFolder + '/' +
                exports.PATHS.uncertaintyPrefix + year;
  var raw = ee.Image(assetId);
 
  var mgha = raw.select([
    'biomass_lower_total', 'biomass_upper_total'
  ]).toFloat().divide(exports.BIOMASS_SCALE);
 
  var dimensionless = raw.select([
    'sd_log_total', 'rel_uncertainty'
  ]).toFloat().divide(UNCERTAINTY_SCALE);
 
  return mgha.addBands(dimensionless);
};
 
/**
 * Convenience: combined uncertainty with the year-specific riparian mask
 * applied (mirrors loadMaskedBiomass).
 *
 * @param {number} year - Calendar year.
 * @returns {ee.Image} Masked combined-uncertainty image.
 */
exports.loadMaskedCombinedUncertainty = function(year) {
  var unc = exports.loadCombinedUncertainty(year);
  return exports.applyYearlyMask(unc, year);
};
  
/**
 * Loads a pre-computed year-over-year biomass change asset.
 *
 * Each asset is masked to pixels valid in both years of the pair — a
 * weaker constraint than the global valid mask used for time series and
 * trend products, because pairwise differences are meaningful wherever
 * both years have data.
 *
 * @param {number} yearFrom - Earlier year of the pair.
 * @param {number} yearTo   - Later year (must equal yearFrom + 1).
 * @returns {ee.Image} 2-band Float32 image:
 *   abs_change_Mg_ha — yearTo biomass minus yearFrom biomass (Mg/ha)
 *   rel_change       — abs_change / yearFrom biomass, masked where prior
 *                       year biomass is below 1 Mg/ha
 */
exports.loadYoYChange = function(yearFrom, yearTo) {
  if (yearTo !== yearFrom + 1) {
    throw new Error('loadYoYChange requires consecutive years; got ' +
                    yearFrom + ' and ' + yearTo);
  }
  var assetId = PATHS.changeFolder + '/' + PATHS.yoyChangePrefix +
                yearFrom + '_to_' + yearTo;
  return ee.Image(assetId);
};
// =============================================================================
// VEGETATION CLASSIFICATION LOADING
// =============================================================================
 
/**
 * Loads the v2 vegetation classification asset.
 *
 * @returns {ee.Image} Multi-band classification image:
 *   class, prob_*, set_size, in_set_*, max_prob, entropy_norm.
 */
exports.loadVegClassification = function() {
  return ee.Image(PATHS.vegClass);
};
 
 
// =============================================================================
// YEARLY MASK LOADING
// =============================================================================
 
/**
 * Loads and cleans the yearly riparian/valley bottom mask for a given year.
 * Returns a binary mask: 1 = valid riparian pixel, 0/masked = upland.
 *
 * Raw mask encoding: 1 = valid, 65535 = also valid (older encoding quirk).
 * Everything else is masked.
 *
 * @param {number} year - Calendar year.
 * @returns {ee.Image} Binary mask (1 = valid), self-masked.
 */
exports.loadYearlyMask = function(year) {
  var maskId = year + 'mask';
  var raw = ee.ImageCollection(PATHS.yearlyMasks)
      .filter(ee.Filter.eq('system:index', maskId))
      .first();
 
  raw = ee.Image(raw).select(0);
  return raw.eq(1).or(raw.eq(65535)).rename('mask').toUint8().selfMask();
};
 
/**
 * Applies the year-specific riparian mask to an image. Pixels outside
 * the valid riparian extent for that year are masked to nodata.
 *
 * @param {ee.Image} image - Image to mask.
 * @param {number} year - Calendar year.
 * @returns {ee.Image} Masked image.
 */
exports.applyYearlyMask = function(image, year) {
  var mask = exports.loadYearlyMask(year);
  return image.updateMask(mask);
};
 
/**
 * Computes the union of all yearly masks — any pixel valid in ANY year.
 * Used by the vegetation classifier (which covers the maximum extent).
 *
 * @param {ee.Geometry} geometry - Area to clip to.
 * @returns {ee.Image} Binary mask (1 = valid in at least one year).
 */
exports.getUnionMask = function(geometry) {
  var years = exports.YEARS;
  var maskImages = years.map(function(year) {
    return exports.loadYearlyMask(year);
  });
 
  return ee.ImageCollection.fromImages(maskImages)
      .reduce(ee.Reducer.anyNonZero())
      .rename('mask')
      .clip(geometry);
};
 
 
// =============================================================================
// CONVENIENCE: LOAD MASKED BIOMASS
// =============================================================================
 
/**
 * Loads biomass for a given year with the year-specific riparian mask
 * applied. This is the most common loading pattern for the app and
 * analysis scripts.
 *
 * @param {number} year - Calendar year.
 * @returns {ee.Image} Masked biomass image (Mg/ha).
 */
exports.loadMaskedBiomass = function(year) {
  var biomass = exports.loadBiomass(year);
  return exports.applyYearlyMask(biomass, year);
};
 
/**
 * Loads CHM for a given year with the year-specific riparian mask applied.
 *
 * @param {number} year - Calendar year.
 * @returns {ee.Image} Masked CHM image (meters).
 */
exports.loadMaskedCHM = function(year) {
  var chm = exports.loadCHM(year);
  return exports.applyYearlyMask(chm, year);
};



// ===========================================================
// global valid mask
//============================================================

exports.loadGlobalValidMask = function() {
  return ee.Image(PATHS.globalValidMask);
};
