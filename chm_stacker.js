/**
 * chm_stacker.js
 * ---------------
 * Produces a per-year aligned training stack: predictor bands + a CanopyHeight
 * band aggregated from native LiDAR resolution (~1 m) to 10 m.
 *
 * The CHM footprint defines where training happens — only pixels with valid
 * CHM data are retained. Predictors are unmasked (filled with 0) to ensure
 * consistent band presence, then the CHM mask controls which pixels enter
 * training.
 *
 * CRS:
 *   Target is EPSG:4326 (WGS84) at 10 m nominal scale.
 *
 *   Previous version used EPSG:6350 (NAD83 CONUS Albers), which produced
 *   angled pixels in the GEE map viewer and forced costly reprojection on
 *   every downstream operation. EPSG:4326 aligns natively with the GEE
 *   basemap and Sentinel-2 display grid.
 *
 *   At Idaho/Montana latitudes (~43-48°N), a 10 m nominal scale in EPSG:4326
 *   produces slightly rectangular pixels (~10 m N-S × ~7 m E-W due to
 *   meridian convergence). This is negligible for the ML pipeline (models
 *   are scale-invariant) and for area calculations (ee.Image.pixelArea()
 *   uses geodetic area regardless of CRS).
 *
 * Reprojection strategy:
 *   The predictor stack (from predictor_stack_module) is NOT reprojected
 *   here. Each predictor band retains its native CRS, and GEE handles
 *   alignment lazily when the consumer samples or exports.
 *
 *   The CHM IS explicitly reprojected via reduceResolution + reproject.
 *   This is structurally necessary: we need to aggregate 1 m LiDAR pixels
 *   into 10 m means, and that requires a defined target grid. Without
 *   setDefaultProjection → reduceResolution → reproject, GEE doesn't know
 *   the native resolution and can't compute the aggregation ratio.
 *
 *   The result is a training stack where the CHM band is in EPSG:4326 at
 *   10 m and the predictor bands are in mixed native CRS. When the main
 *   script calls sample() with scale: 10, GEE aligns everything to the
 *   same 10 m grid at sample locations. This is the GEE-idiomatic way
 *   to handle multi-resolution inputs.
 *
 * Yearly masks:
 *   NOT applied during training. The CHM footprint alone controls which
 *   pixels have training labels. Yearly masks are only applied during
 *   prediction (in the main timeseries script) to constrain output to
 *   the valid riparian/valley bottom extent for each year.
 */
 
var PredictorPrep = require('users/eriknatijohnson/MRRMaidCarbon:finalstack/predictor_stack_module');
 
var TARGET_SCALE = 10;
var TARGET_CRS = 'EPSG:4326';
 
/**
 * Builds the aligned training stack for a single year.
 *
 * @param {ee.Number|number} YEAR - Calendar year (2017-2022 for LiDAR years).
 * @param {ee.Geometry} geometry - Area of interest.
 * @param {string} CHM_COLLECTION_ID - Asset ID of the CHM ImageCollection.
 * @returns {ee.Image} Multi-band training image: 24 predictor bands +
 *   CanopyHeight. Masked to the CHM footprint. CHM band is in EPSG:4326
 *   at 10 m; predictor bands are in mixed native CRS (aligned lazily
 *   at sample time).
 */
exports.getTrainingStack = function(YEAR, geometry, CHM_COLLECTION_ID) {
 
  // -------------------------------------------------------------------------
  // 1. Load predictor stack
  // -------------------------------------------------------------------------
  // Returns 24-band image in mixed native CRS (no explicit reprojection).
  // Bands include spectral, phenology, ecoregion, terrain, and climate.
  var predictorStack = ee.Image(PredictorPrep.getSentinelPredictorStack(YEAR, geometry));
 
  // -------------------------------------------------------------------------
  // 2. Load CHM for this year
  // -------------------------------------------------------------------------
  // LiDAR CHMs are organized as one image per year in the CHM collection.
  // If no CHM exists for this year (shouldn't happen for 2017-2022 but
  // defensive), create a fully masked fallback that produces zero valid
  // training pixels for that year without breaking the pipeline.
  var start = ee.Date.fromYMD(YEAR, 1, 1);
  var end = start.advance(1, 'year');
 
  var chmCollection = ee.ImageCollection(CHM_COLLECTION_ID).filterDate(start, end);
 
  var chmAsset = ee.Image(
    ee.Algorithms.If(
      chmCollection.size().gt(0),
      chmCollection.first(),
      ee.Image.constant(0).rename('CanopyHeight').updateMask(ee.Image(0))
    )
  ).rename('CanopyHeight');
 
  // -------------------------------------------------------------------------
  // 3. Aggregate CHM from native resolution to 10 m
  // -------------------------------------------------------------------------
  // The LiDAR CHMs are typically at ~1 m. We aggregate to 10 m via mean
  // so each training pixel represents the average canopy height over the
  // 10 m footprint — matching what Sentinel-2 "sees" at that pixel.
  //
  // The aggregation chain:
  //   setDefaultProjection → tells GEE the CHM's native CRS and scale,
  //     so reduceResolution knows the source pixel grid.
  //   reduceResolution → computes the mean of native pixels falling within
  //     each target pixel. bestEffort: true allows GEE to adjust if the
  //     native scale doesn't divide evenly into the target scale.
  //   reproject → places the aggregated result onto the EPSG:4326 grid
  //     at 10 m nominal scale.
  //
  // If reproject were called first, GEE would request
  // CHM data in the target grid (potentially with wrong resolution), then
  // reduceResolution would have nothing useful to aggregate.
  var chmNativeProj = chmAsset.projection();
 
  var chmMean10m = chmAsset
      .setDefaultProjection(chmNativeProj)
      .reduceResolution({
        reducer: ee.Reducer.mean(),
        bestEffort: true
      })
      .reproject({crs: TARGET_CRS, scale: TARGET_SCALE})
      .rename('CanopyHeight');
 
  // -------------------------------------------------------------------------
  // 4. Assemble training stack
  // -------------------------------------------------------------------------
  // unmask(0) fills predictor gaps with zeros. This ensures all 24 predictor
  // bands have values everywhere, so the CHM mask alone determines which
  // pixels enter training.
  //
  // Trade-off: unmask introduces artificial zeros where S2 data was missing.
  // Acceptable because (a) two-month median composites rarely have gaps at
  // 10 m in the Intermountain West, and (b) the CHM footprint is typically
  // smaller than S2 coverage, so most unmasked zeros fall outside the
  // training area anyway.
  var unmaskedPredictors = predictorStack.unmask(0);
 
  // Combine predictors + CHM, then mask to CHM footprint.
  // Only pixels with valid LiDAR data are retained for training.
  var finalTrainingStack = unmaskedPredictors
      .addBands(chmMean10m)
      .updateMask(chmMean10m.mask())
      .toFloat();
 
  return finalTrainingStack;
};
