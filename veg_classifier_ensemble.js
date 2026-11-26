/**
 * veg_classifier_ensemble.js
 * --------------------------
 * Four-class vegetation classification with soft classification and
 * conformal prediction sets for uncertainty quantification.
 *
 * Classes:
 *   0 = Evergreen (conifer-dominated, primarily ponderosa pine)
 *   1 = Deciduous (broadleaf-dominated, primarily cottonwood)
 *   2 = Herbaceous (grassland, meadow)
 *   3 = Shrub (primarily willow)
 *
 * ARCHITECTURE (v2):
 *   One tuned RF + one tuned GBT, consistent with the CHM prediction
 *   pipeline. Replaces the previous 10-model ensemble (5 RF + 5 GBT)
 *   with hand-picked hyperparameters. The 2-model architecture makes
 *   LOOCV dramatically cheaper (621 × 2 = 1,242 model fits vs 6,210)
 *   and aligns the ensemble rationale across the full pipeline.
 *
 * PREDICTOR STACK (27 bands):
 *   The classifier uses a different predictor set than the CHM pipeline.
 *   Monthly NDVI (12 bands) captures phenological signatures that separate
 *   deciduous from evergreen and shrub from herbaceous — these seasonal
 *   curves are the single most important discriminator. The CHM itself
 *   is included as a predictor because height separates vegetation
 *   structural classes.
 *
 *   CHM + CHM_StdDev      2 bands   from v2 ensemble mosaic (Int16×100)
 *   Monthly NDVI           12 bands  multi-year median per month (2017-2022)
 *   June-July spectral     6 bands   multi-year median (B2, B3, B4, B5, B8, B11)
 *   Terrain                4 bands   Elevation, Slope, Aspect_sin, Aspect_cos
 *   Climate                3 bands   tmean_annual, ppt_annual, temp_range
 *
 * CONFORMAL PREDICTION:
 *   Following Valle et al. (2023), conformal prediction sets are computed
 *   with a 90% coverage target. The nonconformity score for classification
 *   is 1 − P(true class) — lower scores mean the model was more confident
 *   on the correct class. The conformal threshold is calibrated from
 *   LOOCV held-out predictions, so no data leakage.
 *
 *   A pixel's prediction set includes all classes whose ensemble
 *   probability ≥ threshold. Set size of 1 = confident classification;
 *   set size > 1 = genuine ambiguity between classes.
 *
 * LOOCV:
 *   The evaluation and conformal calibration run in a separate companion
 *   script (loocv_evaluation.js) that holds out each of the ~621 training
 *   samples, retrains the 2-model ensemble on the remaining 620, and
 *   predicts the held-out sample. With only 2 models per iteration,
 *   the computation (621 × 2 = 1,242 fits) completes reliably as a
 *   batch export task. The LOOCV CSV is then processed in R to compute
 *   the confusion matrix, accuracy metrics, and conformal threshold.
 *   That threshold is entered here in CONFORMAL_THRESHOLD.
 *
 * CRS:
 *   EPSG:4326 at 10 m nominal scale. Matches the v2 CHM pipeline.
 *
 * TERMINOLOGY:
 *   Conformal prediction for classification produces "prediction sets"
 *   (not intervals). A prediction set is a subset of classes guaranteed
 *   to contain the true class with probability ≥ 1 − α.
 */
 
// =============================================================================
// 1. MODULE IMPORTS
// =============================================================================
// Terrain and climate predictors come from the shared module to avoid
// duplicating code. The classifier builds its own spectral/NDVI stack
// because it uses multi-year monthly NDVI (12 bands) rather than the
// CHM pipeline's seasonal NDVI/EVI (6 bands).
 
var PredictorPrep = require('users/eriknatijohnson/MRRMaidCarbon:finalstack/predictor_stack_module');
 
// =============================================================================
// 2. CONFIGURATION
// =============================================================================
 
var studyarea = ee.FeatureCollection('projects/ee-eriknatijohnson/assets/studyarea');
var trainingPoints = ee.FeatureCollection('projects/ee-eriknatijohnson/assets/multiclass_points');
// Classes: 0 evergreen, 1 deciduous, 2 herbaceous, 3 shrub
 
var geometry = studyarea.geometry().bounds();
 
// Years used for multi-year compositing. Compositing across all LiDAR years
// produces stable phenological signatures and spectral values that aren't
// biased by any single year's weather or phenological timing.
var COMPOSITE_YEARS = [2017, 2018, 2019, 2020, 2021, 2022];
 
// CHM mosaic priority: 2020 first (best overall LiDAR coverage), then
// gap-filled by other years ranked by mean ensemble agreement.
var CHM_PRIORITY = [2020, 2018, 2019, 2017, 2021, 2022];
 
var SCALE = 10;
var CRS = 'EPSG:4326';
 
// --- Conformal prediction ---
// 90% coverage target following Valle et al. (2023).
// 90% rather than 95% is standard for classification conformal prediction
// because the discrete nature of classes makes 95% overly conservative
// (most pixels would have 2+ classes in the prediction set).
var CONFORMAL_COVERAGE = 0.90;
 
// Conformal threshold from LOOCV held-out predictions.
// This value is computed by loocv_evaluation.js and the companion R script.
var CONFORMAL_THRESHOLD = 0.217617;  // from LOOCV results
 
// --- CHM v2 assets ---
// The v2 CHM assets store height as Int16 × 100 (divide by 100 for meters).
var CHM_FOLDER = 'projects/mrrmaid-carbon/assets/Pred_CHM_v2';
var CHM_PREFIX = 'Predicted_CHM_v2_';
var CHM_SCALE_FACTOR = 100;
 
var YEARLY_MASK_COLLECTION_ID = 'projects/ee-eriknatijohnson/assets/YEARLYMASKS';
 
// --- Cloud Score+ threshold ---
// Using 'cs' band (clear sky probability, 0-1). Pixels below this value
// are masked. 0.65 is slightly more conservative than the CHM pipeline's
// 0.60 because classification is more sensitive to residual haze than
// regression 
var CLOUD_SCORE_THRESH = 0.65;
 
// --- Output paths ---
var DRIVE_FOLDER = 'GEE_VegClassifier_v2_Outputs';
var ASSET_FOLDER = 'projects/mrrmaid-carbon/assets';
 
// --- Class metadata ---
var CLASS_NAMES = {0: 'Evergreen', 1: 'Deciduous', 2: 'Herbaceous', 3: 'Shrub'};
var CLASS_PALETTE = ['#1b7837', '#fdb863', '#d9d9d9', '#762a83'];
 
 
// =============================================================================
// 2B. MODEL HYPERPARAMETERS
// =============================================================================
// Selected via tuning (same framework as the CHM pipeline).
// best RF: Accuracy: 0.8435  Kappa: 0.7836 
// best GBT: Accuracy: 0.8522  Kappa: 0.796
 
var SELECTED_RF = {
  numberOfTrees: 300,
  variablesPerSplit: 5,
  minLeafPopulation: 1,
  bagFraction: 0.5,
  maxNodes: null,
  seed: 42
};
 
var SELECTED_GBT = {
  numberOfTrees: 200,
  shrinkage: 0.05,
  maxNodes: 10,
  samplingRate: 0.8,
  loss: 'LeastSquares',
  seed: 43
};
 
 
// =============================================================================
// 3. SENTINEL-2 CLOUD MASKING
// =============================================================================
 
var csPlus = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED');
 
function maskWithCloudScorePlus(s2Collection, startDate, endDate, bounds) {
  var csPlusFiltered = csPlus
    .filterDate(startDate, endDate)
    .filterBounds(bounds);
 
  var s2Filtered = s2Collection
    .filterDate(startDate, endDate)
    .filterBounds(bounds);
 
  var linked = s2Filtered.linkCollection(csPlusFiltered, ['cs']);
 
  return linked.map(function(img) {
    var csScore = img.select('cs');
    return img.updateMask(csScore.gte(CLOUD_SCORE_THRESH));
  });
}
 
 
// =============================================================================
// 4. MULTI-YEAR MONTHLY NDVI (12 bands)
// =============================================================================
// The monthly NDVI stack is the most important predictor group for
// vegetation classification. Each band is the median NDVI for that
// calendar month across all composite years (2017-2022). This produces
// a stable phenological curve that captures:
//   - Evergreen: flat NDVI year-round (high winter, moderate summer)
//   - Deciduous: strong seasonal swing (low winter, high summer)
//   - Herbaceous: early green-up, mid-summer senescence
//   - Shrub: late green-up, sustained through early fall
//
// No explicit reprojection. S2 is natively ~10 m in EPSG:32611/32612
// (UTM zones covering the study area). GEE resamples lazily at sample time.
 
function getMultiYearMonthlyNDVI(years, geom) {
  var months = ee.List.sequence(1, 12);
  var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED');
 
  var monthlyBands = months.map(function(m) {
    m = ee.Number(m).toInt();
 
    // For each month, compute the median NDVI across all years.
    // Two-level compositing: within-year median (handles cloud gaps),
    // then cross-year median (handles phenological outliers).
    var monthImages = ee.ImageCollection.fromImages(
      years.map(function(year) {
        year = ee.Number(year).toInt();
        var start = ee.Date.fromYMD(year, m, 1);
        var end = start.advance(1, 'month');
 
        var masked = maskWithCloudScorePlus(s2, start, end, geom);
 
        var ndvi = masked.map(function(img) {
          return img.normalizedDifference(['B8', 'B4']).rename('NDVI');
        }).median();
 
        return ndvi;
      })
    );
 
    var monthMedian = monthImages.median();
    var monthStr = ee.String('NDVI_').cat(ee.Number(m).format('%02d'));
    return monthMedian.rename(monthStr);
  });
 
  return ee.ImageCollection.fromImages(monthlyBands).toBands()
    .regexpRename('^[0-9]+_', '');
}
 
 
// =============================================================================
// 5. MULTI-YEAR JUNE-JULY SPECTRAL COMPOSITE (6 bands)
// =============================================================================
// Peak growing season spectral signature. B2 (blue), B3 (green), B4 (red),
// B5 (red edge), B8 (NIR), B11 (SWIR). These capture structural and
// moisture differences between classes that NDVI alone misses — e.g.,
// SWIR separates dry herbaceous from moist shrub canopies.
 
function getMultiYearJuneJulyComposite(years, geom) {
  var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED');
  var bands = ['B2', 'B3', 'B4', 'B5', 'B8', 'B11'];
 
  var yearlyComposites = ee.ImageCollection.fromImages(
    years.map(function(year) {
      year = ee.Number(year).toInt();
      var start = ee.Date.fromYMD(year, 6, 1);
      var end = ee.Date.fromYMD(year, 7, 31).advance(1, 'day');
      var masked = maskWithCloudScorePlus(s2, start, end, geom);
      return masked.select(bands).median();
    })
  );
 
  return yearlyComposites.median()
    .rename(['B2_JJ', 'B3_JJ', 'B4_JJ', 'B5_JJ', 'B8_JJ', 'B11_JJ']);
}
 
 
// =============================================================================
// 6. CHM MOSAIC
// =============================================================================
// The CHM is a predictor for the classifier (not the output). Height is
// the primary separator between herbaceous (short), shrub (medium), and
// tree classes (tall), and helps distinguish deciduous from evergreen
// in combination with phenological signal.
//
// Priority mosaic: 2020 first (best LiDAR coverage), then gap-filled by
// other years. Reversed before mosaic() so the highest-priority year's
// pixels end up on top.
//
// v2 CHM assets are Int16 × 100. We divide by CHM_SCALE_FACTOR to get
// meters before adding to the predictor stack.
 
function getCHMMosaic(priorityYears, geom) {
  var reversedPriority = priorityYears.slice().reverse();
 
  var chmImages = reversedPriority.map(function(year) {
    var assetId = CHM_FOLDER + '/' + CHM_PREFIX + year;
    var img = ee.Image(assetId);
 
    // Rescale from Int16 × 100 to meters
    var chm = img.select('Predicted_CHM')
        .divide(CHM_SCALE_FACTOR)
        .rename('CHM')
        .toFloat();
    var stddev = img.select('Ensemble_StdDev')
        .divide(CHM_SCALE_FACTOR)
        .rename('CHM_StdDev')
        .toFloat();
 
    return chm.addBands(stddev);
  });
 
  return ee.ImageCollection.fromImages(chmImages).mosaic().clip(geom);
}
 
 
// =============================================================================
// 7. UNION MASK
// =============================================================================
// The vegetation classification covers the union of all yearly masks —
// any pixel that was valid in ANY year gets classified. This ensures
// maximum spatial coverage. Individual year biomass products then apply
// year-specific masks to constrain where estimates are reported.
//
// Mask values: 1 = valid, 65535 = valid (older encoding). Everything
// else is masked.
 
function getUnionMask(years, geom) {
  var yearlyMasks = ee.ImageCollection(YEARLY_MASK_COLLECTION_ID);
 
  var maskImages = years.map(function(year) {
    var maskId = year + 'mask';
    var mask = yearlyMasks
      .filter(ee.Filter.eq('system:index', maskId))
      .first();
 
    mask = ee.Image(mask).select(0);
    return mask.eq(1).or(mask.eq(65535)).selfMask();
  });
 
  var unionMask = ee.ImageCollection.fromImages(maskImages)
    .reduce(ee.Reducer.anyNonZero())
    .rename('mask');
 
  return unionMask.clip(geom);
}
 
 
// =============================================================================
// 8. BUILD PREDICTOR STACK (27 bands)
// =============================================================================
 
 
var yearsList = ee.List(COMPOSITE_YEARS);
 
var ndviStack = getMultiYearMonthlyNDVI(yearsList, geometry);
var jjComposite = getMultiYearJuneJulyComposite(yearsList, geometry);
var chmMosaic = getCHMMosaic(CHM_PRIORITY, geometry);
 
// Terrain and climate from shared module — avoids duplicating code.
// These functions return images in their native CRS (no reprojection).
var terrainStack = PredictorPrep.getTerrainStack(geometry).toFloat();
var climateStack = PredictorPrep.getClimateStack(geometry).toFloat();
 
var unionMask = getUnionMask(COMPOSITE_YEARS, geometry);
 
// Assemble. No explicit reproject — GEE aligns bands lazily at sample/export time.
var predictorStack = chmMosaic
  .addBands(ndviStack)
  .addBands(jjComposite)
  .addBands(terrainStack)
  .addBands(climateStack);
 
predictorStack = predictorStack.updateMask(unionMask);
 
var inputBands = predictorStack.bandNames();
 
 
// =============================================================================
// 9. SAMPLE TRAINING POINTS
// =============================================================================
// Training points are pre-curated from multiple sources, filtered
// and verified for the study area. Each point has a 'class' property (0-3).
// We sample the full predictor stack at each point and drop any with null
// values (edge pixels, projection artifacts).
 
var samples = predictorStack.sampleRegions({
  collection: trainingPoints,
  properties: ['class'],
  scale: SCALE,
  geometries: true
}).filter(ee.Filter.notNull(inputBands));
 
var samplesWithId = samples.map(function(f) {
  return f.set('uid', f.id());
});
 
 
 
// =============================================================================
// 10. TRAIN MODELS
// =============================================================================
// Two models: one tuned RF + one tuned GBT, both in MULTIPROBABILITY mode.
// MULTIPROBABILITY returns a 4-element array of class probabilities per pixel,
// which we average across the two models for the ensemble probability.
 
 
// --- RF ---
var rfParams = {
  numberOfTrees: SELECTED_RF.numberOfTrees,
  variablesPerSplit: SELECTED_RF.variablesPerSplit,
  minLeafPopulation: SELECTED_RF.minLeafPopulation,
  bagFraction: SELECTED_RF.bagFraction,
  seed: SELECTED_RF.seed
};
if (SELECTED_RF.maxNodes !== null) {
  rfParams.maxNodes = SELECTED_RF.maxNodes;
}
 
var rfModel = ee.Classifier.smileRandomForest(rfParams)
    .setOutputMode('MULTIPROBABILITY')
    .train({
      features: samplesWithId,
      classProperty: 'class',
      inputProperties: inputBands
    });
 
// --- GBT ---
var gbtModel = ee.Classifier.smileGradientTreeBoost({
  numberOfTrees: SELECTED_GBT.numberOfTrees,
  shrinkage: SELECTED_GBT.shrinkage,
  maxNodes: SELECTED_GBT.maxNodes,
  samplingRate: SELECTED_GBT.samplingRate,
  loss: SELECTED_GBT.loss,
  seed: SELECTED_GBT.seed
})
.setOutputMode('MULTIPROBABILITY')
.train({
  features: samplesWithId,
  classProperty: 'class',
  inputProperties: inputBands
});
 
 
 
// =============================================================================
// 11. COMPUTE ENSEMBLE PROBABILITIES (wall-to-wall)
// =============================================================================
// Classify the entire predictor stack with both models, extract per-class
// probabilities from the array output, and average across the ensemble.
 
var rfProbs = predictorStack.classify(rfModel);
var gbtProbs = predictorStack.classify(gbtModel);
 
// MULTIPROBABILITY output is a 4-element array band named 'classification'.
// Extract each class probability and average across models.
var rfP0 = rfProbs.arrayGet(0);  var gbtP0 = gbtProbs.arrayGet(0);
var rfP1 = rfProbs.arrayGet(1);  var gbtP1 = gbtProbs.arrayGet(1);
var rfP2 = rfProbs.arrayGet(2);  var gbtP2 = gbtProbs.arrayGet(2);
var rfP3 = rfProbs.arrayGet(3);  var gbtP3 = gbtProbs.arrayGet(3);
 
var prob_evergreen  = rfP0.add(gbtP0).divide(2).rename('prob_evergreen');
var prob_deciduous  = rfP1.add(gbtP1).divide(2).rename('prob_deciduous');
var prob_herbaceous = rfP2.add(gbtP2).divide(2).rename('prob_herbaceous');
var prob_shrub      = rfP3.add(gbtP3).divide(2).rename('prob_shrub');
 
var probStack = prob_evergreen
    .addBands(prob_deciduous)
    .addBands(prob_herbaceous)
    .addBands(prob_shrub);
 
 
// =============================================================================
// 12. HARD CLASSIFICATION AND CONFORMAL PREDICTION SETS
// =============================================================================
 
// Hard classification: argmax of ensemble probabilities
var probArray = probStack.toArray();
var hardClass = probArray.arrayArgmax().arrayFlatten([['class']]).toInt();
 
// Conformal prediction sets: include a class if its ensemble probability
// meets or exceeds the calibrated threshold.
// Set size of 1 = confident; 2+ = ambiguous; 0 = very uncertain (rare).
var thresholdImg = ee.Image.constant(CONFORMAL_THRESHOLD);
 
var inSetEvergreen  = prob_evergreen.gte(thresholdImg).rename('in_set_evergreen').toInt();
var inSetDeciduous  = prob_deciduous.gte(thresholdImg).rename('in_set_deciduous').toInt();
var inSetHerbaceous = prob_herbaceous.gte(thresholdImg).rename('in_set_herbaceous').toInt();
var inSetShrub      = prob_shrub.gte(thresholdImg).rename('in_set_shrub').toInt();
 
var setSize = inSetEvergreen.add(inSetDeciduous).add(inSetHerbaceous).add(inSetShrub)
    .rename('set_size');
 
 
// =============================================================================
// 13. UNCERTAINTY METRICS
// =============================================================================
 
// Maximum probability — highest class probability = classification confidence.
var maxProb = probStack.reduce(ee.Reducer.max()).rename('max_prob');
 
// Shannon entropy — measures how spread the probability mass is across classes.
// Higher entropy = more uncertain. Normalized to [0, 1] by dividing by
// log(4) (maximum entropy for 4 equiprobable classes).
var epsilon = 1e-10;  // avoid log(0)
var entropy = probStack.add(epsilon)
    .multiply(probStack.add(epsilon).log())
    .multiply(-1)
    .reduce(ee.Reducer.sum())
    .rename('entropy');
 
var maxEntropy = Math.log(4);
var entropyNorm = entropy.divide(maxEntropy).rename('entropy_norm');
 
 
// =============================================================================
// 14. ASSEMBLE AND EXPORT
// =============================================================================
 
var finalOutput = hardClass
    .addBands(prob_evergreen)
    .addBands(prob_deciduous)
    .addBands(prob_herbaceous)
    .addBands(prob_shrub)
    .addBands(setSize)
    .addBands(inSetEvergreen)
    .addBands(inSetDeciduous)
    .addBands(inSetHerbaceous)
    .addBands(inSetShrub)
    .addBands(maxProb)
    .addBands(entropyNorm)
    .toFloat();
 
finalOutput = finalOutput.updateMask(unionMask);
 
// --- Export to Asset ---
Export.image.toAsset({
  image: finalOutput,
  description: 'Veg_Classification_v2_ToAsset',
  assetId: ASSET_FOLDER + '/Veg_Classification_v2',
  region: geometry,
  scale: SCALE,
  crs: CRS,
  maxPixels: 1e13
});
 
// --- Export to Drive ---
Export.image.toDrive({
  image: finalOutput,
  description: 'Veg_Classification_v2_ToDrive',
  folder: DRIVE_FOLDER,
  fileNamePrefix: 'Veg_Classification_v2',
  region: geometry,
  scale: SCALE,
  crs: CRS,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});
 
// --- Export model configuration ---
var configExport = ee.FeatureCollection([
  ee.Feature(null, {
    model_type: 'RandomForest',
    numberOfTrees: SELECTED_RF.numberOfTrees,
    variablesPerSplit: SELECTED_RF.variablesPerSplit,
    minLeafPopulation: SELECTED_RF.minLeafPopulation,
    bagFraction: SELECTED_RF.bagFraction,
    maxNodes: SELECTED_RF.maxNodes,
    seed: SELECTED_RF.seed
  }),
  ee.Feature(null, {
    model_type: 'GradientBoostedTrees',
    numberOfTrees: SELECTED_GBT.numberOfTrees,
    shrinkage: SELECTED_GBT.shrinkage,
    maxNodes: SELECTED_GBT.maxNodes,
    samplingRate: SELECTED_GBT.samplingRate,
    loss: SELECTED_GBT.loss,
    seed: SELECTED_GBT.seed
  })
]);
 
Export.table.toDrive({
  collection: configExport,
  description: 'VegClass_Model_Configuration_v2',
  folder: DRIVE_FOLDER,
  fileNamePrefix: 'VegClass_Model_Configuration_v2',
  fileFormat: 'CSV'
});
 
// --- Export conformal info ---
Export.table.toDrive({
  collection: ee.FeatureCollection([ee.Feature(null, {
    coverage_target: CONFORMAL_COVERAGE,
    conformal_threshold: CONFORMAL_THRESHOLD,
    n_training_samples: samplesWithId.size(),
    n_models: 2,
    crs: CRS,
    scale: SCALE
  })]),
  description: 'VegClass_Conformal_Info_v2',
  folder: DRIVE_FOLDER,
  fileNamePrefix: 'VegClass_Conformal_Info_v2',
  fileFormat: 'CSV'
});
 
 
// =============================================================================
// 15. VISUALIZATION
// =============================================================================
 
Map.centerObject(geometry, 11);
Map.setOptions('HYBRID');
 
Map.addLayer(hardClass.clip(geometry),
  {min: 0, max: 3, palette: CLASS_PALETTE},
  'Hard Classification', true);
 
Map.addLayer(setSize.clip(geometry),
  {min: 0, max: 4, palette: ['#d73027', '#1a9850', '#fee08b', '#fc8d59', '#d73027']},
  'Prediction Set Size', false);
 
Map.addLayer(entropyNorm.clip(geometry),
  {min: 0, max: 1, palette: ['#ffffb2', '#fecc5c', '#fd8d3c', '#f03b20', '#bd0026']},
  'Normalized Entropy', false);
 
Map.addLayer(prob_evergreen.clip(geometry),
  {min: 0, max: 1, palette: ['#ffffff', '#1b7837']}, 'P(Evergreen)', false);
Map.addLayer(prob_deciduous.clip(geometry),
  {min: 0, max: 1, palette: ['#ffffff', '#fdb863']}, 'P(Deciduous)', false);
Map.addLayer(prob_herbaceous.clip(geometry),
  {min: 0, max: 1, palette: ['#ffffff', '#d9d9d9']}, 'P(Herbaceous)', false);
Map.addLayer(prob_shrub.clip(geometry),
  {min: 0, max: 1, palette: ['#ffffff', '#762a83']}, 'P(Shrub)', false);
 
 
// =============================================================================
// 16. LEGEND
// =============================================================================
 
var legend = ui.Panel({
  style: {position: 'bottom-left', padding: '8px 15px'}
});
legend.add(ui.Label('Vegetation Classes', {
  fontWeight: 'bold', fontSize: '14px', margin: '0 0 8px 0'
}));
 
var classNames = ['Evergreen', 'Deciduous', 'Herbaceous', 'Shrub'];
classNames.forEach(function(name, idx) {
  var colorBox = ui.Label({
    style: {
      backgroundColor: CLASS_PALETTE[idx],
      padding: '10px 16px', margin: '0 8px 4px 0', border: '1px solid #999'
    }
  });
  var label = ui.Label(name, {margin: '4px 0', fontSize: '12px'});
  legend.add(ui.Panel({
    widgets: [colorBox, label],
    layout: ui.Panel.Layout.Flow('horizontal')
  }));
});
 
Map.add(legend);
 
 
