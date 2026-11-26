/**
 * timeseries_and_rf.js
 * ---------------------
 * Multi-year canopy height prediction with conformally calibrated
 * prediction intervals.
 *
 * ARCHITECTURE (v2):
 *   One tuned Random Forest + one tuned Gradient Boosted Tree, selected
 *   via hyperparameter search over a candidate grid. This replaces the
 *   previous 20-model diverse ensemble (10 RF + 10 GBT with varied
 *   hyperparameters), which had two problems:
 *     1. Ensemble std from only two algorithm families was a poor
 *        uncertainty proxy (conformal threshold of 18.1 — unusable).
 *     2. No hyperparameters were optimized; they were hand-picked.
 *
 *   The new approach acknowledges that with two models, ensemble spread
 *   is just |RF - GBT| / 2.
 *   Conformal prediction provides the
 *   coverage guarantee
 *    The two models provide complementary inductive biases:
 *    RF for variance reduction, GBT for bias reduction.
 *
 * CONFORMAL PREDICTION:
 *   Two nonconformity score formulations are computed and compared:
 *     Normalized: score_i = |y_i - ŷ_i| / σ_i  → adaptive intervals
 *     Absolute:   score_i = |y_i - ŷ_i|         → fixed-width intervals
 *
 *   For a 2-model ensemble, absolute scores are likely better because
 *   the denominator in normalized scores (σ = |RF - GBT| / 2) is a
 *   single pairwise difference — too noisy to be informative. The script
 *   computes both so you can compare coverage on the held-out test set
 *   and make an empirical choice.
 *
 *   TERMINOLOGY: conformal prediction produces "prediction intervals"
 *   (PI), not confidence intervals or credible intervals. PIs bound a
 *   future observation; CIs bound a parameter estimate; credible
 *   intervals are Bayesian posterior summaries.
 *
 * DATA SPLITS:
 *   Three-way split with strict separation:
 *     Training (60%)     → model fitting (and internal tuning split)
 *     Calibration (20%)  → conformal threshold estimation
 *     Test (20%)         → coverage verification + accuracy metrics
 *
 *   Previous version used a 70/30 split where the test set did double
 *   duty for both calibration and validation, breaking the conformal
 *   coverage guarantee. That is fixed here.
 *
 * HYPERPARAMETER TUNING:
 *   The training set is further split internally (80/20) for tuning.
 *   Candidate RF and GBT configurations are trained on the 80% and
 *   evaluated on the 20%. The final models are then retrained on the
 *   full 60% training set.
 *
 *   This keeps the calibration and test sets completely independent
 *   of all model selection decisions.
 *
 * OUTPUT FORMAT:
 *   All height-related bands are stored as Int16 × 100 (0.01 m precision).
 *   This halves storage vs Float32 with no practical loss — CHM values
 *   range 0-50 m and 0.01 m precision is well below measurement noise.
 *
 *   Output bands per year:
 *     Predicted_CHM     — ensemble mean (RF + GBT) / 2
 *     Ensemble_StdDev   — std of RF and GBT predictions
 *     PI_Lower          — prediction interval lower bound (conformal)
 *     PI_Upper          — prediction interval upper bound (conformal)
 *     PI_Width          — interval width (upper - lower)
 *     RF_Prediction     — RF model prediction (diagnostic)
 *     GBT_Prediction    — GBT model prediction (diagnostic)
 *
 *   To recover meters: value / 100. The asset_loader module (built in
 *   Phase 4) handles this transparently for downstream consumers.
 *
 * CRS:
 *   EPSG:4326 (WGS84) at 10 m nominal scale. Replaces EPSG:6350
 *   (NAD83 CONUS Albers) which caused angled pixels in the map viewer.
 */
 
// =============================================================================
// 1. MODULE IMPORTS
// =============================================================================
 
var CHMStacker = require('users/eriknatijohnson/MRRMaidCarbon:finalstack/chm_stacker');
var PredictorPrep = require('users/eriknatijohnson/MRRMaidCarbon:finalstack/predictor_stack_module');
 
// =============================================================================
// 2. CONFIGURATION
// =============================================================================
 
// --- Output grid ---
var TARGET_SCALE = 10;
var TARGET_CRS = 'EPSG:4326';
 
// --- Int16 scale factor for height-related output bands ---
// Multiply meters by this before casting to Int16; divide by this to recover meters.
var CHM_SCALE = 100;
 
// --- Sampling ---
var SAMPLE_SIZE = 40000;
var OVERSAMPLE_FACTOR = 10;
 
// --- Three-way data split ---
// Strict separation: training for fitting, calibration for conformal threshold,
// test for coverage verification and accuracy metrics. No set does double duty.
var TRAIN_PROPORTION = 0.60;
var CALIBRATION_PROPORTION = 0.20;
// Test gets the remainder (0.20)
var SPLIT_SEED = 42;
 
// --- Internal tuning split (within training set) ---
// Used only during hyperparameter search. The final models are retrained
// on the full training set after selection.
var TUNE_TRAIN_FRACTION = 0.80;  // 80% of training for candidate fitting
var TUNE_SEED = 99;
 
// --- Conformal prediction ---
var CONFORMAL_ALPHA = 0.05;  // 95% coverage target
 
// --- Conformal score type ---
// 'absolute' → fixed-width PIs: prediction ± threshold (in meters)
// 'normalized' → adaptive PIs: prediction ± threshold × ensemble_std
var CONFORMAL_SCORE_TYPE = 'absolute';
 
// --- Asset paths ---
var AOI_ASSET_ID = 'projects/ee-eriknatijohnson/assets/studyarea';
var AOI = ee.FeatureCollection(AOI_ASSET_ID);
var geometry = AOI.geometry().bounds();
 
var YEARS_TO_PROCESS = [2017, 2018, 2019, 2020, 2021, 2022];
var CHM_COLLECTION_ID = 'projects/ee-eriknatijohnson/assets/LIDAR_YEARS';
var PREDICTION_YEARS = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025];
 
// Output location — update when mrrmaid-carbon project is registered
var ASSET_OUTPUT_FOLDER = 'projects/mrrmaid-carbon/assets/Pred_CHM_v2';
var DRIVE_OUTPUT_FOLDER = 'GEE_CHM_v2_Outputs';
 
var YEARLY_MASK_COLLECTION_ID = 'projects/ee-eriknatijohnson/assets/YEARLYMASKS';
 
// --- Minimum ensemble std for normalized scores ---
// Prevents division by zero when RF and GBT perfectly agree.
// 0.1 m (stored as 10 in Int16 units) is well below measurement noise.
var MIN_STD_EPSILON = 0.1;
 
/* 
// =============================================================================
// 2B. HYPERPARAMETER GRIDS (for tuning phase)
// =============================================================================
// Each array defines candidate configurations to evaluate. Candidates are
// trained on the internal tuning split (train_fit) and scored against
// train_val. Metrics print for all candidates; pick the best and enter
// the config into SELECTED_RF / SELECTED_GBT below.
//
// Grids are designed to sweep the most impactful axes per algorithm:
//   RF:  variablesPerSplit (feature randomness), minLeafPopulation (regularization)
//   GBT: shrinkage × numberOfTrees interaction, maxNodes (tree complexity)
 
var RF_GRID = [
  // --- Sweep feature subset size (most impactful RF hyperparameter) ---
  // Lower = more diversity per tree; higher = more correlated trees
  {id: 'RF_01', nTrees: 300, varsPerSplit: 3,  minLeaf: 1,  bag: 0.5, maxNodes: null, seed: 1},
  {id: 'RF_02', nTrees: 300, varsPerSplit: 5,  minLeaf: 1,  bag: 0.5, maxNodes: null, seed: 2},
  {id: 'RF_03', nTrees: 300, varsPerSplit: 8,  minLeaf: 1,  bag: 0.5, maxNodes: null, seed: 3},
  {id: 'RF_04', nTrees: 300, varsPerSplit: 12, minLeaf: 1,  bag: 0.5, maxNodes: null, seed: 4},
  // --- Sweep leaf regularization ---
  // Higher minLeaf = smoother predictions, less overfit, but potential underfit
  {id: 'RF_05', nTrees: 300, varsPerSplit: 6,  minLeaf: 5,  bag: 0.5, maxNodes: null, seed: 5},
  {id: 'RF_06', nTrees: 300, varsPerSplit: 6,  minLeaf: 10, bag: 0.5, maxNodes: null, seed: 6},
  {id: 'RF_07', nTrees: 300, varsPerSplit: 6,  minLeaf: 20, bag: 0.5, maxNodes: null, seed: 7},
  // --- Sweep bag fraction ---
  // Lower = more diversity between trees; higher = more stable but homogeneous
  {id: 'RF_08', nTrees: 300, varsPerSplit: 6,  minLeaf: 1,  bag: 0.3, maxNodes: null, seed: 8},
  {id: 'RF_09', nTrees: 300, varsPerSplit: 6,  minLeaf: 1,  bag: 0.7, maxNodes: null, seed: 9},
  // --- Tree depth constraint ---
  {id: 'RF_10', nTrees: 300, varsPerSplit: 6,  minLeaf: 1,  bag: 0.5, maxNodes: 50,   seed: 10},
  {id: 'RF_11', nTrees: 300, varsPerSplit: 6,  minLeaf: 1,  bag: 0.5, maxNodes: 100,  seed: 11},
  // --- More trees (diminishing returns test) ---
  {id: 'RF_12', nTrees: 500, varsPerSplit: 6,  minLeaf: 1,  bag: 0.5, maxNodes: null, seed: 12}
];
 
var GBT_GRID = [
  // --- Sweep learning rate (most impactful GBT hyperparameter) ---
  // Lower shrinkage = needs more trees but generalizes better
  {id: 'GBT_01', nTrees: 200, shrinkage: 0.01, maxNodes: 10, sampRate: 0.8, loss: 'LeastSquares', seed: 101},
  {id: 'GBT_02', nTrees: 200, shrinkage: 0.05, maxNodes: 10, sampRate: 0.8, loss: 'LeastSquares', seed: 102},
  {id: 'GBT_03', nTrees: 200, shrinkage: 0.10, maxNodes: 10, sampRate: 0.8, loss: 'LeastSquares', seed: 103},
  {id: 'GBT_04', nTrees: 200, shrinkage: 0.20, maxNodes: 10, sampRate: 0.8, loss: 'LeastSquares', seed: 104},
  // --- Low shrinkage + more trees (the classic GBT sweet spot) ---
  {id: 'GBT_05', nTrees: 500, shrinkage: 0.05, maxNodes: 10, sampRate: 0.8, loss: 'LeastSquares', seed: 105},
  {id: 'GBT_06', nTrees: 500, shrinkage: 0.01, maxNodes: 10, sampRate: 0.8, loss: 'LeastSquares', seed: 106},
  // --- Sweep tree complexity ---
  // maxNodes controls individual tree depth; lower = weaker learners
  {id: 'GBT_07', nTrees: 300, shrinkage: 0.05, maxNodes: 5,  sampRate: 0.8, loss: 'LeastSquares', seed: 107},
  {id: 'GBT_08', nTrees: 300, shrinkage: 0.05, maxNodes: 20, sampRate: 0.8, loss: 'LeastSquares', seed: 108},
  {id: 'GBT_09', nTrees: 300, shrinkage: 0.05, maxNodes: 50, sampRate: 0.8, loss: 'LeastSquares', seed: 109},
  // --- Stochastic GBT (row subsampling) ---
  // samplingRate < 1 adds variance, acts like bagging within boosting
  {id: 'GBT_10', nTrees: 300, shrinkage: 0.05, maxNodes: 10, sampRate: 0.5, loss: 'LeastSquares', seed: 110},
  {id: 'GBT_11', nTrees: 300, shrinkage: 0.05, maxNodes: 10, sampRate: 1.0, loss: 'LeastSquares', seed: 111},
  // --- Robust loss function ---
  // Huber loss is less sensitive to outlier training pixels (e.g., CHM errors)
  {id: 'GBT_12', nTrees: 300, shrinkage: 0.05, maxNodes: 10, sampRate: 0.8, loss: 'Huber', seed: 112}
];

*/


// =============================================================================
// 2C. SELECTED HYPERPARAMETERS (for production — updated after tuning)
// =============================================================================
// secondary tuning done on GBT to test more iterations using Huber loss function
// GBT best parameters: 500 trees, shrinkage=0.05, maxNodes=10, sampRate=0.8, loss=Huber
// GBT best result: RMSE=2.56, MAE=1.39, R2=0.596
// RF already saturated, best parameters are 6vars, 1 minleaf, 0.5 bag, maxnodes null, 500 trees
// RF best result: RMSE: 2.60, MAE: 1.54, R2: 0.584
 
var SELECTED_RF = {
  numberOfTrees: 500,
  variablesPerSplit: 6,
  minLeafPopulation: 1,
  bagFraction: 0.5,
  maxNodes: null,  // null = unlimited
  seed: 42
};
 
var SELECTED_GBT = {
  numberOfTrees: 500,
  shrinkage: 0.05,
  maxNodes: 10,
  samplingRate: 0.8,
  loss: 'Huber',
  seed: 43
};
 
 
// =============================================================================
// 3. CLEAN YEARLY MASKS
// =============================================================================
// Yearly masks from MRRMaid/Landforms define where riparian/valley bottom
// pixels are valid for each year. Applied during prediction, not training
// (training extent is controlled by CHM footprint).
//
// Raw mask values: 1 = valid, 65535 = also valid (encoding quirk).
// Everything else is masked.
 
var rawYearlyMasks = ee.ImageCollection(YEARLY_MASK_COLLECTION_ID);
 
var YEARLY_MASK_COLLECTION = rawYearlyMasks.map(function(img) {
  img = ee.Image(img).select(0);
  var clean = img.eq(1).or(img.eq(65535)).rename('mask').toUint8();
  return clean.copyProperties(img, ['system:time_start', 'system:id', 'system:index']);
});
 
 
// =============================================================================
// 4. BUILD MULTI-YEAR TRAINING STACK
// =============================================================================
// Mosaic all LiDAR years into a single training image. Later years overwrite
// earlier years where they overlap (sort by time, most recent on top).
// This gives each pixel the most recent CHM value available.
 
var yearlyTrainingImages = ee.List(YEARS_TO_PROCESS).map(function(year) {
  year = ee.Number(year).toInt();
  var img = CHMStacker.getTrainingStack(year, geometry, CHM_COLLECTION_ID);
  return ee.Image(img).set('system:time_start', ee.Date.fromYMD(year, 1, 1).millis());
});
 
var yearlyCollection = ee.ImageCollection.fromImages(yearlyTrainingImages);
var sortedCollection = yearlyCollection.sort('system:time_start', false);
var finalTrainingStack = sortedCollection.mosaic();
 
var target = 'CanopyHeight';
var predictors = ee.List(PredictorPrep.getPredictorBandNames());
 
 
// =============================================================================
// 5. SAMPLE AND THREE-WAY SPLIT
// =============================================================================
 
// --- 5A. Sample pixels with valid CHM data ---
// Two-step sampling: first draw candidate locations from the CHM band
// (which is masked to LiDAR footprint), then extract predictor values
// at those locations. This ensures we only sample where training labels exist.
 
var candidateCount = ee.Number(SAMPLE_SIZE).multiply(OVERSAMPLE_FACTOR);
 
var chmCandidates = finalTrainingStack.select(target).sample({
  region: geometry,
  scale: TARGET_SCALE,
  numPixels: candidateCount,
  seed: 42,
  tileScale: 8,
  geometries: true
});
 
var predictorImage = finalTrainingStack.select(predictors);
var sampledWithPredictors = predictorImage.sampleRegions({
  collection: chmCandidates,
  scale: TARGET_SCALE,
  tileScale: 8,
  geometries: true
});
 
// Drop any samples with null values (edge pixels, projection artifacts)
var allBandsList = ee.List([target]).cat(predictors);
var validSamples = sampledWithPredictors
    .filter(ee.Filter.notNull(allBandsList))
    .randomColumn('rand', 12345)
    .sort('rand')
    .limit(SAMPLE_SIZE);
 
// --- 5B. Three-way split ---
// Assign each sample to training, calibration, or test using a random
// column. The boundaries are cumulative: [0, 0.60) = train,
// [0.60, 0.80) = calibration, [0.80, 1.0) = test.
var dataWithSplit = validSamples.randomColumn('split_rand', SPLIT_SEED);
 
var calBoundary = TRAIN_PROPORTION + CALIBRATION_PROPORTION;
 
var training    = dataWithSplit.filter(ee.Filter.lt('split_rand', TRAIN_PROPORTION));
var calibration = dataWithSplit.filter(ee.Filter.and(
    ee.Filter.gte('split_rand', TRAIN_PROPORTION),
    ee.Filter.lt('split_rand', calBoundary)
));
var testing = dataWithSplit.filter(ee.Filter.gte('split_rand', calBoundary));
 
 
/*
// =============================================================================
// 6. HYPERPARAMETER TUNING (optional — comment out after selection)
// =============================================================================
// Splits the training set internally (80/20) to evaluate candidate models.
// Prints RMSE and R² for each candidate. Review the output, pick the best
// RF and GBT configs, and enter them in SELECTED_RF / SELECTED_GBT above.
//
// This section is computationally expensive (trains ~24 models). After you've
// selected hyperparameters, comment it out to skip on subsequent runs.
 
// --- 6A. Internal tuning split ---
var tuneData = training.randomColumn('tune_rand', TUNE_SEED);
var trainFit = tuneData.filter(ee.Filter.lt('tune_rand', TUNE_TRAIN_FRACTION));
var trainVal = tuneData.filter(ee.Filter.gte('tune_rand', TUNE_TRAIN_FRACTION));
 

*/ 
// --- 6B. Helper: compute RMSE from a FeatureCollection with actual/predicted ---
// still need this function for later
function computeTuningMetrics(fc, actualProp, predProp) {
  fc = ee.FeatureCollection(fc).filter(ee.Filter.notNull([actualProp, predProp]));
  var withResiduals = fc.map(function(f) {
    var a = ee.Number(f.get(actualProp));
    var p = ee.Number(f.get(predProp));
    var res = a.subtract(p);
    return f.set('sq_res', res.multiply(res), 'abs_res', res.abs());
  });
  var mse = ee.Number(withResiduals.reduceColumns(ee.Reducer.mean(), ['sq_res']).get('mean'));
  var mae = ee.Number(withResiduals.reduceColumns(ee.Reducer.mean(), ['abs_res']).get('mean'));
  var variance = ee.Number(withResiduals.reduceColumns(ee.Reducer.variance(), [actualProp]).get('variance'));
  var rmse = mse.max(0).sqrt();
  var r2 = ee.Number(1).subtract(mse.divide(variance.max(1e-10)));
  return ee.Dictionary({RMSE: rmse, MAE: mae, R2: r2});
}

/* 
// --- 6C. Evaluate RF candidates ---
// Each candidate is exported as its own task. This isolates computation
// graphs so each task only has to train and evaluate a single model.
// A combined FeatureCollection of all 12 candidates causes GEE to build
// one monolithic graph that exceeds server-side limits.

 
for (var ri = 0; ri < RF_GRID.length; ri++) {
  var rCfg = RF_GRID[ri];
 
  var rfParams = {
    numberOfTrees: rCfg.nTrees,
    variablesPerSplit: rCfg.varsPerSplit,
    minLeafPopulation: rCfg.minLeaf,
    bagFraction: rCfg.bag,
    seed: rCfg.seed
  };
  if (rCfg.maxNodes !== null) {
    rfParams.maxNodes = rCfg.maxNodes;
  }
 
  var rfCandidate = ee.Classifier.smileRandomForest(rfParams)
      .setOutputMode('REGRESSION')
      .train({features: trainFit, classProperty: target, inputProperties: predictors});
 
  var rfPredImage = predictorImage.classify(rfCandidate).rename('pred');
  var rfValPred = rfPredImage.sampleRegions({
    collection: trainVal,
    properties: [target],
    scale: TARGET_SCALE,
    tileScale: 8
  });
 
  var rfMetrics = computeTuningMetrics(rfValPred, target, 'pred');
 
  // One Feature, one FeatureCollection, one export task — fully isolated graph.
  var rfResultFC = ee.FeatureCollection([ee.Feature(null, {
    id: rCfg.id,
    nTrees: rCfg.nTrees,
    varsPerSplit: rCfg.varsPerSplit,
    minLeaf: rCfg.minLeaf,
    bag: rCfg.bag,
    maxNodes: rCfg.maxNodes,
    RMSE: rfMetrics.get('RMSE'),
    MAE: rfMetrics.get('MAE'),
    R2: rfMetrics.get('R2')
  })]);
 
  Export.table.toDrive({
    collection: rfResultFC,
    description: 'RF_Tune_' + rCfg.id,
    folder: DRIVE_OUTPUT_FOLDER,
    fileNamePrefix: 'RF_Tune_' + rCfg.id,
    fileFormat: 'CSV'
  });
}
 
// --- 6D. Evaluate GBT candidates ---
// Same per-candidate export pattern as RF above. Each GBT candidate is
// trained and evaluated in its own isolated export task.

 
for (var gi = 0; gi < GBT_GRID.length; gi++) {
  var gCfg = GBT_GRID[gi];
 
  var gbtCandidate = ee.Classifier.smileGradientTreeBoost({
    numberOfTrees: gCfg.nTrees,
    shrinkage: gCfg.shrinkage,
    maxNodes: gCfg.maxNodes,
    samplingRate: gCfg.sampRate,
    loss: gCfg.loss,
    seed: gCfg.seed
  })
  .setOutputMode('REGRESSION')
  .train({features: trainFit, classProperty: target, inputProperties: predictors});
 
  var gbtPredImage = predictorImage.classify(gbtCandidate).rename('pred');
  var gbtValPred = gbtPredImage.sampleRegions({
    collection: trainVal,
    properties: [target],
    scale: TARGET_SCALE,
    tileScale: 8
  });
 
  var gbtMetrics = computeTuningMetrics(gbtValPred, target, 'pred');
 
  var gbtResultFC = ee.FeatureCollection([ee.Feature(null, {
    id: gCfg.id,
    nTrees: gCfg.nTrees,
    shrinkage: gCfg.shrinkage,
    maxNodes: gCfg.maxNodes,
    sampRate: gCfg.sampRate,
    loss: gCfg.loss,
    RMSE: gbtMetrics.get('RMSE'),
    MAE: gbtMetrics.get('MAE'),
    R2: gbtMetrics.get('R2')
  })]);
 
  Export.table.toDrive({
    collection: gbtResultFC,
    description: 'GBT_Tune_' + gCfg.id,
    folder: DRIVE_OUTPUT_FOLDER,
    fileNamePrefix: 'GBT_Tune_' + gCfg.id,
    fileFormat: 'CSV'
  });
}

*/ 
 
// =============================================================================
// 7. TRAIN FINAL MODELS
// =============================================================================
// Uses SELECTED_RF and SELECTED_GBT hyperparameters on the FULL training
// set (not the internal tuning split). This maximizes the data available
// for the production models.
 
 
// --- RF ---
var rfFinalParams = {
  numberOfTrees: SELECTED_RF.numberOfTrees,
  variablesPerSplit: SELECTED_RF.variablesPerSplit,
  minLeafPopulation: SELECTED_RF.minLeafPopulation,
  bagFraction: SELECTED_RF.bagFraction,
  seed: SELECTED_RF.seed
};
if (SELECTED_RF.maxNodes !== null) {
  rfFinalParams.maxNodes = SELECTED_RF.maxNodes;
}
 
var rfModel = ee.Classifier.smileRandomForest(rfFinalParams)
    .setOutputMode('REGRESSION')
    .train({features: training, classProperty: target, inputProperties: predictors});
 
// --- GBT ---
var gbtModel = ee.Classifier.smileGradientTreeBoost({
  numberOfTrees: SELECTED_GBT.numberOfTrees,
  shrinkage: SELECTED_GBT.shrinkage,
  maxNodes: SELECTED_GBT.maxNodes,
  samplingRate: SELECTED_GBT.samplingRate,
  loss: SELECTED_GBT.loss,
  seed: SELECTED_GBT.seed
})
.setOutputMode('REGRESSION')
.train({features: training, classProperty: target, inputProperties: predictors});
 
 
 
// =============================================================================
// 8. VARIABLE IMPORTANCE
// =============================================================================
// From the RF model only — GBT in GEE doesn't expose variable importance.
 
var explanation = rfModel.explain();
var importances = ee.Dictionary(explanation.get('importance'));
 
var importList = importances.keys().map(function(k) {
  return ee.Feature(null, {variable: k, importance: importances.get(k)});
});
 
var importanceFC = ee.FeatureCollection(importList).sort('importance', false);
 
Export.table.toDrive({
  collection: importanceFC,
  description: 'RF_Variable_Importance',
  folder: DRIVE_OUTPUT_FOLDER,
  fileNamePrefix: 'RF_Variable_Importance',
  fileFormat: 'CSV'
});
 
 
// =============================================================================
// 9. CONFORMAL CALIBRATION (on calibration set — never seen during training)
// =============================================================================
//
// The calibration set was held out before any model fitting or hyperparameter
// selection, preserving the exchangeability assumption required for the
// conformal coverage guarantee.
 
 
// --- 9A. Generate ensemble predictions at calibration locations ---
var rfPredCal = predictorImage.classify(rfModel).rename('rf_pred');
var gbtPredCal = predictorImage.classify(gbtModel).rename('gbt_pred');
 
var ensembleImage = rfPredCal.addBands(gbtPredCal);
var ensembleMeanCal = rfPredCal.add(gbtPredCal).divide(2).rename('ensemble_mean');
var ensembleStdCal = ee.ImageCollection.fromImages([
  rfPredCal.rename('pred'), gbtPredCal.rename('pred')
]).reduce(ee.Reducer.stdDev()).rename('ensemble_std');
 
var ensembleStatsCal = ensembleMeanCal.addBands(ensembleStdCal);
 
var calWithPred = ensembleStatsCal.sampleRegions({
  collection: calibration,
  properties: [target],
  scale: TARGET_SCALE,
  tileScale: 8
});
 
// --- 9B. Compute nonconformity scores (both formulations) ---
var calWithScores = calWithPred.map(function(f) {
  var trueVal = ee.Number(f.get(target));
  var predMean = ee.Number(f.get('ensemble_mean'));
  var predStd = ee.Number(f.get('ensemble_std')).max(MIN_STD_EPSILON);
  var absResidual = trueVal.subtract(predMean).abs();
  var normalizedScore = absResidual.divide(predStd);
  return f.set({
    'abs_score': absResidual,
    'norm_score': normalizedScore
  });
});
 
// --- 9C. Determine thresholds ---
// The conformal threshold is the (1-α)(1 + 1/n) quantile of the calibration
// scores. The finite-sample correction (+1/n) ensures coverage ≥ (1-α) rather
// than just asymptotic equality.
var nCal = calWithScores.size();
 
var adjustedQuantileLevel = nCal.add(1)
    .multiply(ee.Number(1).subtract(CONFORMAL_ALPHA))
    .ceil()
    .divide(nCal)
    .min(1.0);
 
// Absolute threshold (meters — directly interpretable as interval half-width)
var absScores = calWithScores.aggregate_array('abs_score').sort();
var absQuantileIndex = nCal.multiply(adjustedQuantileLevel).ceil()
    .subtract(1).max(0).min(nCal.subtract(1));
var absThreshold = ee.Number(absScores.get(absQuantileIndex));
 
// Normalized threshold (dimensionless — multiply by ensemble_std for interval)
var normScores = calWithScores.aggregate_array('norm_score').sort();
var normQuantileIndex = nCal.multiply(adjustedQuantileLevel).ceil()
    .subtract(1).max(0).min(nCal.subtract(1));
var normThreshold = ee.Number(normScores.get(normQuantileIndex));
 
 
// Diagnostic: ensemble std statistics on calibration set
var stdStats = calWithScores.reduceColumns({
  reducer: ee.Reducer.mean().combine(ee.Reducer.stdDev(), '', true)
      .combine(ee.Reducer.median(), '', true),
  selectors: ['ensemble_std']
});
 
// Select the threshold to use based on CONFORMAL_SCORE_TYPE.
// Client-side conditional — CONFORMAL_SCORE_TYPE is a constant string,
// so the branch is resolved at graph construction time.
var conformalThreshold;
if (CONFORMAL_SCORE_TYPE === 'absolute') {
  conformalThreshold = absThreshold;
} else {
  conformalThreshold = normThreshold;
}
 
 
// =============================================================================
// 10. COVERAGE VERIFICATION (on test set — never seen during anything)
// =============================================================================
// Check that the calibrated threshold achieves the target coverage on data
// that was held out from training, tuning, AND calibration.
 
 
var testWithPred = ensembleStatsCal.sampleRegions({
  collection: testing,
  properties: [target],
  scale: TARGET_SCALE,
  tileScale: 8
});
 
// Evaluate coverage for BOTH score types so you can compare
var testWithCoverage = testWithPred.map(function(f) {
  var trueVal = ee.Number(f.get(target));
  var predMean = ee.Number(f.get('ensemble_mean'));
  var predStd = ee.Number(f.get('ensemble_std')).max(MIN_STD_EPSILON);
 
  // Absolute PI: prediction ± absThreshold
  var absLower = predMean.subtract(absThreshold);
  var absUpper = predMean.add(absThreshold);
  var absCovered = trueVal.gte(absLower).and(trueVal.lte(absUpper));
 
  // Normalized PI: prediction ± normThreshold × ensemble_std
  var normHalfWidth = predStd.multiply(normThreshold);
  var normLower = predMean.subtract(normHalfWidth);
  var normUpper = predMean.add(normHalfWidth);
  var normCovered = trueVal.gte(normLower).and(trueVal.lte(normUpper));
 
  return f.set({
    'abs_covered': absCovered,
    'norm_covered': normCovered,
    'abs_pi_width': absThreshold.multiply(2),
    'norm_pi_width': normHalfWidth.multiply(2)
  });
});
 
var nTest = testWithCoverage.size();
 
// Absolute coverage
var absCoverageSum = testWithCoverage.reduceColumns(ee.Reducer.sum(), ['abs_covered']);
var absCoverage = ee.Number(absCoverageSum.get('sum')).divide(nTest);
var absAvgWidth = testWithCoverage.reduceColumns(ee.Reducer.mean(), ['abs_pi_width']);
 
// Normalized coverage
var normCoverageSum = testWithCoverage.reduceColumns(ee.Reducer.sum(), ['norm_covered']);
var normCoverage = ee.Number(normCoverageSum.get('sum')).divide(nTest);
var normAvgWidth = testWithCoverage.reduceColumns(ee.Reducer.mean(), ['norm_pi_width']);

 
 
// =============================================================================
// 11. VALIDATION METRICS (on test set)
// =============================================================================
// Standard accuracy metrics for the ensemble point prediction.
 
 
var testMetrics = computeTuningMetrics(testWithPred, target, 'ensemble_mean');

 
// Also compute per-model metrics for comparison
var rfPredTest = predictorImage.classify(rfModel).rename('rf_pred');
var testWithRF = rfPredTest.sampleRegions({
  collection: testing, properties: [target], scale: TARGET_SCALE, tileScale: 8
});
var rfTestMetrics = computeTuningMetrics(testWithRF, target, 'rf_pred');

 
var gbtPredTest = predictorImage.classify(gbtModel).rename('gbt_pred');
var testWithGBT = gbtPredTest.sampleRegions({
  collection: testing, properties: [target], scale: TARGET_SCALE, tileScale: 8
});
var gbtTestMetrics = computeTuningMetrics(testWithGBT, target, 'gbt_pred');
 
// Export validation metrics
var validationExport = ee.FeatureCollection([
  ee.Feature(null, {
    model: 'Ensemble',
    RMSE: testMetrics.get('RMSE'),
    MAE: testMetrics.get('MAE'),
    R2: testMetrics.get('R2'),
    abs_coverage: absCoverage,
    abs_threshold: absThreshold,
    abs_mean_pi_width: absAvgWidth.get('mean'),
    norm_coverage: normCoverage,
    norm_threshold: normThreshold,
    norm_mean_pi_width: normAvgWidth.get('mean'),
    conformal_score_type: CONFORMAL_SCORE_TYPE,
    n_training: training.size(),
    n_calibration: nCal,
    n_test: nTest
  })
]);
 
Export.table.toDrive({
  collection: validationExport,
  description: 'Validation_Metrics_v2',
  folder: DRIVE_OUTPUT_FOLDER,
  fileNamePrefix: 'Validation_Metrics_v2',
  fileFormat: 'CSV'
});
 
// Export scatter plot data for external visualization
var scatterData = testWithPred.map(function(f) {
  return ee.Feature(null, {
    actual: f.get(target),
    predicted: f.get('ensemble_mean'),
    ensemble_std: f.get('ensemble_std')
  });
});
 
Export.table.toDrive({
  collection: scatterData,
  description: 'Scatter_Data_v2',
  folder: DRIVE_OUTPUT_FOLDER,
  fileNamePrefix: 'Scatter_Data_v2',
  fileFormat: 'CSV'
});
 
 
// =============================================================================
// 12. BATCH PREDICTION AND EXPORT
// =============================================================================
// Generate annual CHM predictions with calibrated prediction intervals.
// Each year: load Sentinel predictors → apply yearly mask → classify with
// both models → compute ensemble stats → apply conformal threshold → export.
//
// Outputs are Int16 × 100 (divide by 100 for meters). See header for full
// band listing.
 
// --- 12A. Yearly mask helper ---
// Loads and applies the year-specific riparian/valley bottom mask.
// Only valid riparian pixels are retained in the output; everything else
// is masked and exports as nodata.
var applyYearlyMask = function(image, yearNumber) {
  var start = ee.Date.fromYMD(yearNumber, 1, 1);
  var end = start.advance(1, 'year');
 
  var masksByDate = YEARLY_MASK_COLLECTION
      .filterDate(start, end)
      .filterBounds(geometry);
 
  var maskImage = ee.Image(
    ee.Algorithms.If(
      masksByDate.size().gt(0),
      masksByDate.first().select(0),
      ee.Image.constant(1)
    )
  );
 
  return image.updateMask(maskImage.eq(1).selfMask());
};
 
// --- 12B. Prediction interval computation ---
// Encapsulated as a function so the score type switch happens in one place.
var computePredictionInterval = function(ensembleMean, ensembleStd, threshold, scoreType) {
  var halfWidth;
  if (scoreType === 'absolute') {
    // Fixed-width: ± threshold (in meters, same units as CHM)
    halfWidth = ee.Image.constant(threshold);
  } else {
    // Adaptive: ± threshold × ensemble_std
    var clampedStd = ensembleStd.max(MIN_STD_EPSILON);
    halfWidth = clampedStd.multiply(threshold);
  }
 
  var lower = ensembleMean.subtract(halfWidth);
  var upper = ensembleMean.add(halfWidth);
  var width = halfWidth.multiply(2);
 
  return {lower: lower, upper: upper, width: width};
};
 
// --- 12C. Export loop ---
var predictAndExportYear = function(year) {
  var yearNumber = ee.Number(year).toInt();
 
  // Build predictor stack for this year (CRS-agnostic, aligned at export time)
  var yearlyPredictorStack = PredictorPrep.getSentinelPredictorStack(yearNumber, geometry);
  yearlyPredictorStack = applyYearlyMask(yearlyPredictorStack, yearNumber);
  var inputImage = yearlyPredictorStack.select(predictors);
 
  // Classify with both models
  var rfPred = inputImage.classify(rfModel).rename('rf_pred');
  var gbtPred = inputImage.classify(gbtModel).rename('gbt_pred');
 
  // Ensemble statistics
  var ensembleMean = rfPred.add(gbtPred).divide(2).rename('Predicted_CHM');
  var ensembleStd = ee.ImageCollection.fromImages([
    rfPred.rename('pred'), gbtPred.rename('pred')
  ]).reduce(ee.Reducer.stdDev()).rename('Ensemble_StdDev');
 
  // Prediction intervals (conformal-calibrated)
  var pi = computePredictionInterval(
    ensembleMean, ensembleStd, conformalThreshold, CONFORMAL_SCORE_TYPE
  );
 
  // Assemble output bands, scale to Int16
  var outputImage = ensembleMean.multiply(CHM_SCALE).toInt16().rename('Predicted_CHM')
      .addBands(ensembleStd.multiply(CHM_SCALE).toInt16().rename('Ensemble_StdDev'))
      .addBands(pi.lower.multiply(CHM_SCALE).toInt16().rename('PI_Lower'))
      .addBands(pi.upper.multiply(CHM_SCALE).toInt16().rename('PI_Upper'))
      .addBands(pi.width.multiply(CHM_SCALE).toInt16().rename('PI_Width'))
      .addBands(rfPred.multiply(CHM_SCALE).toInt16().rename('RF_Prediction'))
      .addBands(gbtPred.multiply(CHM_SCALE).toInt16().rename('GBT_Prediction'));
 
  // Set metadata properties
  outputImage = outputImage.set({
    'year': year,
    'conformal_alpha': CONFORMAL_ALPHA,
    'conformal_score_type': CONFORMAL_SCORE_TYPE,
    'conformal_threshold': conformalThreshold,
    'chm_scale_factor': CHM_SCALE,
    'rf_nTrees': SELECTED_RF.numberOfTrees,
    'gbt_nTrees': SELECTED_GBT.numberOfTrees,
    'gbt_shrinkage': SELECTED_GBT.shrinkage
  });
 
  var yearStr = year.toString();
 
  Export.image.toAsset({
    image: outputImage,
    description: 'Predicted_CHM_v2_' + yearStr,
    assetId: ASSET_OUTPUT_FOLDER + '/Predicted_CHM_v2_' + yearStr,
    region: geometry,
    scale: TARGET_SCALE,
    crs: TARGET_CRS,
    maxPixels: 1e13
  });
};
 
PREDICTION_YEARS.forEach(predictAndExportYear);
 
 
// =============================================================================
// 13. EXPORT DIAGNOSTICS AND CONFIGURATION
// =============================================================================
 
// --- Model configuration ---
var configExport = ee.FeatureCollection([
  ee.Feature(null, {
    model_type: 'RandomForest',
    numberOfTrees: SELECTED_RF.numberOfTrees,
    variablesPerSplit: SELECTED_RF.variablesPerSplit,
    minLeafPopulation: SELECTED_RF.minLeafPopulation,
    bagFraction: SELECTED_RF.bagFraction,
    maxNodes: SELECTED_RF.maxNodes,
    seed: SELECTED_RF.seed,
    shrinkage: null,
    samplingRate: null,
    loss: null
  }),
  ee.Feature(null, {
    model_type: 'GradientBoostedTrees',
    numberOfTrees: SELECTED_GBT.numberOfTrees,
    variablesPerSplit: null,
    minLeafPopulation: null,
    bagFraction: null,
    maxNodes: SELECTED_GBT.maxNodes,
    seed: SELECTED_GBT.seed,
    shrinkage: SELECTED_GBT.shrinkage,
    samplingRate: SELECTED_GBT.samplingRate,
    loss: SELECTED_GBT.loss
  })
]);
 
Export.table.toDrive({
  collection: configExport,
  description: 'Model_Configuration_v2',
  folder: DRIVE_OUTPUT_FOLDER,
  fileNamePrefix: 'Model_Configuration_v2',
  fileFormat: 'CSV'
});
 
