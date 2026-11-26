/**
 * predictor_stack_module.js
 * -------------------------
 * Generates a single multi-band predictor ee.Image for a given YEAR and AOI.
 *
 * No explicit reproject() calls are made.
 * Each band retains its native CRS internally, and GEE handles alignment
 *
 * The one exception is the ecoregion band, which is rasterized from vector
 * features via reduceToImage. That operation needs a defined resolution
 * context so the rasterization captures boundary detail. We use
 * setDefaultProjection at 100 m (adequate for L3 ecoregions) rather than
 * a hard reproject, so the band can still be lazily re-requested at finer
 * scales without triggering double reprojection.
 *
 * Output bands (24 total):
 *   10 spectral    — Sentinel-2 June-July median, renamed with _06_07 suffix
 *    6 phenology   — NDVI + EVI for Mar-Apr, Jun-Jul, Sep-Oct
 *    1 ecoregion   — EPA Level 3 ecoregion as numeric code
 *    4 terrain     — Elevation, Slope, Aspect_sin, Aspect_cos from SRTM
 *    3 climate     — tmean_annual, ppt_annual, temp_range from PRISM normals
 *
 * Storage format:
 *   All bands are Int16 or Uint16, not Float32. This halves memory per band
 *   with no loss for the ML pipeline — tree models are scale-invariant.
 *   Scale factors are documented in exports.SCALE_FACTORS for any consumer
 *   that needs physical units.
 *
 *   | Band(s)              | Type   | Scale Factor | Physical units          |
 *   |----------------------|--------|--------------|-------------------------|
 *   | Spectral (B2-B12)    | Int16  | 10000        | reflectance (0-1)       |
 *   | NDVI, EVI            | Int16  | 10000        | index (-1 to 1)         |
 *   | Ecoregion_L3         | Uint16 | 1            | integer code            |
 *   | Elevation            | Int16  | 1            | meters                  |
 *   | Slope                | Int16  | 100          | degrees                 |
 *   | Aspect_sin, cos      | Int16  | 10000        | sin/cos (-1 to 1)       |
 *   | tmean_annual         | Int16  | 100          | degrees C               |
 *   | ppt_annual           | Uint16 | 1            | mm                      |
 *   | temp_range           | Int16  | 100          | degrees C               |
 *
 * Terrain note:
 *   Raw aspect is circular (0° = 360° = north), which breaks regression
 *   models because 1° and 359° look very different numerically but point
 *   nearly the same direction. We decompose into sin/cos components:
 *     sin → positive = east-facing,  cos → positive = north-facing
 *   This preserves angular relationships without discontinuity.
 *
 * Climate note:
 *   PRISM 30-year normals (1991-2020) are monthly values. We aggregate to
 *   three annual summaries that capture the climate envelope relevant to
 *   vegetation structure:
 *     tmean_annual  — mean of 12 monthly mean temperatures
 *     ppt_annual    — sum of 12 monthly precipitation values
 *     temp_range    — warmest month mean minus coldest month mean
 *
 * CRS note:
 *   Previous version used EPSG:6350 (NAD83 CONUS Albers) with explicit
 *   reproject() on every sub-stack. That caused angled pixels in the GEE
 *   map viewer (Web Mercator basemap) and added unnecessary reprojection
 *   overhead. 
 */
 
// =============================================================================
// GLOBAL CONSTANTS
// =============================================================================
 
var S2_COLLECTION = 'COPERNICUS/S2_SR_HARMONIZED';
var CS_PLUS_COLLECTION = 'GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED';
 
var SRTM_ASSET = 'USGS/SRTMGL1_003';
var PRISM_NORMALS_ASSET = 'OREGONSTATE/PRISM/Norm91m';
 
// Cloud Score+ threshold — pixels below this cs_cdf value are masked.
// 0.60 is moderately aggressive; captures most thin cloud and haze while
// retaining enough pixels for stable median composites in the Intermountain West.
var CLOUD_FILTER = 80;         // scene-level: skip scenes above 80% cloud cover
var CS_CLEAR_THRESHOLD = 0.60; // pixel-level: mask where cs_cdf < 0.60
 
var FIXED_BANDS = ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B11', 'B12'];
 
 
// =============================================================================
// 1. SENTINEL-2 UTILITIES
// =============================================================================
 
/**
 * Applies pixel-level cloud mask via Cloud Score+ cs_cdf band.
 * Sentinel-2 SR reflectance is already stored as Int16 scaled by 10000
 * natively, so no additional scaling is needed — just select and mask.
 */
var maskCloudsCSPlus = function(image) {
  var csCdf = image.select('cs_cdf');
  var clearMask = csCdf.gte(CS_CLEAR_THRESHOLD);
  return image.updateMask(clearMask)
              .select(['B.*'])
              .toInt16()
              .copyProperties(image, ['system:time_start']);
};
 
/**
 * Computes NDVI and EVI from Sentinel-2 bands.
 *
 * Input bands are Int16 scaled by 10000 (native S2 SR format).
 * Output indices are computed in float, then rescaled to Int16 * 10000.
 *
 * EVI denominator: the standard formula uses a constant of 1.0 in the
 * denominator (NIR + 6*Red - 7.5*Blue + 1). Since our inputs are scaled
 * by 10000, the constant becomes 10000 to maintain unit consistency.
 */
var addIndices = function(image) {
  var b2 = image.select('B2').toFloat();
  var b4 = image.select('B4').toFloat();
  var b8 = image.select('B8').toFloat();
 
  var ndvi = b8.subtract(b4).divide(b8.add(b4))
      .multiply(10000)
      .toInt16()
      .rename('NDVI');
 
  var evi = b8.subtract(b4)
      .divide(b8.add(b4.multiply(6)).subtract(b2.multiply(7.5)).add(10000))
      .multiply(2.5)
      .multiply(10000)
      .toInt16()
      .rename('EVI');
 
  return image.addBands([ndvi, evi]);
};
 
/**
 * Filters and masks a Sentinel-2 collection for a given year and month range.
 *
 * Applies two levels of cloud filtering:
 *   1. Scene-level: CLOUDY_PIXEL_PERCENTAGE < CLOUD_FILTER (80%)
 *   2. Pixel-level: Cloud Score+ cs_cdf >= CS_CLEAR_THRESHOLD (0.60)
 *
 * Cloud Score+ is linked (not joined) — if cs_cdf is unavailable for a
 * given scene, that scene is still included but without pixel-level masking.
*/

var filterAndMaskCSPlus = function(year, startMonth, endMonth, geometry) {
  year = ee.Number(year);
  startMonth = ee.Number(startMonth);
  endMonth = ee.Number(endMonth);
 
  // Build date strings for server-side evaluation.
  // Using string concatenation avoids ee.Date.fromYMD issues in module context.
  var startStr = year.format('%d').cat('-')
      .cat(startMonth.format('%02d')).cat('-01');
  var endStr = year.format('%d').cat('-')
      .cat(endMonth.format('%02d')).cat('-01');
 
  var start = ee.Date(startStr);
  var end = ee.Date(endStr).advance(1, 'month');
 
  var s2 = ee.ImageCollection(S2_COLLECTION)
      .filterBounds(geometry)
      .filterDate(start, end)
      .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', CLOUD_FILTER));
 
  var csPlus = ee.ImageCollection(CS_PLUS_COLLECTION)
      .filterBounds(geometry)
      .filterDate(start, end);
 
  return s2.linkCollection(csPlus, ['cs_cdf']).map(maskCloudsCSPlus);
};
 
 
// =============================================================================
// 2. TERRAIN (SRTM)
// =============================================================================
 
/**
 * Generates terrain predictor bands from SRTM 30m DEM.
 *
 *
 * Returns 4 bands:
 *   Elevation   — meters, Int16
 *   Slope       — degrees × 100, Int16 (divide by 100 for degrees)
 *   Aspect_sin  — sin(aspect) × 10000, Int16
 *   Aspect_cos  — cos(aspect) × 10000, Int16
 */
var getTerrainStack = function(geometry) {
  var srtm = ee.Image(SRTM_ASSET);
  var terrain = ee.Terrain.products(srtm);
 
  var elevation = terrain.select('elevation').toInt16().rename('Elevation');
 
  // Slope: 0-90°, multiply by 100 for 0.01° precision in Int16
  var slope = terrain.select('slope')
      .multiply(100)
      .toInt16()
      .rename('Slope');
 
  // Aspect: decompose circular variable into sin/cos.
  // Scale by 10000 for Int16 storage of -1 to 1 range.
  var aspect = terrain.select('aspect');
  var aspectRadians = aspect.multiply(Math.PI).divide(180);
  var aspectSin = aspectRadians.sin()
      .multiply(10000).toInt16().rename('Aspect_sin');
  var aspectCos = aspectRadians.cos()
      .multiply(10000).toInt16().rename('Aspect_cos');
 
  return elevation
      .addBands(slope)
      .addBands(aspectSin)
      .addBands(aspectCos);
};
 
 
// =============================================================================
// 3. CLIMATE (PRISM 30-year normals)
// =============================================================================
 
/**
 * Generates climate predictor bands from PRISM 1991-2020 normals (~4 km).
 *
 * No explicit reprojection. PRISM's native resolution is ~800m (stored
 * at ~4km in GEE's asset). GEE resamples lazily at request time.
 * Climate variables are spatially smooth, so bilinear resampling to 10m
 * produces sensible values — unlike categorical data, where nearest-neighbor
 * is required.
 *
 * Returns 3 bands:
 *   tmean_annual — mean annual temperature × 100, Int16 (÷100 for °C)
 *   ppt_annual   — total annual precipitation, Uint16 (mm)
 *   temp_range   — warmest − coldest month mean × 100, Int16 (÷100 for °C)
 */
var getClimateStack = function(geometry) {
  var prismNormals = ee.ImageCollection(PRISM_NORMALS_ASSET);
 
  var tmeanAnnual = prismNormals
      .select('tmean')
      .mean()
      .multiply(100)
      .toInt16()
      .rename('tmean_annual');
 
  // Total annual precipitation. Uint16 (0-65535 mm) is sufficient;
  // the wettest location in the Intermountain West is well under 2000 mm/yr.
  var pptAnnual = prismNormals
      .select('ppt')
      .sum()
      .toUint16()
      .rename('ppt_annual');
 
  var tmeanMax = prismNormals.select('tmean').max();
  var tmeanMin = prismNormals.select('tmean').min();
  var tempRange = tmeanMax.subtract(tmeanMin)
      .multiply(100)
      .toInt16()
      .rename('temp_range');
 
  return tmeanAnnual
      .addBands(pptAnnual)
      .addBands(tempRange);
};
 
 
// =============================================================================
// 4. MAIN EXPORT: getSentinelPredictorStack
// =============================================================================
 
/**
 * Builds the complete 24-band predictor stack for one year over the given
 * geometry. No explicit reprojection — the consumer controls CRS and scale
 * at sampling or export time.
 *
 * @param {ee.Number|number} YEAR - Calendar year.
 * @param {ee.Geometry} geometry - Area of interest (used to filter S2 scenes
 *   and clip the output; does not set CRS).
 * @returns {ee.Image} 24-band predictor image in mixed native CRS. The
 *   consumer should specify scale and CRS when sampling or exporting.
 */
exports.getSentinelPredictorStack = function(YEAR, geometry) {
 
  // -------------------------------------------------------------------------
  // 4.1 Fixed spectral composite (June-July median)
  // -------------------------------------------------------------------------
  // June-July captures peak growing season in the Intermountain West.
  // Median compositing over two months produces a cloud-free, spectrally
  // representative pixel for most of the study area.
  var fixedComposite = filterAndMaskCSPlus(YEAR, 6, 7, geometry);
 
  var FIXED_NAMES = ee.List(FIXED_BANDS.map(function(b) { return b + '_06_07'; }));
 
  // Fallback: fully masked image with correct band names. If no S2 scenes
  // pass filters (rare but possible in heavily clouded years), the pipeline
  // continues with this — those pixels will be masked in the final training
  // stack and won't enter the model.
  var EMPTY_COMPOSITE = ee.Image.constant(ee.List.repeat(0, FIXED_BANDS.length))
      .rename(FIXED_NAMES)
      .updateMask(ee.Image(0))
      .toInt16();
 
  var composite = ee.Image(
    ee.Algorithms.If(
      fixedComposite.size().gt(0),
      fixedComposite.median().select(FIXED_BANDS, FIXED_NAMES).toInt16(),
      EMPTY_COMPOSITE
    )
  );
 
  // -------------------------------------------------------------------------
  // 4.2 Phenology stack (3 seasonal periods × 2 indices = 6 bands)
  // -------------------------------------------------------------------------
  // Three phenological windows capture the green-up, peak, and senescence
  // signal. This is especially important for separating deciduous (strong
  // seasonal NDVI swing) from evergreen (flat) and shrub from herbaceous
  // (different green-up timing).
  var phenologyPeriods = [
    {label: 'MarApr', start: 3, end: 4},
    {label: 'JunJul', start: 6, end: 7},
    {label: 'SepOct', start: 9, end: 10}
  ];
 
  var periodImages = phenologyPeriods.map(function(p) {
    var col = filterAndMaskCSPlus(YEAR, p.start, p.end, geometry);
 
    var targetNames = ee.List([
      ee.String('NDVI_').cat(p.label),
      ee.String('EVI_').cat(p.label)
    ]);
 
    // Same masked-zero fallback as the spectral composite.
    var zeroImg = ee.Image.constant([0, 0])
        .rename(['NDVI', 'EVI'])
        .updateMask(ee.Image(0))
        .toInt16();
 
    var indicesImage = ee.Image(
      ee.Algorithms.If(
        col.size().gt(0),
        addIndices(col.median()).select(['NDVI', 'EVI']).rename(targetNames),
        zeroImg.rename(targetNames)
      )
    );
 
    return indicesImage;
  });
 
  var phenologyStack = ee.Image(periodImages[0])
      .addBands(ee.Image(periodImages[1]))
      .addBands(ee.Image(periodImages[2]));
 
  // -------------------------------------------------------------------------
  // 4.3 Ecoregion band (EPA Level 3)
  // -------------------------------------------------------------------------
  // reduceToImage creates a raster from vector features. Unlike the other
  // data sources, there's no pre-existing raster with a native projection,
  // so we use setDefaultProjection to tell GEE the intended resolution.
  // 100 m is plenty for L3 ecoregion boundaries, which are drawn at
  // ~1:250,000 scale. This avoids rasterizing at 1-degree default.
  var ecoregions = ee.FeatureCollection('EPA/Ecoregions/2013/L3')
      .filterBounds(geometry)
      .map(function(f) {
        return f.set('l3_numeric_code', ee.Number.parse(f.get('us_l3code')));
      })
      .reduceToImage({
        properties: ['l3_numeric_code'],
        reducer: ee.Reducer.first()
      })
      .toUint16()
      .rename('Ecoregion_L3')
      .setDefaultProjection({crs: 'EPSG:4326', scale: 100});
 
  // -------------------------------------------------------------------------
  // 4.4 Terrain stack (SRTM)
  // -------------------------------------------------------------------------
  var terrainStack = getTerrainStack(geometry);
 
  // -------------------------------------------------------------------------
  // 4.5 Climate stack (PRISM 30-year normals)
  // -------------------------------------------------------------------------
  var climateStack = getClimateStack(geometry);
 
  // -------------------------------------------------------------------------
  // 4.6 Assemble and clip
  // -------------------------------------------------------------------------
  // clip() constrains the spatial extent but does not set a CRS.
  var fullStack = composite
      .addBands(phenologyStack)
      .addBands(ecoregions)
      .addBands(terrainStack)
      .addBands(climateStack)
      .clip(geometry);
 
  return fullStack;
};
 
 
// =============================================================================
// 5. UTILITY EXPORTS
// =============================================================================
 
exports.getTerrainStack = getTerrainStack;
exports.getClimateStack = getClimateStack;
 
/**
 * Scale factors: divide stored Int16 value by this to get physical units.
 * Only needed if a consumer wants human-readable values (e.g., for display
 * or equation application). The ML pipeline operates on raw scaled values.
 */
exports.SCALE_FACTORS = {
  'B2_06_07': 10000, 'B3_06_07': 10000, 'B4_06_07': 10000,
  'B5_06_07': 10000, 'B6_06_07': 10000, 'B7_06_07': 10000,
  'B8_06_07': 10000, 'B8A_06_07': 10000, 'B11_06_07': 10000, 'B12_06_07': 10000,
  'NDVI_MarApr': 10000, 'EVI_MarApr': 10000,
  'NDVI_JunJul': 10000, 'EVI_JunJul': 10000,
  'NDVI_SepOct': 10000, 'EVI_SepOct': 10000,
  'Ecoregion_L3': 1,
  'Elevation': 1,
  'Slope': 100,
  'Aspect_sin': 10000,
  'Aspect_cos': 10000,
  'tmean_annual': 100,
  'ppt_annual': 1,
  'temp_range': 100
};
 
/**
 * Returns the ordered list of 24 predictor band names.
 * Matches the band order in getSentinelPredictorStack output.
 */
exports.getPredictorBandNames = function() {
  return [
    'B2_06_07', 'B3_06_07', 'B4_06_07', 'B5_06_07', 'B6_06_07',
    'B7_06_07', 'B8_06_07', 'B8A_06_07', 'B11_06_07', 'B12_06_07',
    'NDVI_MarApr', 'EVI_MarApr', 'NDVI_JunJul', 'EVI_JunJul',
    'NDVI_SepOct', 'EVI_SepOct',
    'Ecoregion_L3',
    'Elevation', 'Slope', 'Aspect_sin', 'Aspect_cos',
    'tmean_annual', 'ppt_annual', 'temp_range'
  ];
};
 
/**
 * Returns data type for each band (Int16 or Uint16).
 */
exports.getDataTypes = function() {
  return {
    'B2_06_07': 'Int16', 'B3_06_07': 'Int16', 'B4_06_07': 'Int16',
    'B5_06_07': 'Int16', 'B6_06_07': 'Int16', 'B7_06_07': 'Int16',
    'B8_06_07': 'Int16', 'B8A_06_07': 'Int16', 'B11_06_07': 'Int16',
    'B12_06_07': 'Int16',
    'NDVI_MarApr': 'Int16', 'EVI_MarApr': 'Int16',
    'NDVI_JunJul': 'Int16', 'EVI_JunJul': 'Int16',
    'NDVI_SepOct': 'Int16', 'EVI_SepOct': 'Int16',
    'Ecoregion_L3': 'Uint16',
    'Elevation': 'Int16', 'Slope': 'Int16',
    'Aspect_sin': 'Int16', 'Aspect_cos': 'Int16',
    'tmean_annual': 'Int16', 'ppt_annual': 'Uint16', 'temp_range': 'Int16'
  };
};
