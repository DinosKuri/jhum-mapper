//==============================================================
// MIZORAM JHUM — PART I: ANNUAL SERIES, FALLOW AGE AND CYCLE
//
// Builds the v6 jhum map for EVERY year, then derives, fully
// automatically:
//   - jhum_count        how many times each pixel was cleared
//   - last_clearing     year of the most recent clearing
//   - fallow_age        years since the last clearing (right-censored)
//   - first_clearing    year of the earliest detected clearing
//   - cycle_mean        mean years between successive clearings
//   - cycle_last        length of the most recent completed cycle
//
// Post-processing applied to every annual map (all automatic):
//   1. Minimum mapping unit: patches smaller than MIN_PIXELS are
//      dropped (speckle from leaf-off is small and scattered; jhum
//      plots are ~1 ha = ~11 Landsat pixels).
//   2. Optional recovery check: a true clearing greens up again the
//      following dry season; a permanently bare surface does not.
//
// Rule (calibrated on 400 points, validated on 310):
//   jhum if (NBR_post_wide <= 0.26 and dNBR_wide >= 0.20)
//        or (NBR_post_narrow <= 0.14 and dNBR_narrow >= 0.46)
//==============================================================

//--------------------------------------------------------------
// 0. SETTINGS
//--------------------------------------------------------------
var START_YEAR = 1989;        // 1988 has no usable "before" imagery — excluded
var END_YEAR   = 2025;
var VERSION    = 'v6';
// Area of interest: use YOUR OWN boundary shapefile if you have one
// (upload the shp as a GEE table asset); otherwise the GAUL state boundary.
var ROI_ASSET  = '';          // e.g. 'projects/YOUR_PROJECT/assets/Mizoram_boundary'
var ROI_NAME   = 'Mizoram';   // used only when ROI_ASSET is ''
var DRIVE_DIR  = 'Jhum_series_v6';
var ASSET_DIR  = 'projects/YOUR_PROJECT/assets/';   // for the annual stack

// Years used for the CYCLE metrics. Coverage is incomplete before 1998
// (1989, 1991-92, 1995-97 are 65-90% complete), so cycle lengths computed
// from those years would be biased long (clearings were simply unseen).
// Fallow age and last_clearing still use ALL years.
var METRIC_START = 1998;

// post-processing
var MIN_PIXELS = 1;        // 1 = no minimum mapping unit (the validated configuration).
                           // 5 (≈0.45 ha) was evaluated and removes ~5 % of area; it is
                           // kept as an option but the reported series uses 1.
var RECOVERY_CHECK = false; // true = require greening in the following year
var RECOVERY_MIN   = 0.10;  // minimum NBR increase next year if used

// outputs
var EXPORT_STACK_ASSET = true;   // one multi-band asset: jhum_<year> per band
var EXPORT_SUMMARY     = true;   // age / count / cycle rasters to Drive
var EXPORT_AREA_TABLE  = true;

// rule and compositing (locked)
var RULE = {wide: {post: 0.26, dnbr: 0.20}, narrow: {post: 0.14, dnbr: 0.46}, useNarrow: true};
var W_PRE_S = [-2, 11, 1], W_PRE_E = [-1, 6, 1];
var W_POST_S = [-1, 11, 1], W_POST_E = [0, 6, 1];
var N_PRE_S = [-1, 3, 1], N_PRE_E = [-1, 6, 1];
var N_POST_S = [0, 3, 1], N_POST_E = [0, 6, 1];
var PCT = 25;

var YEARS = [];
for (var yy = START_YEAR; yy <= END_YEAR; yy++) YEARS.push(yy);

//--------------------------------------------------------------
// 1. AREA AND HARMONIZED LANDSAT
//--------------------------------------------------------------
var roi = ROI_ASSET ? ee.FeatureCollection(ROI_ASSET)
        : ee.FeatureCollection('FAO/GAUL/2015/level1')
            .filter(ee.Filter.eq('ADM1_NAME', ROI_NAME));
var ROI_MASK = ee.Image.constant(1).clip(roi).mask();   // hard clip for every output
Map.centerObject(roi, 8);
var WATER_FREE = ee.Image('JRC/GSW1_4/GlobalSurfaceWater')
                   .select('occurrence').unmask(0).lt(50);

function prepL(img, nir, swir2, harmonize) {
  img = ee.Image(img);
  var qa = img.select('QA_PIXEL');
  var clear = qa.bitwiseAnd(1 << 1).eq(0).and(qa.bitwiseAnd(1 << 3).eq(0))
              .and(qa.bitwiseAnd(1 << 4).eq(0));
  var out = img.select([nir, swir2], ['NIR', 'SWIR2'])
               .multiply(0.0000275).add(-0.2).updateMask(clear);
  if (harmonize) out = out.multiply([0.8462, 0.9071]).add([0.0412, 0.0172]);
  return ee.Image(out.toFloat().copyProperties(img, ['system:time_start']));
}
var LS = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2').filterBounds(roi)
           .map(function (i) { return prepL(i, 'SR_B4', 'SR_B7', true); })
  .merge(ee.ImageCollection('LANDSAT/LE07/C02/T1_L2').filterBounds(roi)
           .map(function (i) { return prepL(i, 'SR_B4', 'SR_B7', true); }))
  .merge(ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filterBounds(roi)
           .map(function (i) { return prepL(i, 'SR_B5', 'SR_B7', false); }))
  .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2').filterBounds(roi)
           .map(function (i) { return prepL(i, 'SR_B5', 'SR_B7', false); }));

function win(y, w) { return ee.Date.fromYMD(y + w[0], w[1], w[2]); }
function nbrOf(i) { return ee.Image(i).normalizedDifference(['NIR', 'SWIR2']).rename('NBR'); }
function pctImg(col) {
  return ee.Image(ee.Algorithms.If(col.size().gt(0),
           col.map(nbrOf).reduce(ee.Reducer.percentile([PCT])),
           ee.Image.constant(0).updateMask(0))).rename('NBR').toFloat();
}
// Safe scene count: an EMPTY collection returns an image with NO bands,
// which breaks unmask(). Early years (1988-89) can have empty windows.
function safeCount(col) {
  return ee.Image(ee.Algorithms.If(col.size().gt(0),
           col.select('NIR').count(), ee.Image.constant(0))).rename('n').unmask(0);
}
function postWide(y)  { return pctImg(LS.filterDate(win(y, W_POST_S), win(y, W_POST_E))); }
function preWide(y)   { return pctImg(LS.filterDate(win(y, W_PRE_S),  win(y, W_PRE_E))); }
function postNarrow(y){ return pctImg(LS.filterDate(win(y, N_POST_S), win(y, N_POST_E))); }
function preNarrow(y) { return pctImg(LS.filterDate(win(y, N_PRE_S),  win(y, N_PRE_E))); }

//--------------------------------------------------------------
// 2. ANNUAL JHUM MAP + AUTOMATIC POST-PROCESSING
//--------------------------------------------------------------
function rawJhum(y) {
  var pW = postWide(y), dW = preWide(y).subtract(pW);
  var pN = postNarrow(y), dN = preNarrow(y).subtract(pN);
  var A = pW.lte(RULE.wide.post).and(dW.gte(RULE.wide.dnbr));
  var B = pN.lte(RULE.narrow.post).and(dN.gte(RULE.narrow.dnbr));
  var j = RULE.useNarrow ? A.or(B) : A;
  var valid = safeCount(LS.filterDate(win(y, W_POST_S), win(y, W_POST_E))).gt(0)
              .and(safeCount(LS.filterDate(win(y, W_PRE_S), win(y, W_PRE_E))).gt(0));
  return j.updateMask(valid).updateMask(WATER_FREE).updateMask(ROI_MASK)
          .rename('jhum').clip(roi);
}

function cleanJhum(y) {
  var j = rawJhum(y);
  if (RECOVERY_CHECK && y < END_YEAR) {                 // greening the next year
    var rec = postWide(y + 1).subtract(postWide(y)).gte(RECOVERY_MIN);
    j = j.and(rec.unmask(1));
  }
  if (MIN_PIXELS > 1) {                                  // minimum mapping unit
    var n = j.selfMask().connectedPixelCount(MIN_PIXELS + 1, true).rename('n');
    j = j.and(n.gte(MIN_PIXELS).unmask(0, false));
  }
  return j.updateMask(ROI_MASK).rename('jhum_' + y).toByte().clip(roi);
}

var annualList = YEARS.map(function (y) {
  return cleanJhum(y).unmask(0).updateMask(ROI_MASK).clip(roi).set('year', y); });
var stack = ee.Image.cat(annualList).updateMask(ROI_MASK).clip(roi);   // one band per year
var annual = ee.ImageCollection.fromImages(
  annualList.map(function (img, i) {
    return ee.Image(img).rename('jhum').set('year', YEARS[i]); }));

//--------------------------------------------------------------
// 3. TEMPORAL METRICS (fallow age, cycle length, frequency)
//--------------------------------------------------------------
var count = annual.sum().rename('jhum_count').toByte().updateMask(ROI_MASK).clip(roi);

// year images, masked to clearing years
var yearImgs = annual.map(function (img) {
  return ee.Image.constant(ee.Number(img.get('year'))).toInt16()
           .updateMask(img.eq(1)).rename('yr');
});
var lastClearing  = yearImgs.max().rename('last_clearing');
var firstClearing = yearImgs.min().rename('first_clearing');
var fallowAge     = ee.Image.constant(END_YEAR).subtract(lastClearing)
                      .rename('fallow_age').toInt16();

// --- cycle metrics use only the full-coverage years (METRIC_START onward) ---
var annualFull = annual.filter(ee.Filter.gte('year', METRIC_START));
var countFull  = annualFull.sum().rename('jhum_count_full').toByte();
var yearImgsFull = annualFull.map(function (img) {
  return ee.Image.constant(ee.Number(img.get('year'))).toInt16()
           .updateMask(img.eq(1)).rename('yr'); });
var lastFull  = yearImgsFull.max();
var firstFull = yearImgsFull.min();

// mean cycle length = (last - first) / (number of intervals)
var cycleMean = lastFull.subtract(firstFull)
                  .divide(countFull.subtract(1).max(1))
                  .updateMask(countFull.gte(2)).rename('cycle_mean').toFloat();

// most recent completed cycle = last clearing - previous clearing
var prevClearing = yearImgsFull.map(function (im) {
    return im.updateMask(im.lt(lastFull)); }).max().rename('prev_clearing');
var cycleLast = lastFull.subtract(prevClearing).rename('cycle_last').toInt16();

// All bands must share one data type for a GeoTIFF export -> cast to Float32.
// cycle_mean is fractional, so Float is the only common type.
var summary = ee.Image.cat([count, lastClearing, firstClearing,
                            fallowAge, cycleMean, cycleLast, countFull])
                .toFloat()
                .updateMask(ROI_MASK).clip(roi);

//--------------------------------------------------------------
// 4. DISPLAY
//--------------------------------------------------------------
Map.addLayer(cleanJhum(END_YEAR).selfMask(), {palette: ['#ff00ff']}, 'Jhum ' + END_YEAR);
Map.addLayer(fallowAge, {min: 0, max: 25,
  palette: ['#d73027', '#fc8d59', '#fee08b', '#d9ef8b', '#91cf60', '#1a9850', '#006837']},
  'Fallow age (years since last clearing)', false);
Map.addLayer(count, {min: 0, max: 6,
  palette: ['#ffffff', '#fee5d9', '#fcae91', '#fb6a4a', '#de2d26', '#a50f15']},
  'Number of clearings ' + START_YEAR + '-' + END_YEAR, false);
Map.addLayer(cycleMean, {min: 3, max: 15,
  palette: ['#b2182b', '#ef8a62', '#fddbc7', '#d1e5f0', '#67a9cf', '#2166ac']},
  'Mean cycle length (years)', false);

//--------------------------------------------------------------
// 5. AREA TABLE (one aggregation per year — safe for the Console)
//--------------------------------------------------------------
var px = ee.Image.pixelArea().divide(1e6);
function areaRow(y) {
  var j = cleanJhum(y);
  var raw = rawJhum(y);
  var img = ee.Image.cat([
    px.updateMask(j.selfMask()).rename('jhum_km2'),
    px.updateMask(raw.selfMask()).rename('jhum_km2_before_filter'),
    px.updateMask(j.mask()).rename('valid_km2')]);
  return ee.Feature(null, img.reduceRegion({
    reducer: ee.Reducer.sum(), geometry: roi.geometry(),
    scale: 30, maxPixels: 1e13, tileScale: 4})).set('year', y);
}
var areaTable = ee.FeatureCollection(YEARS.map(areaRow));
print('Area for ' + END_YEAR + ' (check before exporting the whole series)', areaRow(END_YEAR));

//--------------------------------------------------------------
// 6. EXPORTS
//--------------------------------------------------------------
if (EXPORT_STACK_ASSET) {
  Export.image.toAsset({
    image: stack, description: 'MIZ_Jhum_AnnualStack_' + START_YEAR + '_' + END_YEAR + '_' + VERSION,
    assetId: ASSET_DIR + 'MIZ_Jhum_AnnualStack_' + START_YEAR + '_' + END_YEAR + '_' + VERSION,
    region: roi.geometry(), scale: 30, crs: 'EPSG:32646', maxPixels: 1e13
  });
}
if (EXPORT_SUMMARY) {
  Export.image.toDrive({
    image: summary, description: 'MIZ_Jhum_FallowMetrics_' + START_YEAR + '_' + END_YEAR + '_' + VERSION,
    fileNamePrefix: 'MIZ_Jhum_FallowMetrics_' + START_YEAR + '_' + END_YEAR + '_' + VERSION,
    folder: DRIVE_DIR, region: roi.geometry(), scale: 30, crs: 'EPSG:32646',
    maxPixels: 1e13, fileFormat: 'GeoTIFF', formatOptions: {cloudOptimized: true}
  });
}
if (EXPORT_AREA_TABLE) {
  Export.table.toDrive({
    collection: areaTable,
    description: 'MIZ_Jhum_AnnualArea_' + START_YEAR + '_' + END_YEAR + '_' + VERSION,
    fileNamePrefix: 'MIZ_Jhum_AnnualArea_' + START_YEAR + '_' + END_YEAR + '_' + VERSION,
    folder: DRIVE_DIR, fileFormat: 'CSV',
    selectors: ['year', 'jhum_km2', 'jhum_km2_before_filter', 'valid_km2']
  });
}
