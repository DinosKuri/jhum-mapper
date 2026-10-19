//==============================================================
// MIZORAM JHUM DETECTION — PART D: STRESS TEST
//
// Question: the rule "NBR_post <= 0.16" separates jhum from random
// background very well. But does it also (wrongly) flag land that is
// bare for OTHER reasons — settlements, riverbeds/landslides, fields?
// And does adding a change condition (dNBR >= X) fix that without
// losing real jhum?
//
// RULE 1 : NBR_post <= CUT_POST                     (state only)
// RULE 2 : NBR_post <= CUT_POST  AND  dNBR >= X     (state + change)
//==============================================================

//--------------------------------------------------------------
// 0. SETTINGS
//--------------------------------------------------------------
var SAMPLE_ASSET = 'projects/phdmrsei89448/assets/MIZ_Jhum_Samples_Landsat_1990_2025_v3';   // PRIVATE ASSET — see README, Private Earth Engine assets
var CUT_POST  = 0.10;   // <-- NBR_post leave-one-year-out threshold from Part C (v2 = 0.10; update after v3)
var X_LIST    = [0, 0.10, 0.20, 0.25, 0.30, 0.35, 0.40]; // dNBR change bars tested
var YEARS     = [1990, 2000, 2005, 2010, 2015, 2020, 2025];
var DRIVE_DIR = 'Jhum_tests_v3';   // GEE uses a folder NAME, not a path
var VERSION   = 'v3';

// Confuser land-cover from ESA WorldCover (only available for 2020 & 2021)
var TEST_YEARS = {2020: 'ESA/WorldCover/v100', 2025: 'ESA/WorldCover/v200'};
var CONFUSERS  = {50: 'BuiltUp', 60: 'Bare_riverbed_landslide', 40: 'Cropland'};
var N_PER_CLASS = 1000;

// Windows and percentiles — MUST match Part A  ([yearOffset, month, day])
var PRE_START = [-1, 3, 1], PRE_END = [-1, 6, 1], POST_START = [0, 3, 1], POST_END = [0, 6, 1];
var PRE_PCT = 75, POST_PCT = 25;
function win(y, w) { return ee.Date.fromYMD(y + w[0], w[1], w[2]); }

//--------------------------------------------------------------
// 1. SAME HARMONIZED LANDSAT AS PART A
//--------------------------------------------------------------
var roi = ee.FeatureCollection('FAO/GAUL/2015/level1')
            .filter(ee.Filter.eq('ADM1_NAME', 'Mizoram'));
Map.centerObject(roi, 8);

function prepL(img, nir, swir2) {
  img = ee.Image(img);
  var qa = img.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(1 << 1).eq(0).and(qa.bitwiseAnd(1 << 3).eq(0))
             .and(qa.bitwiseAnd(1 << 4).eq(0));
  return ee.Image(img.select([nir, swir2], ['NIR', 'SWIR2'])
            .multiply(0.0000275).add(-0.2).updateMask(mask).toFloat()
            .copyProperties(img, ['system:time_start']));
}
function toOLI(img) {
  img = ee.Image(img);
  return ee.Image(img.multiply([0.8462, 0.9071]).add([0.0412, 0.0172]).toFloat()
            .copyProperties(img, ['system:time_start']));
}
var LS = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2').filterBounds(roi)
           .map(function (i) { return toOLI(prepL(i, 'SR_B4', 'SR_B7')); })
  .merge(ee.ImageCollection('LANDSAT/LE07/C02/T1_L2').filterBounds(roi)
           .map(function (i) { return toOLI(prepL(i, 'SR_B4', 'SR_B7')); }))
  .merge(ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filterBounds(roi)
           .map(function (i) { return prepL(i, 'SR_B5', 'SR_B7'); }))
  .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2').filterBounds(roi)
           .map(function (i) { return prepL(i, 'SR_B5', 'SR_B7'); }));

function nbrPercentile(col, pct) {
  return col.map(function (i) {
    return ee.Image(i).normalizedDifference(['NIR', 'SWIR2']).rename('NBR').toFloat();
  }).reduce(ee.Reducer.percentile([pct]));
}
function yearFeatures(y) {
  var pre = nbrPercentile(LS.filterDate(win(y, PRE_START), win(y, PRE_END)), PRE_PCT)
              .rename('NBR_pre');
  var post = nbrPercentile(LS.filterDate(win(y, POST_START), win(y, POST_END)), POST_PCT)
              .rename('NBR_post');
  return ee.Image.cat([pre, post, pre.subtract(post).rename('dNBR')]);
}

//--------------------------------------------------------------
// 2. HOW MUCH REAL JHUM DO WE KEEP? (from your sample table)
//--------------------------------------------------------------
var S = ee.FeatureCollection(SAMPLE_ASSET);

var keep = ee.FeatureCollection(YEARS.map(function (y) {
  var j  = S.filter(ee.Filter.eq('year', y)).filter(ee.Filter.eq('jhum', 1));
  var bg = S.filter(ee.Filter.eq('year', y)).filter(ee.Filter.eq('jhum', 0));
  var props = {year_label: 'Y' + y};
  X_LIST.forEach(function (x) {
    var rule = ee.Filter.and(ee.Filter.lte('NBR_post', CUT_POST), ee.Filter.gte('dNBR', x));
    var k = 'dNBR_ge_' + String(x).replace('.', 'p');
    props['jhum_kept_' + k] = ee.Number(j.filter(rule).size()).divide(j.size()).multiply(100);
    props['bg_flagged_' + k] = ee.Number(bg.filter(rule).size()).divide(bg.size()).multiply(100);
  });
  return ee.Feature(null, props);
}));

var keepCols = X_LIST.map(function (x) { return 'jhum_kept_dNBR_ge_' + String(x).replace('.', 'p'); });
var bgCols   = X_LIST.map(function (x) { return 'bg_flagged_dNBR_ge_' + String(x).replace('.', 'p'); });

print('OUTPUT 1 – % of REAL JHUM kept (higher = better). dNBR_ge_0 = Rule 1 alone',
      ui.Chart.feature.byFeature(keep, 'year_label', keepCols).setChartType('Table'));
print('OUTPUT 2 – % of random BACKGROUND wrongly flagged (lower = better)',
      ui.Chart.feature.byFeature(keep, 'year_label', bgCols).setChartType('Table'));

//--------------------------------------------------------------
// 3. HOW MUCH CONFUSER LAND IS WRONGLY FLAGGED? (2020, 2025)
//--------------------------------------------------------------
var confRows = ee.FeatureCollection([]);
var confPoints = ee.FeatureCollection([]);   // raw confuser pixels, exported for Colab
Object.keys(TEST_YEARS).forEach(function (yStr) {
  var y  = Number(yStr);
  var lc = ee.ImageCollection(TEST_YEARS[yStr]).first().select('Map').rename('lc');
  var img = yearFeatures(y).addBands(lc);

  var rawPts = img.stratifiedSample({
    numPoints: 0, classBand: 'lc', region: roi.geometry(), scale: 30,
    classValues: Object.keys(CONFUSERS).map(Number),
    classPoints: Object.keys(CONFUSERS).map(function () { return N_PER_CLASS; }),
    seed: 7, tileScale: 8, geometries: false
  });
  var pts = rawPts.filter(ee.Filter.notNull(['NBR_post', 'dNBR']));
  confPoints = confPoints.merge(pts.map(function (f) { return f.set('year', y); }));

  Object.keys(CONFUSERS).forEach(function (code) {
    var c = pts.filter(ee.Filter.eq('lc', Number(code)));
    var props = {row: CONFUSERS[code] + '_Y' + y, n: c.size()};
    X_LIST.forEach(function (x) {
      var rule = ee.Filter.and(ee.Filter.lte('NBR_post', CUT_POST), ee.Filter.gte('dNBR', x));
      props['flagged_dNBR_ge_' + String(x).replace('.', 'p')] =
          ee.Number(c.filter(rule).size()).divide(c.size().max(1)).multiply(100);
    });
    confRows = confRows.merge(ee.FeatureCollection([ee.Feature(null, props)]));
  });
});

var confCols = ['n'].concat(X_LIST.map(function (x) {
  return 'flagged_dNBR_ge_' + String(x).replace('.', 'p'); }));
print('OUTPUT 3 – % of CONFUSER land wrongly called jhum (lower = better)',
      ui.Chart.feature.byFeature(confRows, 'row', confCols).setChartType('Table'));

//--------------------------------------------------------------
// 4. VISUAL CHECK (2025): Rule 1 vs Rule 2
//--------------------------------------------------------------
var f25 = yearFeatures(2025);
Map.addLayer(f25.select('NBR_post').lte(CUT_POST).selfMask().clip(roi),
             {palette: ['orange']}, 'Rule 1: NBR_post only (2025)', false);
Map.addLayer(f25.select('NBR_post').lte(CUT_POST).and(f25.select('dNBR').gte(0.20))
             .selfMask().clip(roi),
             {palette: ['red']}, 'Rule 2: NBR_post + dNBR>=0.20 (2025)');

//--------------------------------------------------------------
// 5. EXPORTS
//    D1, D2 = result tables (backup if the Console times out)
//    D3     = confuser pixels (built-up, bare, cropland) for the Colab notebook
//--------------------------------------------------------------
Export.table.toDrive({collection: keep, description: 'MIZ_Jhum_D1_JhumKept_' + VERSION,
  fileNamePrefix: 'MIZ_Jhum_D1_JhumKept_' + VERSION, folder: DRIVE_DIR, fileFormat: 'CSV'});
Export.table.toDrive({collection: confRows, description: 'MIZ_Jhum_D2_ConfuserFlagged_' + VERSION,
  fileNamePrefix: 'MIZ_Jhum_D2_ConfuserFlagged_' + VERSION, folder: DRIVE_DIR, fileFormat: 'CSV'});
Export.table.toDrive({collection: confPoints,
  description: 'MIZ_Jhum_D3_ConfuserSamples_WorldCover_2020_2025_' + VERSION,
  fileNamePrefix: 'MIZ_Jhum_D3_ConfuserSamples_WorldCover_2020_2025_' + VERSION,
  folder: DRIVE_DIR, fileFormat: 'CSV',
  selectors: ['year', 'lc', 'NBR_pre', 'NBR_post', 'dNBR']});
