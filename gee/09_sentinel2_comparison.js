//==============================================================
// MIZORAM JHUM — PART J: SENTINEL-2 TEST (2025)
//
// Repeats the v6 design with Sentinel-2 instead of Landsat, at the
// SAME labelled points, so the two sensors can be compared directly:
//   Landsat 30 m  vs  Sentinel-2 20 m (B8A / B12)
//
// Same compositing logic as v6:
//   wide   : 1 Nov – 31 May, and the SAME months one year earlier
//   narrow : 1 Mar – 31 May, and the same months one year earlier
//   25th percentile of per-image NBR on BOTH sides
//
// Cloud masking uses Cloud Score+ (QA60 is unreliable after 2022).
//
// OUTPUT: one CSV per point set, with S2 features next to your labels.
// Thresholds are then fitted on the CALIBRATION set and tested on the
// VALIDATION set, exactly as for Landsat.
//==============================================================

//--------------------------------------------------------------
// 0. SETTINGS
//--------------------------------------------------------------
var YEAR      = 2025;
var DRIVE_DIR = 'Jhum_sentinel2_v7';
var VERSION   = 'v7s2c';
var ROI_NAME  = 'Mizoram';

// Point sets to sample. Upload each labelled CSV as a GEE table asset
// (columns lon, lat, ref_jhum, …). Geometry is rebuilt from lon/lat.
var POINT_SETS = [
  {tag: 'calib400', asset: 'projects/YOUR_PROJECT/assets/Label_2025'},
  {tag: 'valid200', asset: 'projects/YOUR_PROJECT/assets/Blind_25'},
  {tag: 'valid_far40', asset: 'projects/YOUR_PROJECT/assets/XI_BLIND_25'}
  // add the 70 near-stratum points once uploaded:
  // ,{tag: 'valid_near70', asset: 'projects/YOUR_PROJECT/assets/N1_BLIND_2025'}
];

// windows and percentile — identical to the Landsat v6 rule
var W_PRE_S = [-2, 11, 1], W_PRE_E = [-1, 6, 1];
var W_POST_S = [-1, 11, 1], W_POST_E = [0, 6, 1];
var N_PRE_S = [-1, 3, 1], N_PRE_E = [-1, 6, 1];
var N_POST_S = [0, 3, 1], N_POST_E = [0, 6, 1];
var PCT = 25;

var CS_THRESH   = 0.60;   // Cloud Score+ cs_cdf; higher = stricter
var SAMPLE_SCALE = 20;    // B8A and B12 native resolution

//--------------------------------------------------------------
// 1. SENTINEL-2, SEVERAL CLOUD/COMPOSITE VARIANTS
//    Round 1 showed Sentinel-2 performing WORSE than Landsat, with
//    "after" composites too green on known jhum plots. The likely
//    cause is Cloud Score+ masking smoke, haze and dark burned
//    surfaces during the burning season. These variants test that:
//      cs06_p25 : cs_cdf >= 0.60, 25th percentile   (round 1)
//      cs40_p25 : cs_cdf >= 0.40, 25th percentile   (looser mask)
//      cs06_p10 : cs_cdf >= 0.60, 10th percentile   (deeper minimum)
//      cs_p25   : 'cs' band >= 0.50, 25th percentile (different score)
//--------------------------------------------------------------
var roi = ee.FeatureCollection('FAO/GAUL/2015/level1')
            .filter(ee.Filter.eq('ADM1_NAME', ROI_NAME));
Map.centerObject(roi, 9);

var csPlus = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED');
var S2raw = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
              .filterBounds(roi)
              .linkCollection(csPlus, ['cs', 'cs_cdf']);

function masked(band, thresh) {
  return S2raw.map(function (img) {
    var clear = img.select(band).gte(thresh);
    return img.select(['B8A', 'B12'], ['NIR', 'SWIR2'])
              .divide(10000).updateMask(clear).toFloat()
              .copyProperties(img, ['system:time_start']);
  });
}

function win(y, w) { return ee.Date.fromYMD(y + w[0], w[1], w[2]); }
function nbrOf(i) { return ee.Image(i).normalizedDifference(['NIR', 'SWIR2']).rename('NBR'); }
function pctImg(col, pct) {
  return ee.Image(ee.Algorithms.If(col.size().gt(0),
           col.map(nbrOf).reduce(ee.Reducer.percentile([pct])),
           ee.Image.constant(0).updateMask(0))).rename('NBR').toFloat();
}

var VARIANTS = [
  {id: 'cs06_p25', band: 'cs_cdf', th: 0.60, pct: 25},
  {id: 'cs40_p25', band: 'cs_cdf', th: 0.40, pct: 25},
  {id: 'cs06_p10', band: 'cs_cdf', th: 0.60, pct: 10},
  {id: 'cs_p25',   band: 'cs',     th: 0.50, pct: 25}
];

var bands = [];
VARIANTS.forEach(function (v) {
  var col = masked(v.band, v.th);
  var wPre  = pctImg(col.filterDate(win(YEAR, W_PRE_S),  win(YEAR, W_PRE_E)),  v.pct);
  var wPost = pctImg(col.filterDate(win(YEAR, W_POST_S), win(YEAR, W_POST_E)), v.pct);
  var nPre  = pctImg(col.filterDate(win(YEAR, N_PRE_S),  win(YEAR, N_PRE_E)),  v.pct);
  var nPost = pctImg(col.filterDate(win(YEAR, N_POST_S), win(YEAR, N_POST_E)), v.pct);
  bands.push(wPre.rename('pre_wide_' + v.id));          // kept: needed to test
  bands.push(wPost.rename('post_wide_' + v.id));        // whether the BEFORE image
  bands.push(wPre.subtract(wPost).rename('dnbr_wide_' + v.id));   // is too dark
  bands.push(nPre.rename('pre_narrow_' + v.id));
  bands.push(nPost.rename('post_narrow_' + v.id));
  bands.push(nPre.subtract(nPost).rename('dnbr_narrow_' + v.id));
  if (v.id === 'cs40_p25') {
    // how many clear scenes feed each composite (Landsat typically has ~5x fewer)
    var cnt = function (col) {
      return ee.Image(ee.Algorithms.If(col.size().gt(0),
               col.select('NIR').count(), ee.Image.constant(0))).rename('n').unmask(0);
    };
    bands.push(cnt(col.filterDate(win(YEAR, W_PRE_S),  win(YEAR, W_PRE_E))).rename('n_pre_wide'));
    bands.push(cnt(col.filterDate(win(YEAR, W_POST_S), win(YEAR, W_POST_E))).rename('n_post_wide'));
  }
  if (v.id === 'cs06_p25') {
    print('scenes kept (cs_cdf>=0.6) — wide before / wide after / narrow after',
          col.filterDate(win(YEAR, W_PRE_S), win(YEAR, W_PRE_E)).size(),
          col.filterDate(win(YEAR, W_POST_S), win(YEAR, W_POST_E)).size(),
          col.filterDate(win(YEAR, N_POST_S), win(YEAR, N_POST_E)).size());
  }
});
var features = ee.Image.cat(bands).toFloat();
print('S2 feature bands', features.bandNames());

//--------------------------------------------------------------
// 2. SAMPLE EACH POINT SET
//--------------------------------------------------------------
function rebuild(fc) {
  return fc.map(function (f) {
    var x = ee.Number.parse(ee.Algorithms.String(f.get('lon')));
    var y = ee.Number.parse(ee.Algorithms.String(f.get('lat')));
    return ee.Feature(ee.Geometry.Point([x, y]), f.toDictionary());
  });
}

var KEEP = ['point_id', 'ref_jhum', 'ref_confidence', 'ref_landcover', 'notes'];
var COLS = KEEP.slice();
VARIANTS.forEach(function (v) {
  COLS = COLS.concat(['pre_wide_' + v.id, 'post_wide_' + v.id, 'dnbr_wide_' + v.id,
                      'pre_narrow_' + v.id, 'post_narrow_' + v.id, 'dnbr_narrow_' + v.id]);
});
COLS = COLS.concat(['n_pre_wide', 'n_post_wide']);

POINT_SETS.forEach(function (ps) {
  var pts = rebuild(ee.FeatureCollection(ps.asset));
  var tbl = features.sampleRegions({
    collection: pts, properties: KEEP, scale: SAMPLE_SCALE,
    tileScale: 8, geometries: false
  });
  print(ps.tag + ': points sampled', tbl.size());
  Export.table.toDrive({
    collection: tbl,
    description: 'MIZ_Jhum_J1_S2features_' + ps.tag + '_' + YEAR + '_' + VERSION,
    fileNamePrefix: 'MIZ_Jhum_J1_S2features_' + ps.tag + '_' + YEAR + '_' + VERSION,
    folder: DRIVE_DIR, fileFormat: 'CSV', selectors: COLS
  });
});

//--------------------------------------------------------------
// 3. QUICK LOOK (first variant only)
//--------------------------------------------------------------
Map.addLayer(features.select('dnbr_wide_cs06_p25').clip(roi),
  {min: -0.1, max: 0.8, palette: ['#1a9850', '#ffffbf', '#d73027']}, 'S2 dNBR wide (cs06_p25)', false);
print('Thresholds must be refitted per variant on the calibration points.');
