//==============================================================
// MIZORAM JHUM — PART G: RE-CALIBRATION FEATURE EXTRACTION
//
// The 400 verified check points showed the v3 map over-calls jhum
// (user's accuracy 0.44) and misses late/early burns (producer's 0.30).
// Two suspected causes:
//   (1) the "before" image uses the 75th percentile (greenest) while the
//       "after" image uses the 25th (most burned). That asymmetry creates
//       artificial dNBR on deciduous forest and bamboo that dries each
//       dry season -> false positives.
//   (2) the burn window (Mar–May) misses plots burned in Nov–Feb or
//       in late May–June -> false negatives.
//
// This script extracts SEVERAL composite variants at your 400 labelled
// points, so the equation can be refitted on verified data.
//
// BEFORE RUNNING: upload Label_2025_for_GEE.csv as a GEE table asset
//   Assets -> NEW -> CSV file, Advanced: X column = lon, Y column = lat,
//   CRS = EPSG:4326, asset name Label_2025
//==============================================================

//--------------------------------------------------------------
// 0. SETTINGS
//--------------------------------------------------------------
var POINTS_ASSET = 'projects/YOUR_PROJECT/assets/Label_2025';  // <-- your uploaded points
var YEAR         = 2025;
var DRIVE_DIR    = 'Jhum_recalibration_v5';
var VERSION      = 'v5';

// Composite variants tested: [preWindow, postWindow, prePct, postPct, label]
// Windows are [yearOffset, month, day] pairs.
var VARIANTS = [
  {id: 'v3',   preS: [-1, 3, 1], preE: [-1, 6, 1], postS: [0, 3, 1], postE: [0, 6, 1], pPre: 75, pPost: 25},
  {id: 'sym25', preS: [-1, 3, 1], preE: [-1, 6, 1], postS: [0, 3, 1], postE: [0, 6, 1], pPre: 25, pPost: 25},
  {id: 'sym50', preS: [-1, 3, 1], preE: [-1, 6, 1], postS: [0, 3, 1], postE: [0, 6, 1], pPre: 50, pPost: 50},
  {id: 'wide25', preS: [-2, 11, 1], preE: [-1, 6, 1], postS: [-1, 11, 1], postE: [0, 6, 1], pPre: 25, pPost: 25},
  {id: 'wide50', preS: [-2, 11, 1], preE: [-1, 6, 1], postS: [-1, 11, 1], postE: [0, 6, 1], pPre: 50, pPost: 50}
];

//--------------------------------------------------------------
// 1. HARMONIZED LANDSAT (same as Part A)
//--------------------------------------------------------------
var REBUILD_GEOMETRY = true;   // true = build points from the lon/lat columns
                               // (needed if the CSV was uploaded without X/Y set)

var pts = ee.FeatureCollection(POINTS_ASSET);
if (REBUILD_GEOMETRY) {
  pts = pts.map(function (f) {
    var x = ee.Number.parse(ee.Algorithms.String(f.get('lon')));
    var y = ee.Number.parse(ee.Algorithms.String(f.get('lat')));
    return ee.Feature(ee.Geometry.Point([x, y]), f.toDictionary());
  });
}
print('Points loaded', pts.size(), pts.first());
print('First geometry (must be a Point with coordinates)', pts.first().geometry());

var roi = ee.FeatureCollection('FAO/GAUL/2015/level1')
            .filter(ee.Filter.eq('ADM1_NAME', 'Mizoram'));

function prepL(img, nir, swir2, harmonize) {
  img = ee.Image(img);
  var qa = img.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(1 << 1).eq(0).and(qa.bitwiseAnd(1 << 3).eq(0))
             .and(qa.bitwiseAnd(1 << 4).eq(0));
  var out = img.select([nir, swir2], ['NIR', 'SWIR2']).multiply(0.0000275).add(-0.2)
               .updateMask(mask);
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
function pctImg(col, p) {
  return ee.Image(ee.Algorithms.If(col.size().gt(0),
           col.map(nbrOf).reduce(ee.Reducer.percentile([p])),
           ee.Image.constant(0).updateMask(0)));
}

//--------------------------------------------------------------
// 2. FEATURE IMAGE: all variants as separate bands
//--------------------------------------------------------------
var bands = [];
VARIANTS.forEach(function (v) {
  var pre  = pctImg(LS.filterDate(win(YEAR, v.preS),  win(YEAR, v.preE)),  v.pPre);
  var post = pctImg(LS.filterDate(win(YEAR, v.postS), win(YEAR, v.postE)), v.pPost);
  bands.push(pre.rename('pre_' + v.id));
  bands.push(post.rename('post_' + v.id));
  bands.push(pre.subtract(post).rename('dnbr_' + v.id));
});

// extra context: seasonal amplitude of the PREVIOUS year (how much this
// pixel naturally browns each dry season) and terrain
var prevSeason = LS.filterDate(win(YEAR, [-1, 11, 1]), win(YEAR, [0, 6, 1])).map(nbrOf);
bands.push(prevSeason.reduce(ee.Reducer.percentile([90]))
             .subtract(prevSeason.reduce(ee.Reducer.percentile([10])))
             .rename('season_amplitude'));
bands.push(ee.Terrain.products(ee.Image('USGS/SRTMGL1_003'))
             .select(['elevation', 'slope']));

var featureImage = ee.Image.cat(bands).toFloat();
print('Bands extracted', featureImage.bandNames());

//--------------------------------------------------------------
// 3. SAMPLE AT THE LABELLED POINTS AND EXPORT
//--------------------------------------------------------------
var keep = ['point_id', 'ref_jhum', 'ref_confidence', 'ref_landcover', 'notes'];
var table = featureImage.sampleRegions({
  collection: pts, properties: keep, scale: 30, tileScale: 8, geometries: false
});

var cols = keep.slice();
VARIANTS.forEach(function (v) {
  cols = cols.concat(['pre_' + v.id, 'post_' + v.id, 'dnbr_' + v.id]);
});
cols = cols.concat(['season_amplitude', 'elevation', 'slope']);

Export.table.toDrive({
  collection: table,
  description: 'MIZ_Jhum_G1_RecalibrationFeatures_' + YEAR + '_' + VERSION,
  fileNamePrefix: 'MIZ_Jhum_G1_RecalibrationFeatures_' + YEAR + '_' + VERSION,
  folder: DRIVE_DIR, fileFormat: 'CSV', selectors: cols
});

//--------------------------------------------------------------
// 4. QUICK LOOK: medians of each variant by reference class
//--------------------------------------------------------------
// works whether GEE read ref_jhum as a number or as text
var j = table.filter(ee.Filter.inList('ref_jhum', [1, '1']));
var n = table.filter(ee.Filter.inList('ref_jhum', [0, '0']));
var summary = ee.FeatureCollection(VARIANTS.map(function (v) {
  var med = function (c, b) { return c.aggregate_array(b).reduce(ee.Reducer.median()); };
  return ee.Feature(null, {
    variant: v.id,
    jhum_dNBR: med(j, 'dnbr_' + v.id), notjhum_dNBR: med(n, 'dnbr_' + v.id),
    jhum_post: med(j, 'post_' + v.id), notjhum_post: med(n, 'post_' + v.id)
  });
}));
print('Variant comparison (bigger gap between jhum and not-jhum = better)',
      ui.Chart.feature.byFeature(summary, 'variant',
        ['jhum_dNBR', 'notjhum_dNBR', 'jhum_post', 'notjhum_post']).setChartType('Table'));
