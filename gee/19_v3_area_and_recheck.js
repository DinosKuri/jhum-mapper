//==============================================================
// v3 PRELIMINARY MAP — AREA RECOVERY AND CHECK-POINT RE-RUN
//
// PURPOSE
//   The v3 map is the preliminary 2025 map that was rejected and
//   replaced by v6. It appears in the paper only in Table 4 and
//   Table S1. Two things about it are missing:
//
//     1. Its mapped AREA was never retained. Without the area the
//        producer's accuracy cannot be estimated, which is why
//        Table 4 carries a dash in that cell.
//     2. Its map class at the 400 check points lives only in
//        MIZ_Jhum_E2_CheckPoints_KEY_2025_M3_NBRpost_dNBR_v3.csv,
//        with no second copy.
//
//   This script rebuilds the v3 map exactly and fixes both. It
//   needs NO new interpretation: the 400 reference labels already
//   exist. Run it, read two numbers off the Console, and the dash
//   in Table 4 can be filled.
//
// WHAT IT PRINTS
//   A  mapped area of the v3 map, km2           -> fills Table 4
//   B  unmapped area, km2
//   C  the v3 class at each of the 400 points   -> verifies the
//      'stratum' column of the E2 key, point by point
//
// THE v3 RULE (recovered from gee/01, gee/11 and the E2 key; the
// logistic coefficients below reproduce the stored P_jhum_pct of
// all 400 check points to within 0.015, and the cut of 0.2164
// falls exactly in the gap between the two classes, P = 0.210 to
// 0.220, so this is the rule that made that file)
//   windows  1 Mar - 31 May, against the same months one year
//            earlier  (ANNIVERSARY, but ASYMMETRIC percentiles)
//   before   75th percentile of per-scene NBR   <- the greenest
//   after    25th percentile of per-scene NBR   <- the darkest
//   jhum if  P(jhum) >= 0.2164  AND  dNBR >= 0.30
//
// DO NOT use this rule for anything you publish. It is the design
// the paper rejects. The published rule is in gee/00 and gee/05.
//==============================================================

//--------------------------------------------------------------
// 0. SETTINGS
//--------------------------------------------------------------
var YEAR = 2025;

// Leave '' to only recover the area (section 4). To also re-run the
// check points (section 5), upload calibration_points_400_2025.csv
// from the dataset record as a table asset and put its id here.
var POINTS_ASSET = '';            // e.g. 'projects/xxx/assets/calibration_points_400_2025'

var DRIVE_DIR = 'MIZ_Jhum_v3_recheck';

// The v3 equation, as fitted in Colab and used in gee/01.
var EQ = {
  b0: -3.484372,
  coef: {NBR_post: -11.261907, dNBR: 9.137461},
  cut: 0.2164,
  min_dNBR: 0.30
};

// v3 windows and percentiles. The asymmetry is the point: it is
// what the paper shows manufactures change on seasonal surfaces.
var PRE_START = [-1, 3, 1], PRE_END = [-1, 6, 1];
var POST_START = [0, 3, 1], POST_END = [0, 6, 1];
var PRE_PCT = 75, POST_PCT = 25;

//--------------------------------------------------------------
// 1. ROI AND HARMONIZED LANDSAT  (identical to gee/00 and gee/04)
//--------------------------------------------------------------
var roi = ee.FeatureCollection('FAO/GAUL/2015/level1')
            .filter(ee.Filter.eq('ADM1_NAME', 'Mizoram'));

function prepL(img, nir, swir2, harmonize) {
  img = ee.Image(img);
  var qa = img.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(1 << 1).eq(0).and(qa.bitwiseAnd(1 << 3).eq(0))
             .and(qa.bitwiseAnd(1 << 4).eq(0));
  var out = img.select([nir, swir2], ['NIR', 'SWIR2'])
               .multiply(0.0000275).add(-0.2).updateMask(mask);
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

var WATER = ee.Image('JRC/GSW1_4/GlobalSurfaceWater').select('occurrence');
var WATER_FREE = WATER.unmask(0).lt(50);

function win(y, w) { return ee.Date.fromYMD(y + w[0], w[1], w[2]); }
function nbrOf(i) { return ee.Image(i).normalizedDifference(['NIR', 'SWIR2']).rename('NBR'); }
function pctImg(col, p) {
  return ee.Image(ee.Algorithms.If(col.size().gt(0),
           col.map(nbrOf).reduce(ee.Reducer.percentile([p])),
           ee.Image.constant(0).updateMask(0)));
}

//--------------------------------------------------------------
// 2. THE v3 MAP
//--------------------------------------------------------------
var preCol  = LS.filterDate(win(YEAR, PRE_START),  win(YEAR, PRE_END));
var postCol = LS.filterDate(win(YEAR, POST_START), win(YEAR, POST_END));
print('scenes in the before window', preCol.size());
print('scenes in the after window',  postCol.size());

var pre  = pctImg(preCol,  PRE_PCT ).rename('NBR_pre');
var post = pctImg(postCol, POST_PCT).rename('NBR_post');
var dnbr = pre.subtract(post).rename('dNBR');

var z = ee.Image.constant(EQ.b0)
          .add(post.multiply(EQ.coef.NBR_post))
          .add(dnbr.multiply(EQ.coef.dNBR));
var P = z.multiply(-1).exp().add(1).pow(-1).rename('P_jhum');

var valid = preCol.size().gt(0) ? pre.mask().and(post.mask()) : pre.mask();
var jhumV3 = P.gte(EQ.cut).and(dnbr.gte(EQ.min_dNBR))
              .updateMask(valid).updateMask(WATER_FREE)
              .rename('jhum_v3').clip(roi);

Map.centerObject(roi, 8);
Map.addLayer(jhumV3.selfMask(), {palette: ['#b2182b']}, 'v3 map ' + YEAR);

//--------------------------------------------------------------
// 3. THE NUMBER THAT FILLS TABLE 4
//--------------------------------------------------------------
var PROJ = ee.Projection('EPSG:32646').atScale(30);
var PX   = ee.Image.pixelArea().divide(1e6);

var areaMapped = PX.updateMask(jhumV3.eq(1)).reduceRegion({
  reducer: ee.Reducer.sum(), geometry: roi.geometry(),
  crs: 'EPSG:32646', scale: 30, maxPixels: 1e13, tileScale: 4}).get('area');

var areaFrame = PX.updateMask(jhumV3.mask()).reduceRegion({
  reducer: ee.Reducer.sum(), geometry: roi.geometry(),
  crs: 'EPSG:32646', scale: 30, maxPixels: 1e13, tileScale: 4}).get('area');

print('A  v3 MAPPED AREA, km2   (stratum 1)', areaMapped);
print('B  v3 FRAME AREA, km2    (mapped + unmapped, valid pixels)', areaFrame);
print('   unmapped area = B - A');
print('For scale, the published v6 map is 528.2 km2 on the same grid.');

//--------------------------------------------------------------
// 4. HOW TO FILL THE DASH ONCE YOU HAVE A AND B
//--------------------------------------------------------------
print(
 'PRODUCER\'S ACCURACY OF v3, from the areas above and the labels\n' +
 'you already have (no new interpretation needed):\n' +
 '   W1 = A / B                 weight of the mapped stratum\n' +
 '   W0 = (B - A) / B           weight of the unmapped stratum\n' +
 '   p1 = 50 / 100  = 0.500     reference jhum in the mapped stratum\n' +
 '   p0 =  2 / 300  = 0.00667   reference jhum in the unmapped stratum\n' +
 '\n' +
 '   TP        = W1 * p1\n' +
 '   omission  = W0 * p0\n' +
 '   producer\'s accuracy = TP / (TP + omission)\n' +
 '   user\'s accuracy     = p1 = 0.50   (already known)\n' +
 '\n' +
 'p1 and p0 come from joining the E2 check-point key to\n' +
 'calibration_points_400_2025.csv on point_id.');

//--------------------------------------------------------------
// 5. OPTIONAL — re-run the class at the 400 existing check points
//    This verifies the stratum column of the E2 key. It needs no
//    new labelling: the labels are already in the dataset record.
//--------------------------------------------------------------
if (POINTS_ASSET) {
  var pts = ee.FeatureCollection(POINTS_ASSET).map(function (f) {
    var x = ee.Number.parse(ee.Algorithms.String(f.get('lon')));
    var y = ee.Number.parse(ee.Algorithms.String(f.get('lat')));
    return ee.Feature(ee.Geometry.Point([x, y]), f.toDictionary());
  });
  print('check points loaded', pts.size());

  var sampled = ee.Image.cat([jhumV3.unmask(0).rename('map_v3'), P, post, dnbr])
                  .reduceRegions({collection: pts, reducer: ee.Reducer.first(), scale: 30});

  Export.table.toDrive({
    collection: sampled,
    description: 'MIZ_Jhum_v3_recheck_400pts_' + YEAR,
    fileNamePrefix: 'MIZ_Jhum_v3_recheck_400pts_' + YEAR,
    folder: DRIVE_DIR, fileFormat: 'CSV',
    selectors: ['point_id', 'ref_jhum', 'map_v3', 'P_jhum', 'NBR_post', 'dNBR']});

  print('Export queued. When it finishes, check that map_v3 reproduces\n' +
        'the stratum column of the E2 key: 100 ones and 300 zeros, and\n' +
        '50 of the ones labelled ref_jhum = 1.');
} else {
  print('POINTS_ASSET is empty, so section 5 is skipped. The area in\n' +
        'section 3 is the part that fills Table 4; the re-run is only a\n' +
        'cross-check of a file you already have.');
}

//--------------------------------------------------------------
// 6. NOTE ON DRAWING A FRESH SAMPLE
//--------------------------------------------------------------
print(
 'A FRESH SAMPLE OF THE v3 MAP IS NOT RECOMMENDED.\n' +
 'The existing 100 mapped-jhum check points are already a valid\n' +
 'probability sample of that map, interpreted blind on PlanetScope,\n' +
 'and they give a user\'s accuracy of 0.50. A new sample would cost\n' +
 'roughly 100 fresh interpretations and would change one cell of one\n' +
 'table describing a design the paper rejects. The published map is\n' +
 'unaffected either way: it is assessed on the 520-point dual-map\n' +
 'sample and reports 0.672 and 0.614.');
