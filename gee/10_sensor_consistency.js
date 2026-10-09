//==============================================================
// MIZORAM JHUM — CROSS-SENSOR NBR CONSISTENCY
//
// Reviewer item 10: show that harmonized NBR is comparable across
// Landsat 5 TM, 7 ETM+, 8 OLI and 9 OLI-2.
//
// Method: sample NBR over STABLE FOREST (no clearing in the record,
// no water, slope and elevation unrestricted) in years when two
// sensors overlap, then compare distributions.
//
//   TM  vs ETM+ : 1999–2011
//   ETM+ vs OLI : 2013–2021
//   OLI vs OLI-2: 2022–2025
//
// Outputs: per-sensor NBR statistics, paired differences (mean,
// median, RMSE) and a histogram overlap figure for the supplement.
//==============================================================

//--------------------------------------------------------------
// 0. SETTINGS
//--------------------------------------------------------------
var ROI_NAME  = 'Mizoram';
var DRIVE_DIR = 'Jhum_sensor_check';
var N_POINTS  = 2000;        // stable-forest sample points
var SEED      = 99;
var SEASON    = [[3, 1], [6, 1]];   // dry-season window, [month, day]

// Overlap periods: [startYear, endYear, sensorA, sensorB]
var PAIRS = [
  {a: 'L5', b: 'L7', y0: 1999, y1: 2011},
  {a: 'L7', b: 'L8', y0: 2013, y1: 2021},
  {a: 'L8', b: 'L9', y0: 2022, y1: 2025}
];

//--------------------------------------------------------------
// 1. AREA, STABLE-FOREST MASK
//--------------------------------------------------------------
var roi = ee.FeatureCollection('FAO/GAUL/2015/level1')
            .filter(ee.Filter.eq('ADM1_NAME', ROI_NAME));
Map.centerObject(roi, 8);

var gfc = ee.Image('UMD/hansen/global_forest_change_2023_v1_11');
var stable = gfc.select('treecover2000').gte(60)       // forest in 2000
               .and(gfc.select('loss').eq(0))          // never lost
               .and(gfc.select('gain').eq(0))          // never gained
               .and(ee.Image('JRC/GSW1_4/GlobalSurfaceWater')
                      .select('occurrence').unmask(0).lt(10));
Map.addLayer(stable.selfMask().clip(roi), {palette: ['#1b7837']}, 'stable forest', false);

var pts = stable.selfMask().stratifiedSample({
  numPoints: N_POINTS, classBand: 'treecover2000', region: roi.geometry(),
  scale: 30, seed: SEED, tileScale: 8, geometries: true
});
print('stable-forest sample points', pts.size());

//--------------------------------------------------------------
// 2. PER-SENSOR NBR
//--------------------------------------------------------------
function prep(img, nir, swir2, harmonize) {
  img = ee.Image(img);
  var qa = img.select('QA_PIXEL');
  var clear = qa.bitwiseAnd(1 << 1).eq(0).and(qa.bitwiseAnd(1 << 3).eq(0))
              .and(qa.bitwiseAnd(1 << 4).eq(0));
  var out = img.select([nir, swir2], ['NIR', 'SWIR2'])
               .multiply(0.0000275).add(-0.2).updateMask(clear);
  if (harmonize) out = out.multiply([0.8462, 0.9071]).add([0.0412, 0.0172]);
  return ee.Image(out.toFloat().copyProperties(img, ['system:time_start']));
}
var COL = {
  L5: ee.ImageCollection('LANDSAT/LT05/C02/T1_L2').filterBounds(roi)
        .map(function (i) { return prep(i, 'SR_B4', 'SR_B7', true); }),
  L7: ee.ImageCollection('LANDSAT/LE07/C02/T1_L2').filterBounds(roi)
        .map(function (i) { return prep(i, 'SR_B4', 'SR_B7', true); }),
  L8: ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filterBounds(roi)
        .map(function (i) { return prep(i, 'SR_B5', 'SR_B7', false); }),
  L9: ee.ImageCollection('LANDSAT/LC09/C02/T1_L2').filterBounds(roi)
        .map(function (i) { return prep(i, 'SR_B5', 'SR_B7', false); })
};
function nbrMedian(sensor, y0, y1) {
  var col = COL[sensor].filterDate(ee.Date.fromYMD(y0, SEASON[0][0], SEASON[0][1]),
                                   ee.Date.fromYMD(y1, SEASON[1][0], SEASON[1][1]));
  return ee.Image(ee.Algorithms.If(col.size().gt(0),
           col.map(function (i) {
             return ee.Image(i).normalizedDifference(['NIR', 'SWIR2']).rename('NBR'); })
             .median(),
           ee.Image.constant(0).updateMask(0))).rename('NBR_' + sensor).toFloat();
}

//--------------------------------------------------------------
// 3. PAIRED COMPARISON AT THE SAMPLE POINTS
//--------------------------------------------------------------
PAIRS.forEach(function (p) {
  var a = nbrMedian(p.a, p.y0, p.y1);
  var b = nbrMedian(p.b, p.y0, p.y1);
  var img = ee.Image.cat([a, b, a.subtract(b).rename('diff')]);

  var tbl = img.sampleRegions({collection: pts, scale: 30, tileScale: 8, geometries: false})
               .filter(ee.Filter.notNull(['diff']));

  var stats = tbl.reduceColumns({
    reducer: ee.Reducer.mean().combine(ee.Reducer.median(), '', true)
               .combine(ee.Reducer.stdDev(), '', true)
               .combine(ee.Reducer.count(), '', true),
    selectors: ['diff']});
  print(p.a + ' minus ' + p.b + ' over stable forest, ' + p.y0 + '-' + p.y1, stats);

  print(ui.Chart.feature.histogram(tbl, 'NBR_' + p.a, 40)
        .setOptions({title: 'NBR over stable forest, ' + p.a + ' (' + p.y0 + '-' + p.y1 + ')'}));
  print(ui.Chart.feature.histogram(tbl, 'NBR_' + p.b, 40)
        .setOptions({title: 'NBR over stable forest, ' + p.b + ' (' + p.y0 + '-' + p.y1 + ')'}));

  Export.table.toDrive({
    collection: tbl,
    description: 'MIZ_SensorCheck_' + p.a + '_vs_' + p.b,
    fileNamePrefix: 'MIZ_SensorCheck_' + p.a + '_vs_' + p.b,
    folder: DRIVE_DIR, fileFormat: 'CSV',
    selectors: ['NBR_' + p.a, 'NBR_' + p.b, 'diff']});
});

print('Report mean and median difference, standard deviation and RMSE per pair.',
      'RMSE = sqrt(mean(diff^2)); compute it from the exported CSV or in the notebook.');
