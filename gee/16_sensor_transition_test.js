//==============================================================
// MIZORAM JHUM — SENSOR TRANSITION TEST
//
// Why: the stable-forest check found a residual offset between
// harmonized ETM+ and OLI of about -0.067 NBR units, which is a
// third of the 0.20 dNBR threshold. If a year's "before" composite
// is dominated by one sensor and its "after" composite by another,
// that offset could shift the mapped area.
//
// Test: for the transition years, rebuild the map three ways and
// compare the area:
//   (1) all sensors, as published
//   (2) excluding the newly arrived sensor
//   (3) excluding the departing sensor
// A year is robust if the three areas agree within a few per cent.
//
// Transition years tested: 1999 and 2000 (ETM+ enters), 2012 and
// 2013 (TM leaves, OLI enters), 2022 and 2023 (OLI-2 enters).
//==============================================================

var ROI_NAME  = 'Mizoram';
var DRIVE_DIR = 'Jhum_sensor_transition';

var TESTS = [
  {year: 1999, sets: {all: ['L5','L7'], noNew: ['L5'],       noOld: ['L7']}},
  {year: 2000, sets: {all: ['L5','L7'], noNew: ['L5'],       noOld: ['L7']}},
  {year: 2012, sets: {all: ['L5','L7'], noNew: ['L5'],       noOld: ['L7']}},
  {year: 2013, sets: {all: ['L7','L8'], noNew: ['L7'],       noOld: ['L8']}},
  {year: 2014, sets: {all: ['L7','L8'], noNew: ['L7'],       noOld: ['L8']}},
  {year: 2022, sets: {all: ['L8','L9'], noNew: ['L8'],       noOld: ['L9']}},
  {year: 2023, sets: {all: ['L8','L9'], noNew: ['L8'],       noOld: ['L9']}}
];

var RULE = {wide: {post: 0.26, dnbr: 0.20}, narrow: {post: 0.14, dnbr: 0.46}};
var W_PRE_S = [-2, 11, 1], W_PRE_E = [-1, 6, 1];
var W_POST_S = [-1, 11, 1], W_POST_E = [0, 6, 1];
var N_PRE_S = [-1, 3, 1], N_PRE_E = [-1, 6, 1];
var N_POST_S = [0, 3, 1], N_POST_E = [0, 6, 1];
var PCT = 25;

var roi = ee.FeatureCollection('FAO/GAUL/2015/level1')
            .filter(ee.Filter.eq('ADM1_NAME', ROI_NAME));
var ROI_MASK = ee.Image.constant(1).clip(roi).mask();
var WATER_FREE = ee.Image('JRC/GSW1_4/GlobalSurfaceWater')
                   .select('occurrence').unmask(0).lt(50);
Map.centerObject(roi, 8);

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
var BASE = {
  L5: ee.ImageCollection('LANDSAT/LT05/C02/T1_L2').filterBounds(roi)
        .map(function (i) { return prep(i, 'SR_B4', 'SR_B7', true); }),
  L7: ee.ImageCollection('LANDSAT/LE07/C02/T1_L2').filterBounds(roi)
        .map(function (i) { return prep(i, 'SR_B4', 'SR_B7', true); }),
  L8: ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filterBounds(roi)
        .map(function (i) { return prep(i, 'SR_B5', 'SR_B7', false); }),
  L9: ee.ImageCollection('LANDSAT/LC09/C02/T1_L2').filterBounds(roi)
        .map(function (i) { return prep(i, 'SR_B5', 'SR_B7', false); })
};
function stack(list) {
  var col = BASE[list[0]];
  for (var i = 1; i < list.length; i++) col = col.merge(BASE[list[i]]);
  return col;
}
function win(y, w) { return ee.Date.fromYMD(y + w[0], w[1], w[2]); }
function pctImg(col) {
  return ee.Image(ee.Algorithms.If(col.size().gt(0),
    col.map(function (i) {
      return ee.Image(i).normalizedDifference(['NIR', 'SWIR2']).rename('NBR'); })
      .reduce(ee.Reducer.percentile([PCT])),
    ee.Image.constant(0).updateMask(0))).rename('NBR').toFloat();
}
function safeCount(col) {
  return ee.Image(ee.Algorithms.If(col.size().gt(0),
    col.select('NIR').count(), ee.Image.constant(0))).rename('n').unmask(0);
}
function jhum(y, list) {
  var LS = stack(list);
  var wPreC = LS.filterDate(win(y, W_PRE_S), win(y, W_PRE_E));
  var wPostC = LS.filterDate(win(y, W_POST_S), win(y, W_POST_E));
  var pW = pctImg(wPostC), dW = pctImg(wPreC).subtract(pW);
  var pN = pctImg(LS.filterDate(win(y, N_POST_S), win(y, N_POST_E)));
  var dN = pctImg(LS.filterDate(win(y, N_PRE_S), win(y, N_PRE_E))).subtract(pN);
  var A = pW.lte(RULE.wide.post).and(dW.gte(RULE.wide.dnbr));
  var B = pN.lte(RULE.narrow.post).and(dN.gte(RULE.narrow.dnbr));
  var valid = safeCount(wPreC).gt(0).and(safeCount(wPostC).gt(0));
  return A.or(B).updateMask(valid).updateMask(WATER_FREE).updateMask(ROI_MASK).clip(roi);
}

var px = ee.Image.pixelArea().divide(1e6);
var rows = [];
TESTS.forEach(function (t) {
  var keys = ['all', 'noNew', 'noOld'];
  var bands = keys.map(function (k) {
    return px.updateMask(jhum(t.year, t.sets[k]).selfMask()).rename('km2_' + k);
  });
  var valid = keys.map(function (k) {
    return px.updateMask(jhum(t.year, t.sets[k]).mask()).rename('valid_' + k);
  });
  rows.push(ee.Feature(null, ee.Image.cat(bands.concat(valid)).reduceRegion({
    reducer: ee.Reducer.sum(), geometry: roi.geometry(), scale: 30,
    maxPixels: 1e13, tileScale: 8})).set('year', t.year)
      .set('sensors_all', t.sets.all.join('+'))
      .set('sensors_noNew', t.sets.noNew.join('+'))
      .set('sensors_noOld', t.sets.noOld.join('+')));
});
var tbl = ee.FeatureCollection(rows);
print('First test year (check before exporting them all)', rows[0]);

Export.table.toDrive({
  collection: tbl, description: 'MIZ_Jhum_SensorTransitionTest',
  fileNamePrefix: 'MIZ_Jhum_SensorTransitionTest', folder: DRIVE_DIR, fileFormat: 'CSV',
  selectors: ['year', 'sensors_all', 'km2_all', 'valid_all',
              'sensors_noNew', 'km2_noNew', 'valid_noNew',
              'sensors_noOld', 'km2_noOld', 'valid_noOld']});
