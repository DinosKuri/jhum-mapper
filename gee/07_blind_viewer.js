//==============================================================
// MIZORAM JHUM — BLIND VIEWER (use this while labelling)
//
// Shows ONLY Landsat imagery: true colour and SWIR, for the "before"
// and "after" dry seasons of a chosen year.
//
// It deliberately contains NO jhum map, NO index layers, NO probability
// and NO click inspector. Seeing any of those during labelling destroys
// the independence of the reference data.
//
// Inspect map and index values (e.g. in 01_final_map.js, click inspector)
// only AFTER all labelling is finished.
//==============================================================

//--------------------------------------------------------------
// 0. SETTINGS
//--------------------------------------------------------------
var YEAR     = 2010;                     // year being labelled
var ROI_NAME = 'Mizoram';
var YEARS    = [1990, 2000, 2005, 2010, 2015, 2018, 2020, 2022, 2025];

// dry-season windows (same as the mapping rule, but used here only to
// build pictures, never to classify anything)
var PRE_S  = [-1, 11, 1], PRE_E  = [0, 6, 1];   // previous dry season
var POST_S = [0, 11, 1],  POST_E = [1, 6, 1];   // dry season of YEAR

//--------------------------------------------------------------
// 1. AREA AND HARMONIZED LANDSAT (6 bands, for true colour)
//--------------------------------------------------------------
var roi = ee.FeatureCollection('FAO/GAUL/2015/level1')
            .filter(ee.Filter.eq('ADM1_NAME', ROI_NAME));
Map.centerObject(roi, 8);

var NAMES  = ['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2'];
var SLOPES = [0.8474, 0.8483, 0.9047, 0.8462, 0.8937, 0.9071];
var INTER  = [0.0003, 0.0088, 0.0061, 0.0412, 0.0254, 0.0172];
var OLD = ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7'];
var NEW = ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7'];

function prep(img, bands, harmonize) {
  img = ee.Image(img);
  var qa = img.select('QA_PIXEL');
  var clear = qa.bitwiseAnd(1 << 1).eq(0)
              .and(qa.bitwiseAnd(1 << 3).eq(0))
              .and(qa.bitwiseAnd(1 << 4).eq(0));
  var out = img.select(bands, NAMES).multiply(0.0000275).add(-0.2).updateMask(clear);
  if (harmonize) out = out.multiply(SLOPES).add(INTER);
  return ee.Image(out.toFloat().copyProperties(img, ['system:time_start']));
}
var LS = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2').filterBounds(roi)
           .map(function (i) { return prep(i, OLD, true); })
  .merge(ee.ImageCollection('LANDSAT/LE07/C02/T1_L2').filterBounds(roi)
           .map(function (i) { return prep(i, OLD, true); }))
  .merge(ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filterBounds(roi)
           .map(function (i) { return prep(i, NEW, false); }))
  .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2').filterBounds(roi)
           .map(function (i) { return prep(i, NEW, false); }));

function win(y, w) { return ee.Date.fromYMD(y + w[0], w[1], w[2]); }

//--------------------------------------------------------------
// 2. PICTURES ONLY
//--------------------------------------------------------------
var TRUE_VIS = {bands: ['RED', 'GREEN', 'BLUE'], min: 0.02, max: 0.25, gamma: 1.3};
var BURN_VIS = {bands: ['SWIR2', 'NIR', 'RED'],  min: 0.02, max: 0.40, gamma: 1.1};

function draw(y) {
  Map.layers().reset();
  var before = LS.filterDate(win(y, PRE_S),  win(y, PRE_E)).median().clip(roi);
  var after  = LS.filterDate(win(y, POST_S), win(y, POST_E)).median().clip(roi);
  Map.addLayer(before, TRUE_VIS, (y - 1) + '/' + y + ' true colour — BEFORE', false);
  Map.addLayer(after,  TRUE_VIS, y + '/' + (y + 1) + ' true colour — AFTER', false);
  Map.addLayer(before, BURN_VIS, (y - 1) + '/' + y + ' SWIR — BEFORE', false);
  Map.addLayer(after,  BURN_VIS, y + '/' + (y + 1) + ' SWIR — AFTER (burns = orange-red)', true);
  status.setValue('Year ' + y + '. Tick one layer at a time and compare BEFORE with AFTER.');
}

//--------------------------------------------------------------
// 3. MINIMAL PANEL (year selector only)
//--------------------------------------------------------------
var status = ui.Label('loading…');
var sel = ui.Select({items: YEARS.map(String), value: String(YEAR),
                     onChange: function (v) { draw(Number(v)); }});
Map.add(ui.Panel({
  widgets: [ui.Label('Blind viewer — imagery only', {fontWeight: 'bold'}),
            ui.Panel([ui.Label('Year:'), sel], ui.Panel.Layout.flow('horizontal')),
            status,
            ui.Label('No map, no indices, no pixel values are shown here. That is deliberate.',
                     {fontSize: '11px', color: '555555'})],
  style: {position: 'top-left', width: '320px', padding: '8px'}}));

draw(YEAR);
