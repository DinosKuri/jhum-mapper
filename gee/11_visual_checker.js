//==============================================================
// MIZORAM JHUM — PART F: VISUAL CHECKER
//
// Shows, for any year:
//   - TRUE COLOUR Landsat before (previous-year Mar–May) and after (Mar–May)
//   - FALSE COLOUR (SWIR2/NIR/Red) where burns are unmistakable
//   - NBR_post, dNBR, P(jhum) and the final jhum map
//   - Click anywhere: prints all values and the equation result for that pixel
//
// Same harmonization, windows and equation as Part A/B (v3).
//==============================================================

//--------------------------------------------------------------
// 0. SETTINGS
//--------------------------------------------------------------
var YEAR = 2025;                                  // start year (changeable in the panel)
var POLYGON_ASSET = 'projects/phdmrsei89448/assets/Jhum_1990toJhum_2025';   // PRIVATE ASSET — see README, Private Earth Engine assets
var CHECKPOINT_ASSET = '';                        // optional: asset of your check points, else ''
var YEARS = [1990, 2000, 2005, 2010, 2015, 2018, 2020, 2022, 2025];

var EQ = {b0: -3.484372, coef: {NBR_post: -11.261907, dNBR: 9.137461},
          cut: 0.2164, min_dNBR: 0.30};

var PRE_START = [-1, 3, 1], PRE_END = [-1, 6, 1];   // "before": same season, previous year
var POST_START = [0, 3, 1], POST_END = [0, 6, 1];   // "after": burn season
var PRE_PCT = 75, POST_PCT = 25;

//--------------------------------------------------------------
// 1. AREA AND HARMONIZED LANDSAT (all 6 bands, for true colour)
//--------------------------------------------------------------
var roi = ee.FeatureCollection('FAO/GAUL/2015/level1')
            .filter(ee.Filter.eq('ADM1_NAME', 'Mizoram'));
var ref = ee.FeatureCollection(POLYGON_ASSET);
var WATER_FREE = ee.Image('JRC/GSW1_4/GlobalSurfaceWater')
                   .select('occurrence').unmask(0).lt(50);

var NAMES = ['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2'];
// Roy et al. (2016) TM/ETM+ -> OLI
var SLOPES = [0.8474, 0.8483, 0.9047, 0.8462, 0.8937, 0.9071];
var INTER  = [0.0003, 0.0088, 0.0061, 0.0412, 0.0254, 0.0172];

function prepL(img, bands, harmonize) {
  img = ee.Image(img);
  var qa = img.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(1 << 1).eq(0)
             .and(qa.bitwiseAnd(1 << 3).eq(0))
             .and(qa.bitwiseAnd(1 << 4).eq(0));
  var out = img.select(bands, NAMES).multiply(0.0000275).add(-0.2).updateMask(mask);
  if (harmonize) { out = out.multiply(SLOPES).add(INTER); }
  return ee.Image(out.toFloat().copyProperties(img, ['system:time_start']));
}

var OLD = ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7'];   // TM / ETM+
var NEW = ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7'];   // OLI / OLI-2

var LS = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2').filterBounds(roi)
           .map(function (i) { return prepL(i, OLD, true); })
  .merge(ee.ImageCollection('LANDSAT/LE07/C02/T1_L2').filterBounds(roi)
           .map(function (i) { return prepL(i, OLD, true); }))
  .merge(ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filterBounds(roi)
           .map(function (i) { return prepL(i, NEW, false); }))
  .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2').filterBounds(roi)
           .map(function (i) { return prepL(i, NEW, false); }));

function win(y, w) { return ee.Date.fromYMD(y + w[0], w[1], w[2]); }
function preCol(y)  { return LS.filterDate(win(y, PRE_START),  win(y, PRE_END)); }
function postCol(y) { return LS.filterDate(win(y, POST_START), win(y, POST_END)); }

//--------------------------------------------------------------
// 2. COMPOSITES, INDICES AND THE JHUM MAP FOR ONE YEAR
//--------------------------------------------------------------
function nbrOf(img) { return img.normalizedDifference(['NIR', 'SWIR2']).rename('NBR'); }

function build(y) {
  var pc = preCol(y), qc = postCol(y);

  // true-colour / false-colour composites (median = a normal-looking picture)
  var preRGB  = pc.median().clip(roi);
  var postRGB = qc.median().clip(roi);

  // index composites (percentiles, exactly as the equation uses them)
  var pre  = ee.Image(ee.Algorithms.If(pc.size().gt(0),
               pc.map(nbrOf).reduce(ee.Reducer.percentile([PRE_PCT])),
               ee.Image.constant(0).updateMask(0))).rename('NBR_pre');
  var post = ee.Image(ee.Algorithms.If(qc.size().gt(0),
               qc.map(nbrOf).reduce(ee.Reducer.percentile([POST_PCT])),
               ee.Image.constant(0).updateMask(0))).rename('NBR_post');
  var dnbr = pre.subtract(post).rename('dNBR');

  var z = ee.Image.constant(EQ.b0)
            .add(post.multiply(EQ.coef.NBR_post))
            .add(dnbr.multiply(EQ.coef.dNBR));
  var P = z.multiply(-1).exp().add(1).pow(-1).rename('P_jhum');
  var valid = pc.select('NIR').count().unmask(0).gt(0)
              .and(qc.select('NIR').count().unmask(0).gt(0));
  var jhum = P.gte(EQ.cut).and(dnbr.gte(EQ.min_dNBR))
              .updateMask(valid).updateMask(WATER_FREE).rename('Jhum');

  return {preRGB: preRGB, postRGB: postRGB,
          idx: ee.Image.cat([pre, post, dnbr, P]).clip(roi),
          jhum: jhum.clip(roi), nScenes: [pc.size(), qc.size()]};
}

//--------------------------------------------------------------
// 3. VISUALISATION SETTINGS
//--------------------------------------------------------------
var TRUE_VIS  = {bands: ['RED', 'GREEN', 'BLUE'], min: 0.02, max: 0.25, gamma: 1.3};
var BURN_VIS  = {bands: ['SWIR2', 'NIR', 'RED'],  min: 0.02, max: 0.40, gamma: 1.1};
var NBR_VIS   = {min: -0.3, max: 0.7, palette: ['#67001f', '#f7f7f7', '#1a9850']};
var DNBR_VIS  = {min: -0.1, max: 0.8, palette: ['#1a9850', '#ffffbf', '#d73027']};
var P_VIS     = {min: 0, max: 1, palette: ['#ffffff', '#fdae61', '#d73027']};

function draw(y) {
  Map.layers().reset();
  var b = build(y);
  Map.addLayer(b.preRGB,  TRUE_VIS, (y - 1) + ' true colour (before, Mar–May)', false);
  Map.addLayer(b.postRGB, TRUE_VIS, y + ' true colour (after burning, Mar–May)', true);
  Map.addLayer(b.preRGB,  BURN_VIS, (y - 1) + ' SWIR (before)', false);
  Map.addLayer(b.postRGB, BURN_VIS, y + ' SWIR (burn scars = bright red/orange)', false);
  Map.addLayer(b.idx.select('NBR_post'), NBR_VIS, y + ' NBR_post', false);
  Map.addLayer(b.idx.select('dNBR'), DNBR_VIS, y + ' dNBR', false);
  Map.addLayer(b.idx.select('P_jhum'), P_VIS, y + ' P(jhum)', false);
  Map.addLayer(b.jhum.selfMask(), {palette: ['#ff00ff']}, y + ' JHUM MAP', true);
  Map.addLayer(ee.Image().byte().paint(ref.filter(ee.Filter.eq('jhum_year', y)), 1, 2),
               {palette: ['#00ffff']}, y + ' training polygons', false);
  if (CHECKPOINT_ASSET) {
    Map.addLayer(ee.FeatureCollection(CHECKPOINT_ASSET).filter(ee.Filter.eq('year', y)),
                 {color: 'yellow'}, y + ' check points', false);
  }
  info.setValue('Year ' + y + ' — layers loaded. Click the map to inspect a pixel.');
  print('Scenes ' + y + ' (before / after)', b.nScenes[0], b.nScenes[1]);
  current = b;
}

//--------------------------------------------------------------
// 4. PANEL + CLICK INSPECTOR
//--------------------------------------------------------------
var current = null;
var info = ui.Label('Loading…');
var out  = ui.Label('');

var yearSelect = ui.Select({
  items: YEARS.map(String), value: String(YEAR), placeholder: 'year',
  onChange: function (v) { draw(Number(v)); }
});

var panel = ui.Panel({
  widgets: [
    ui.Label('Jhum visual checker (v3 equation)', {fontWeight: 'bold', fontSize: '15px'}),
    ui.Panel([ui.Label('Year:'), yearSelect], ui.Panel.Layout.flow('horizontal')),
    info, out,
    ui.Label('Tip: tick one layer at a time. SWIR view shows burns most clearly.',
             {fontSize: '11px', color: '555555'})
  ],
  style: {position: 'top-left', width: '330px', padding: '8px'}
});
Map.add(panel);

Map.style().set('cursor', 'crosshair');
Map.onClick(function (coords) {
  out.setValue('Reading pixel…');
  var pt = ee.Geometry.Point([coords.lon, coords.lat]);
  var y = Number(yearSelect.getValue());
  var b = current || build(y);
  ee.Image.cat([b.idx, b.jhum.unmask(0)]).reduceRegion({
    reducer: ee.Reducer.first(), geometry: pt, scale: 30
  }).evaluate(function (v) {
    if (!v || v.dNBR === null) { out.setValue('No valid data at this point.'); return; }
    var verdict = v.Jhum === 1 ? 'JHUM' : 'not jhum';
    out.setValue(
      'lon ' + coords.lon.toFixed(5) + ', lat ' + coords.lat.toFixed(5) + '\n' +
      'NBR_pre  ' + v.NBR_pre.toFixed(3) + '\n' +
      'NBR_post ' + v.NBR_post.toFixed(3) + '\n' +
      'dNBR     ' + v.dNBR.toFixed(3) + '\n' +
      'P(jhum)  ' + (v.P_jhum * 100).toFixed(1) + ' %\n' +
      'rule: P >= ' + EQ.cut + ' and dNBR >= ' + EQ.min_dNBR + '\n' +
      '=> ' + verdict);
  });
});

//--------------------------------------------------------------
// 5. START
//--------------------------------------------------------------
Map.centerObject(roi, 8);
draw(YEAR);
