//==============================================================
//  MIZORAM JHUM MAP — FINAL VALIDATED SCRIPT
//  Harmonized Landsat 5 / 7 / 8 / 9, 30 m, 1989 to 2025
//
//  VALIDATED ACCURACY (independent, blind, stratified samples)
//    2025 (OLI/OLI-2, 310 points): OA 0.995, UA 0.886, PA 0.901, F1 0.893
//    2010 (TM/ETM+,   160 points): OA 0.996, UA 0.983, PA 0.812, F1 0.890
//    Error-adjusted area 2025: 518 ± 108 km²;  2010: 454 ± 168 km²
//
//  THE RULE  (a pixel is current-year jhum if EITHER branch is true)
//    A  wide  : NBR_post ≤ 0.26  AND  dNBR ≥ 0.20
//    B  narrow: NBR_post ≤ 0.14  AND  dNBR ≥ 0.46
//
//    NBR  = (NIR − SWIR2) / (NIR + SWIR2)
//    dNBR = NBR_pre − NBR_post
//
//    wide pair   : 1 Nov – 31 May   vs the SAME months one year earlier
//    narrow pair : 1 Mar – 31 May   vs the SAME months one year earlier
//    both composites: 25th percentile of per-scene NBR, on BOTH sides
//
//  WHY THE SAME PERCENTILE ON BOTH SIDES
//    An asymmetric pair (greenest before, darkest after) manufactures
//    change on deciduous forest that never burned. In development this
//    halved user's accuracy (0.44 against 0.89).
//
//  HOW TO USE
//    1. Set YEARS below (one year, or the full series).
//    2. Press Run. Layers appear; the Console prints the area.
//    3. Open the Tasks tab and run the exports.
//
//  Khiangte, L. (2026). jhum-mapper v1.0.0.
//==============================================================

//--------------------------------------------------------------
// 0. SETTINGS  (edit only this block)
//--------------------------------------------------------------
var YEARS = [2025];              // e.g. [2025] or [1990, 2000, 2010, 2020, 2025]
                                 // full series: use the loop below instead
// var YEARS = []; for (var y = 1989; y <= 2025; y++) YEARS.push(y);

var ROI_NAME  = 'Mizoram';       // FAO GAUL level-1 name
var ROI_ASSET = '';              // or your own boundary asset; '' uses GAUL

var EXPORT_MAPS  = true;         // GeoTIFF per year
var EXPORT_AREA  = true;         // CSV with area and data quality
var DRIVE_FOLDER = 'Jhum_final';
var EXPORT_CRS   = 'EPSG:32646'; // UTM 46N

// Decision rule (validated; change only if you recalibrate)
var RULE = {
  wide:   {post: 0.26, dnbr: 0.20},
  narrow: {post: 0.14, dnbr: 0.46},
  useNarrowBranch: true
};

// Composite windows [yearOffset, month, day] and percentile
var W_PRE_S  = [-2, 11, 1], W_PRE_E  = [-1, 6, 1];
var W_POST_S = [-1, 11, 1], W_POST_E = [ 0, 6, 1];
var N_PRE_S  = [-1,  3, 1], N_PRE_E  = [-1, 6, 1];
var N_POST_S = [ 0,  3, 1], N_POST_E = [ 0, 6, 1];
var PCT = 25;

//--------------------------------------------------------------
// 1. AREA AND MASKS
//--------------------------------------------------------------
var roi = ROI_ASSET ? ee.FeatureCollection(ROI_ASSET)
        : ee.FeatureCollection('FAO/GAUL/2015/level1')
            .filter(ee.Filter.eq('ADM1_NAME', ROI_NAME));
var ROI_MASK = ee.Image.constant(1).clip(roi).mask();
Map.centerObject(roi, 8);

var WATER_FREE = ee.Image('JRC/GSW1_4/GlobalSurfaceWater')
                   .select('occurrence').unmask(0).lt(50);

//--------------------------------------------------------------
// 2. HARMONIZED LANDSAT 5 / 7 / 8 / 9  (Collection 2, Level-2)
//    TM and ETM+ are put on the OLI scale with Roy et al. (2016).
//--------------------------------------------------------------
var NAMES  = ['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2'];
var SLOPES = [0.8474, 0.8483, 0.9047, 0.8462, 0.8937, 0.9071];
var INTER  = [0.0003, 0.0088, 0.0061, 0.0412, 0.0254, 0.0172];
var OLD = ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7'];
var NEW = ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7'];

function prep(img, bands, harmonize) {
  img = ee.Image(img);
  var qa = img.select('QA_PIXEL');
  var clear = qa.bitwiseAnd(1 << 1).eq(0)        // dilated cloud
              .and(qa.bitwiseAnd(1 << 3).eq(0))  // cloud
              .and(qa.bitwiseAnd(1 << 4).eq(0)); // cloud shadow
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
function nbrOf(i)  { return ee.Image(i).normalizedDifference(['NIR', 'SWIR2']).rename('NBR'); }

// 25th percentile of per-scene NBR; an empty window returns a masked image
function pctImg(col) {
  return ee.Image(ee.Algorithms.If(col.size().gt(0),
           col.map(nbrOf).reduce(ee.Reducer.percentile([PCT])),
           ee.Image.constant(0).updateMask(0))).rename('NBR').toFloat();
}
// an empty collection has no bands, which breaks unmask(): guard it
function safeCount(col) {
  return ee.Image(ee.Algorithms.If(col.size().gt(0),
           col.select('NIR').count(), ee.Image.constant(0))).rename('n').unmask(0);
}

//--------------------------------------------------------------
// 3. THE JHUM MAP FOR ONE YEAR
//--------------------------------------------------------------
function jhumFor(y) {
  var wPreCol  = LS.filterDate(win(y, W_PRE_S),  win(y, W_PRE_E));
  var wPostCol = LS.filterDate(win(y, W_POST_S), win(y, W_POST_E));
  var nPreCol  = LS.filterDate(win(y, N_PRE_S),  win(y, N_PRE_E));
  var nPostCol = LS.filterDate(win(y, N_POST_S), win(y, N_POST_E));

  var preW  = pctImg(wPreCol).rename('NBR_pre_wide');
  var postW = pctImg(wPostCol).rename('NBR_post_wide');
  var dW    = preW.subtract(postW).rename('dNBR_wide');
  var preN  = pctImg(nPreCol).rename('NBR_pre_narrow');
  var postN = pctImg(nPostCol).rename('NBR_post_narrow');
  var dN    = preN.subtract(postN).rename('dNBR_narrow');

  var A = postW.lte(RULE.wide.post).and(dW.gte(RULE.wide.dnbr));
  var B = postN.lte(RULE.narrow.post).and(dN.gte(RULE.narrow.dnbr));
  var jhum = RULE.useNarrowBranch ? A.or(B) : A;

  var valid = safeCount(wPreCol).gt(0).and(safeCount(wPostCol).gt(0));

  return {
    jhum:    jhum.updateMask(valid).updateMask(WATER_FREE).updateMask(ROI_MASK)
                 .rename('Jhum').clip(roi),
    branchA: A.updateMask(valid).updateMask(ROI_MASK).clip(roi),
    branchB: B.updateMask(valid).updateMask(ROI_MASK).clip(roi),
    idx:     ee.Image.cat([preW, postW, dW, preN, postN, dN])
               .updateMask(valid).updateMask(ROI_MASK).clip(roi),
    valid:   valid.updateMask(ROI_MASK).clip(roi),
    beforeRGB: wPreCol.median().clip(roi),
    afterRGB:  wPostCol.median().clip(roi),
    nScenes: [wPreCol.size(), wPostCol.size()]
  };
}

//--------------------------------------------------------------
// 4. DISPLAY, AREA AND EXPORTS
//--------------------------------------------------------------
var TRUE_VIS = {bands: ['RED', 'GREEN', 'BLUE'], min: 0.02, max: 0.25, gamma: 1.3};
var BURN_VIS = {bands: ['SWIR2', 'NIR', 'RED'],  min: 0.02, max: 0.40, gamma: 1.1};
var px = ee.Image.pixelArea().divide(1e6);
var rows = [];
var last = null;

YEARS.forEach(function (y) {
  var m = jhumFor(y);
  last = m;

  Map.addLayer(m.beforeRGB, TRUE_VIS, (y - 1) + ' true colour (before)', false);
  Map.addLayer(m.afterRGB,  BURN_VIS, y + ' SWIR (burns are orange-red)', false);
  Map.addLayer(m.idx.select('dNBR_wide'),
    {min: -0.1, max: 0.8, palette: ['#1a9850', '#ffffbf', '#d73027']}, y + ' dNBR', false);
  Map.addLayer(m.jhum.selfMask(), {palette: ['#ff00ff']}, y + ' JHUM MAP',
               y === YEARS[YEARS.length - 1]);

  // one aggregation per year keeps the Console within its limits
  var stack = ee.Image.cat([
    px.updateMask(m.jhum.selfMask()).rename('jhum_km2'),
    px.updateMask(m.branchA.selfMask()).rename('jhum_km2_wide_branch'),
    px.updateMask(m.branchB.and(m.branchA.not()).selfMask()).rename('jhum_km2_narrow_only'),
    px.updateMask(m.valid.selfMask()).rename('valid_km2')]);
  rows.push(ee.Feature(null, stack.reduceRegion({
    reducer: ee.Reducer.sum(), geometry: roi.geometry(),
    scale: 30, maxPixels: 1e13, tileScale: 4})).set('year', y));

  if (EXPORT_MAPS) {
    var NAME = 'MIZ_Jhum_' + y + '_Landsat30m_final';
    Export.image.toDrive({
      image: ee.Image.cat([m.jhum.unmask(0).toByte().rename('Jhum'),
                           m.branchA.unmask(0).toByte().rename('branch_wide'),
                           m.branchB.unmask(0).toByte().rename('branch_narrow')]).clip(roi),
      description: NAME, fileNamePrefix: NAME, folder: DRIVE_FOLDER,
      region: roi.geometry(), scale: 30, crs: EXPORT_CRS, maxPixels: 1e13,
      fileFormat: 'GeoTIFF', formatOptions: {cloudOptimized: true}
    });
  }
});

var areaTable = ee.FeatureCollection(rows);
print('Jhum area for ' + YEARS[YEARS.length - 1] + ' (km²)', rows[rows.length - 1]);
print('Scenes used, before / after', last.nScenes[0], last.nScenes[1]);
print('If more than about five years are listed, read the full table from the export task.');

if (EXPORT_AREA) {
  Export.table.toDrive({
    collection: areaTable,
    description: 'MIZ_Jhum_Area_final_' + YEARS[0] + '_' + YEARS[YEARS.length - 1],
    fileNamePrefix: 'MIZ_Jhum_Area_final_' + YEARS[0] + '_' + YEARS[YEARS.length - 1],
    folder: DRIVE_FOLDER, fileFormat: 'CSV',
    selectors: ['year', 'jhum_km2', 'jhum_km2_wide_branch', 'jhum_km2_narrow_only', 'valid_km2']
  });
}

//--------------------------------------------------------------
// 5. CLICK INSPECTOR (last year in the list)
//--------------------------------------------------------------
var out = ui.Label('Click any pixel to see its values.');
Map.add(ui.Panel({
  widgets: [ui.Label('Mizoram jhum map, ' + YEARS[YEARS.length - 1], {fontWeight: 'bold'}), out],
  style: {position: 'top-left', width: '320px', padding: '8px'}}));
Map.style().set('cursor', 'crosshair');
Map.onClick(function (c) {
  out.setValue('reading…');
  ee.Image.cat([last.idx, last.jhum.unmask(0)]).reduceRegion({
    reducer: ee.Reducer.first(), geometry: ee.Geometry.Point([c.lon, c.lat]), scale: 30
  }).evaluate(function (v) {
    if (!v || v.dNBR_wide === null) { out.setValue('No valid data here.'); return; }
    out.setValue('lon ' + c.lon.toFixed(5) + ', lat ' + c.lat.toFixed(5) +
      '\nWIDE    post ' + v.NBR_post_wide.toFixed(3)   + '   dNBR ' + v.dNBR_wide.toFixed(3) +
      '\nNARROW  post ' + v.NBR_post_narrow.toFixed(3) + '   dNBR ' + v.dNBR_narrow.toFixed(3) +
      '\nrule A: post ≤ ' + RULE.wide.post   + ' and dNBR ≥ ' + RULE.wide.dnbr +
      '\nrule B: post ≤ ' + RULE.narrow.post + ' and dNBR ≥ ' + RULE.narrow.dnbr +
      '\n=> ' + (v.Jhum === 1 ? 'JHUM' : 'not jhum'));
  });
});
