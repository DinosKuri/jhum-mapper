//==============================================================
// FIGURES 1 AND 9: EXPORT LAYERS FOR QGIS   (v2, 7 Oct 2026)
//
// Changes from v1
//   1. NODATA FIX. Earth Engine can write masked pixels as 0. For
//      fallow age, 0 means "cleared this year", so every never-cleared
//      pixel would have painted as freshly cleared. All rasters are now
//      unmasked to an explicit sentinel and tagged as noData.
//   2. PAGE-SCALE PANELS. At 183 mm a 30 m pixel prints at about one
//      fifth of a 300 dpi dot, so neither nearest-neighbour rendering
//      nor a 1-pixel dilation makes single plots visible statewide.
//      The state panels now use 1 km summaries; 30 m data is kept for
//      the detail windows, where pixels are visible.
//   3. ONE COMPUTATION. The 37-year stack is exported once as a
//      multi-band image instead of being recomputed for each layer.
//   4. AREA CHECK. The Console prints the 2025 mapped area, which must
//      read 527.4 km2 (the validated map) before you export anything.
//   5. Detail windows: Landsat 8 + 9 (v1 used 9 only, so more cloud
//      gaps), a SWIR false-colour version where burns stand out, jhum
//      outlines as polygons, and the window frames as a shapefile.
//   6. Hillshade covers a box 15 km beyond the state, so the
//      surrounding terrain can be shown in grey instead of black.
//   7. Neighbouring states and countries exported for labels/context.
//
// Rule and windows are identical to 00_MIZORAM_FINAL.js (v6).//
// NOT THE ROUTE USED FOR THE PUBLISHED FIGURES. Figure 9 and the
// graphical abstract in the paper were drawn in Matplotlib by
// python/fig9/, on a Copernicus GLO-30 hillshade. This QGIS route
// and its SRTM hillshade are kept because they were built first and
// still work; nothing in the paper is rendered from them.
//
// The 2025 mapped area printed here reads 527.4 km2. The paper
// reports 528.2 km2, which is the same map measured on the declared
// EPSG:32646 30 m grid instead of the scene projection. The 0.8 km2
// difference is grid declaration, not a different map.

// All rasters: EPSG:32646. Files land in Drive folder DRIVE_DIR.
//==============================================================

//--------------------------------------------------------------
// 0. SETTINGS
//--------------------------------------------------------------
var YEAR_MAP  = 2025;
var START     = 1989;
var END       = 2025;
var ROI_NAME  = 'Mizoram';
var DRIVE_DIR = 'Jhum_fig9';
var CRS       = 'EPSG:32646';
var AGG_M     = 1000;      // cell size for the state-scale panels (m)
var NODATA    = -1;        // sentinel for every raster except elevation

// Detail windows, centred on validation points a and b of Figure 3
// so the same places appear in both figures. 2.5 km half-width = 5 km.
var INSETS = [
  {name: 'c1', lon: 93.277870, lat: 23.323005, halfWidthM: 2500},
  {name: 'c2', lon: 92.909740, lat: 22.668672, halfWidthM: 2500}
];

var RULE = {wide: {post: 0.26, dnbr: 0.20}, narrow: {post: 0.14, dnbr: 0.46}};
var W_PRE_S  = [-2, 11, 1], W_PRE_E  = [-1, 6, 1];
var W_POST_S = [-1, 11, 1], W_POST_E = [ 0, 6, 1];
var N_PRE_S  = [-1,  3, 1], N_PRE_E  = [-1, 6, 1];
var N_POST_S = [ 0,  3, 1], N_POST_E = [ 0, 6, 1];
var PCT = 25;

//--------------------------------------------------------------
// 1. AREA, MASKS, TERRAIN
//--------------------------------------------------------------
var roi = ee.FeatureCollection('FAO/GAUL/2015/level1')
            .filter(ee.Filter.eq('ADM1_NAME', ROI_NAME));
var ROI_MASK = ee.Image.constant(1).clip(roi).mask();
var WATER_FREE = ee.Image('JRC/GSW1_4/GlobalSurfaceWater')
                   .select('occurrence').unmask(0).lt(50);
var CONTEXT_BOX = roi.geometry().buffer(15000).bounds();   // state + 15 km
Map.centerObject(roi, 8);

var SRTM = ee.Image('USGS/SRTMGL1_003');
var HILLSHADE = ee.Terrain.hillshade(SRTM.multiply(1.6), 315, 45).rename('hillshade');

//--------------------------------------------------------------
// 2. HARMONIZED LANDSAT (identical to the final script)
//--------------------------------------------------------------
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
function nbrOf(i)  { return ee.Image(i).normalizedDifference(['NIR', 'SWIR2']).rename('NBR'); }
function pctImg(col) {
  return ee.Image(ee.Algorithms.If(col.size().gt(0),
           col.map(nbrOf).reduce(ee.Reducer.percentile([PCT])),
           ee.Image.constant(0).updateMask(0))).rename('NBR').toFloat();
}
function safeCount(col) {
  return ee.Image(ee.Algorithms.If(col.size().gt(0),
           col.select('NIR').count(), ee.Image.constant(0))).rename('n').unmask(0);
}
function jhumFor(y) {
  var wPreC  = LS.filterDate(win(y, W_PRE_S),  win(y, W_PRE_E));
  var wPostC = LS.filterDate(win(y, W_POST_S), win(y, W_POST_E));
  var pW = pctImg(wPostC), dW = pctImg(wPreC).subtract(pW);
  var pN = pctImg(LS.filterDate(win(y, N_POST_S), win(y, N_POST_E)));
  var dN = pctImg(LS.filterDate(win(y, N_PRE_S),  win(y, N_PRE_E))).subtract(pN);
  var A = pW.lte(RULE.wide.post).and(dW.gte(RULE.wide.dnbr));
  var B = pN.lte(RULE.narrow.post).and(dN.gte(RULE.narrow.dnbr));
  var valid = safeCount(wPreC).gt(0).and(safeCount(wPostC).gt(0));
  return A.or(B).updateMask(valid).updateMask(WATER_FREE).updateMask(ROI_MASK)
          .rename('jhum').clip(roi);
}

//--------------------------------------------------------------
// 3. ANNUAL STACK AND TEMPORAL METRICS (30 m)
//--------------------------------------------------------------
var YEARS = [];
for (var y = START; y <= END; y++) YEARS.push(y);

var annual = ee.ImageCollection.fromImages(YEARS.map(function (yy) {
  return jhumFor(yy).unmask(0).updateMask(ROI_MASK).toByte().rename('jhum').set('year', yy);
}));

var jhumMap = jhumFor(YEAR_MAP).unmask(0).updateMask(ROI_MASK).toByte().rename('jhum_2025');

var count = annual.sum().rename('jhum_count').toInt16().updateMask(ROI_MASK);

function yearImage(img) {
  return ee.Image.constant(ee.Number(img.get('year'))).toInt16()
           .updateMask(img.eq(1)).rename('yr');
}
var lastClearing = annual.map(yearImage).max();
var fallowAge = ee.Image.constant(END).subtract(lastClearing).toInt16()
                  .rename('fallow_age');                 // masked where never cleared

var annualFull = annual.filter(ee.Filter.gte('year', 1998));   // complete coverage only
var countFull  = annualFull.sum();
var yrFull     = annualFull.map(yearImage);
var cycleMean  = yrFull.max().subtract(yrFull.min())
                   .divide(countFull.subtract(1).max(1))
                   .updateMask(countFull.gte(2)).toFloat().rename('cycle_mean');

//--------------------------------------------------------------
// 4. AREA CHECK  (must read 527.4 km2 before you export)
//--------------------------------------------------------------
var PX_KM2 = ee.Image.pixelArea().divide(1e6);
print('CHECK: mapped jhum ' + YEAR_MAP + ', km2 (validated map = 527.4)',
  PX_KM2.updateMask(jhumMap.eq(1)).reduceRegion({
    reducer: ee.Reducer.sum(), geometry: roi.geometry(),
    scale: 30, maxPixels: 1e13, tileScale: 4}).get('area'));
print('CHECK: elevation range in the state, m (legend ends)',
  SRTM.reduceRegion({reducer: ee.Reducer.minMax(), geometry: roi.geometry(),
    scale: 90, maxPixels: 1e10}));

//--------------------------------------------------------------
// 5. 1 km SUMMARIES FOR THE STATE-SCALE PANELS
//    Input is pinned to 30 m first: composites have no native
//    projection, and reduceResolution needs one.
//--------------------------------------------------------------
function toKm(img) {
  return img.toFloat().reproject({crs: CRS, scale: 30})
            .reduceResolution({reducer: ee.Reducer.mean(), maxPixels: 2048})
            .reproject({crs: CRS, scale: AGG_M});
}
var jhumPct    = toKm(jhumMap).multiply(100).rename('jhum_pct_2025');       // % of cell cleared in 2025
var everPct    = toKm(count.gt(0)).multiply(100).rename('ever_cleared_pct'); // % of cell cleared at least once
var fallowMean = toKm(fallowAge).rename('fallow_age_mean');                  // mean, over cleared pixels only
var countMean  = toKm(count.updateMask(count.gt(0))).rename('clearings_mean'); // mean clearings, cleared pixels

var stack1km = ee.Image.cat([jhumPct, everPct, fallowMean, countMean])
                 .updateMask(ROI_MASK).unmask(NODATA).toFloat().clip(roi);

//--------------------------------------------------------------
// 6. QUICK LOOK
//--------------------------------------------------------------
Map.addLayer(HILLSHADE.clip(CONTEXT_BOX), {min: 0, max: 255}, 'hillshade', false);
Map.addLayer(jhumPct.selfMask().clip(roi), {min: 0, max: 20,
  palette: ['#fde0dd', '#fa9fb5', '#e34a33', '#b30000', '#67000d']}, '2025 jhum, % of 1 km cell');
Map.addLayer(fallowMean.clip(roi), {min: 0, max: 25,
  palette: ['#67000d', '#cb181d', '#fb6a4a', '#fdd0a2', '#c7e9c0', '#41ab5d', '#00441b']},
  'mean years since last clearing (1 km)', false);
Map.addLayer(jhumMap.selfMask().clip(roi), {palette: ['#C1292E']}, 'jhum 2025, 30 m', false);

//--------------------------------------------------------------
// 7. EXPORTS
//--------------------------------------------------------------
function exportImg(img, name, region, scale, noData) {
  Export.image.toDrive({
    image: img, description: name, fileNamePrefix: name, folder: DRIVE_DIR,
    region: region, scale: scale, crs: CRS, maxPixels: 1e13,
    fileFormat: 'GeoTIFF', formatOptions: {cloudOptimized: true, noData: noData}
  });
}

// (a) state panels, 1 km. Bands: jhum_pct_2025, ever_cleared_pct,
//     fallow_age_mean, clearings_mean. -1 = outside the state.
exportImg(stack1km, 'FIG9_state_1km_' + START + '_' + END, roi.geometry(), AGG_M, NODATA);

// (b) 30 m stack, one computation. Bands: jhum_2025 (0/1),
//     fallow_age (years), jhum_count, cycle_mean (years).
//     -1 = outside the state OR never cleared (fallow_age) OR fewer
//     than two clearings since 1998 (cycle_mean).
var stack30 = ee.Image.cat([jhumMap.toFloat(), fallowAge.toFloat(),
                            count.toFloat(), cycleMean])
                .unmask(NODATA).updateMask(ROI_MASK).unmask(NODATA).toFloat().clip(roi);
exportImg(stack30, 'FIG9_stack30m_' + START + '_' + END, roi.geometry(), 30, NODATA);

// (c) terrain over state + 15 km, for Figures 1 and 9
exportImg(HILLSHADE.clip(CONTEXT_BOX).toByte(), 'FIG_hillshade_context', CONTEXT_BOX, 30, 0);
exportImg(SRTM.clip(roi).unmask(-9999).clip(CONTEXT_BOX).toInt16(),
          'FIG_elevation_state', CONTEXT_BOX, 30, -9999);

// (d) boundaries
Export.table.toDrive({collection: roi, description: 'FIG_mizoram_boundary',
  fileNamePrefix: 'FIG_mizoram_boundary', folder: DRIVE_DIR, fileFormat: 'SHP'});
Export.table.toDrive({
  collection: ee.FeatureCollection('FAO/GAUL/2015/level2')
                .filter(ee.Filter.eq('ADM1_NAME', ROI_NAME)),
  description: 'FIG_mizoram_districts_GAUL2015', fileNamePrefix: 'FIG_mizoram_districts_GAUL2015',
  folder: DRIVE_DIR, fileFormat: 'SHP'});
// neighbours, for edge labels (Assam, Manipur, Tripura, Bangladesh, Myanmar)
var NEIGHBOURS = ee.FeatureCollection('FAO/GAUL/2015/level1')
  .filterBounds(roi.geometry().bounds().buffer(40000))
  .filter(ee.Filter.neq('ADM1_NAME', ROI_NAME));
Export.table.toDrive({collection: NEIGHBOURS, description: 'FIG_neighbours',
  fileNamePrefix: 'FIG_neighbours', folder: DRIVE_DIR, fileFormat: 'SHP'});

// (e) detail windows
var frames = [];
INSETS.forEach(function (w) {
  var box = ee.Geometry.Point([w.lon, w.lat]).buffer(w.halfWidthM).bounds();
  frames.push(ee.Feature(box, {name: w.name}));

  // 30 m layers for the window
  exportImg(stack30.clip(box), 'FIG9_' + w.name + '_stack30m', box, 30, NODATA);

  // mapped jhum as polygons, for a clean outline in QGIS
  var polys = jhumMap.selfMask().reproject({crs: CRS, scale: 30}).reduceToVectors({
    geometry: box, crs: CRS, scale: 30, geometryType: 'polygon',
    eightConnected: true, labelProperty: 'jhum', maxPixels: 1e10});
  Export.table.toDrive({collection: polys, description: 'FIG9_' + w.name + '_jhum_polygons',
    fileNamePrefix: 'FIG9_' + w.name + '_jhum_polygons', folder: DRIVE_DIR, fileFormat: 'SHP'});

  // Landsat 8 + 9, 1 Mar - 31 May of the map year (post-burn, pre-monsoon), median
  function rgbPrep(i) {
    var qa = ee.Image(i).select('QA_PIXEL');
    var clear = qa.bitwiseAnd(1 << 1).eq(0).and(qa.bitwiseAnd(1 << 3).eq(0))
                  .and(qa.bitwiseAnd(1 << 4).eq(0));
    return ee.Image(i).select(['SR_B7', 'SR_B5', 'SR_B4', 'SR_B3', 'SR_B2'],
                              ['SWIR2', 'NIR', 'RED', 'GREEN', 'BLUE'])
             .multiply(0.0000275).add(-0.2).updateMask(clear);
  }
  var l89 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
              .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2'))
              .filterBounds(box)
              .filterDate(ee.Date.fromYMD(YEAR_MAP, 3, 1), ee.Date.fromYMD(YEAR_MAP, 6, 1))
              .map(rgbPrep).median().multiply(10000).toInt16();
  exportImg(l89.select(['RED', 'GREEN', 'BLUE']).clip(box).unmask(-9999),
            'FIG9_' + w.name + '_landsat_truecolour', box, 30, -9999);
  exportImg(l89.select(['SWIR2', 'NIR', 'RED']).clip(box).unmask(-9999),
            'FIG9_' + w.name + '_landsat_swir', box, 30, -9999);

  Map.addLayer(ee.Image().byte().paint(ee.FeatureCollection([ee.Feature(box)]), 1, 2),
               {palette: ['#1A1A1A']}, w.name + ' frame');
});
Export.table.toDrive({collection: ee.FeatureCollection(frames),
  description: 'FIG9_detail_frames', fileNamePrefix: 'FIG9_detail_frames',
  folder: DRIVE_DIR, fileFormat: 'SHP'});

print('1. Confirm the area CHECK above reads 527.4 km2.');
print('2. Tasks tab: run every task. The two stack exports are the slow ones.');
print('3. In QGIS, -1 is noData for FIG9_* rasters; set it by hand if not detected.');
