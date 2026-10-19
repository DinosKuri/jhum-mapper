//==============================================================
// FIGURE 9 — EXPORT LAYERS FOR QGIS
//
// Produces every raster Figure 9 needs, ready to style in QGIS:
//   1. jhum_2025            1 = mapped jhum, 0 = not, nodata outside
//   2. fallow_age           years since last clearing (1989-2025)
//   3. jhum_count           number of clearings per pixel
//   4. cycle_mean           mean years between clearings (1998-2025)
//   5. hillshade            SRTM hillshade for the base
//   6. elevation            SRTM, for an optional relief tint
//   7. two detail windows   the same layers clipped to 1:50,000 insets
//
// Everything is exported in EPSG:32646 at 30 m so the layers align
// in QGIS without reprojection.
//
// NOT THE ROUTE USED FOR THE PUBLISHED FIGURE 9. The published
// figure was drawn in Matplotlib by python/fig9/, on a Copernicus
// GLO-30 hillshade. The SRTM hillshade and elevation layers here
// belong to this QGIS route only and are not in the paper.
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

// Detail windows for the insets. Replace with areas you know well;
// these are centred on validation points a and b of Figure 3.
var INSETS = [
  {name: 'inset_A', lon: 93.277870, lat: 23.323005, halfWidthM: 2500},
  {name: 'inset_B', lon: 92.909740, lat: 22.668672, halfWidthM: 2500}
];

var RULE = {wide: {post: 0.26, dnbr: 0.20}, narrow: {post: 0.14, dnbr: 0.46}};
var W_PRE_S = [-2, 11, 1], W_PRE_E = [-1, 6, 1];
var W_POST_S = [-1, 11, 1], W_POST_E = [0, 6, 1];
var N_PRE_S = [-1, 3, 1], N_PRE_E = [-1, 6, 1];
var N_POST_S = [0, 3, 1], N_POST_E = [0, 6, 1];
var PCT = 25;

//--------------------------------------------------------------
// 1. AREA, TERRAIN, LANDSAT
//--------------------------------------------------------------
var roi = ee.FeatureCollection('FAO/GAUL/2015/level1')
            .filter(ee.Filter.eq('ADM1_NAME', ROI_NAME));
var ROI_MASK = ee.Image.constant(1).clip(roi).mask();
var WATER_FREE = ee.Image('JRC/GSW1_4/GlobalSurfaceWater')
                   .select('occurrence').unmask(0).lt(50);
Map.centerObject(roi, 8);

var SRTM = ee.Image('USGS/SRTMGL1_003');
var HILLSHADE = ee.Terrain.hillshade(SRTM.multiply(1.6), 315, 45).rename('hillshade');

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
var LS = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2').filterBounds(roi)
           .map(function (i) { return prep(i, 'SR_B4', 'SR_B7', true); })
  .merge(ee.ImageCollection('LANDSAT/LE07/C02/T1_L2').filterBounds(roi)
           .map(function (i) { return prep(i, 'SR_B4', 'SR_B7', true); }))
  .merge(ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filterBounds(roi)
           .map(function (i) { return prep(i, 'SR_B5', 'SR_B7', false); }))
  .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2').filterBounds(roi)
           .map(function (i) { return prep(i, 'SR_B5', 'SR_B7', false); }));

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
function jhumFor(y) {
  var wPreC = LS.filterDate(win(y, W_PRE_S), win(y, W_PRE_E));
  var wPostC = LS.filterDate(win(y, W_POST_S), win(y, W_POST_E));
  var pW = pctImg(wPostC), dW = pctImg(wPreC).subtract(pW);
  var pN = pctImg(LS.filterDate(win(y, N_POST_S), win(y, N_POST_E)));
  var dN = pctImg(LS.filterDate(win(y, N_PRE_S), win(y, N_PRE_E))).subtract(pN);
  var A = pW.lte(RULE.wide.post).and(dW.gte(RULE.wide.dnbr));
  var B = pN.lte(RULE.narrow.post).and(dN.gte(RULE.narrow.dnbr));
  var valid = safeCount(wPreC).gt(0).and(safeCount(wPostC).gt(0));
  return A.or(B).updateMask(valid).updateMask(WATER_FREE).updateMask(ROI_MASK)
          .rename('jhum').clip(roi);
}

//--------------------------------------------------------------
// 2. ANNUAL STACK AND TEMPORAL METRICS
//--------------------------------------------------------------
var YEARS = [];
for (var y = START; y <= END; y++) YEARS.push(y);

var annual = ee.ImageCollection.fromImages(YEARS.map(function (yy) {
  return jhumFor(yy).unmask(0).updateMask(ROI_MASK).rename('jhum').set('year', yy);
}));

var count = annual.sum().rename('jhum_count').toInt16().updateMask(ROI_MASK).clip(roi);

var yearImgs = annual.map(function (img) {
  return ee.Image.constant(ee.Number(img.get('year'))).toInt16()
           .updateMask(img.eq(1)).rename('yr');
});
var lastClearing = yearImgs.max();
var fallowAge = ee.Image.constant(END).subtract(lastClearing)
                  .rename('fallow_age').toInt16().updateMask(ROI_MASK).clip(roi);

// cycle metrics use 1998 onward only (complete observation)
var annualFull = annual.filter(ee.Filter.gte('year', 1998));
var countFull = annualFull.sum();
var yrFull = annualFull.map(function (img) {
  return ee.Image.constant(ee.Number(img.get('year'))).toInt16()
           .updateMask(img.eq(1)).rename('yr'); });
var cycleMean = yrFull.max().subtract(yrFull.min())
                  .divide(countFull.subtract(1).max(1))
                  .updateMask(countFull.gte(2))
                  .rename('cycle_mean').toFloat().updateMask(ROI_MASK).clip(roi);

var jhum2025 = jhumFor(YEAR_MAP).unmask(0).toByte().updateMask(ROI_MASK).clip(roi);

//--------------------------------------------------------------
// 3. DISPLAY (for a quick look before exporting)
//--------------------------------------------------------------
Map.addLayer(HILLSHADE.clip(roi), {min: 0, max: 255}, 'hillshade', false);
Map.addLayer(jhum2025.selfMask(), {palette: ['#C1292E']}, 'jhum ' + YEAR_MAP);
Map.addLayer(fallowAge, {min: 0, max: 25,
  palette: ['#67000d', '#cb181d', '#fb6a4a', '#fdd0a2', '#c7e9c0', '#41ab5d', '#00441b']},
  'fallow age', false);
Map.addLayer(count, {min: 1, max: 6,
  palette: ['#fee5d9', '#fcae91', '#fb6a4a', '#de2d26', '#a50f15', '#67000d']},
  'number of clearings', false);

//--------------------------------------------------------------
// 4. EXPORTS, STATEWIDE
//--------------------------------------------------------------
function exportImg(img, name, dtype) {
  Export.image.toDrive({
    image: dtype === 'byte' ? img.toByte() : img,
    description: name, fileNamePrefix: name, folder: DRIVE_DIR,
    region: roi.geometry(), scale: 30, crs: CRS, maxPixels: 1e13,
    fileFormat: 'GeoTIFF', formatOptions: {cloudOptimized: true}
  });
}
exportImg(jhum2025, 'FIG9_jhum_' + YEAR_MAP, 'byte');
exportImg(fallowAge, 'FIG9_fallow_age_' + START + '_' + END);
exportImg(count, 'FIG9_jhum_count_' + START + '_' + END);
exportImg(cycleMean, 'FIG9_cycle_mean_1998_' + END);
exportImg(HILLSHADE.clip(roi), 'FIG9_hillshade', 'byte');
exportImg(SRTM.updateMask(ROI_MASK).clip(roi).toInt16(), 'FIG9_elevation');

// state boundary, for the frame and the locator
Export.table.toDrive({
  collection: roi, description: 'FIG9_mizoram_boundary',
  fileNamePrefix: 'FIG9_mizoram_boundary', folder: DRIVE_DIR, fileFormat: 'SHP'});

// district boundaries, for labels
var districts = ee.FeatureCollection('FAO/GAUL/2015/level2')
                  .filter(ee.Filter.eq('ADM1_NAME', ROI_NAME));
Export.table.toDrive({
  collection: districts, description: 'FIG9_mizoram_districts',
  fileNamePrefix: 'FIG9_mizoram_districts', folder: DRIVE_DIR, fileFormat: 'SHP'});

//--------------------------------------------------------------
// 5. DETAIL WINDOWS FOR THE INSETS
//--------------------------------------------------------------
INSETS.forEach(function (w) {
  var box = ee.Geometry.Point([w.lon, w.lat]).buffer(w.halfWidthM).bounds();
  Map.addLayer(ee.Image().byte().paint(ee.FeatureCollection([ee.Feature(box)]), 1, 3),
               {palette: ['#1A1A1A']}, w.name + ' frame', false);

  Export.image.toDrive({
    image: jhum2025.toByte().clip(box), description: 'FIG9_' + w.name + '_jhum',
    fileNamePrefix: 'FIG9_' + w.name + '_jhum', folder: DRIVE_DIR,
    region: box, scale: 30, crs: CRS, maxPixels: 1e10, fileFormat: 'GeoTIFF'});

  Export.image.toDrive({
    image: fallowAge.clip(box), description: 'FIG9_' + w.name + '_fallow_age',
    fileNamePrefix: 'FIG9_' + w.name + '_fallow_age', folder: DRIVE_DIR,
    region: box, scale: 30, crs: CRS, maxPixels: 1e10, fileFormat: 'GeoTIFF'});

  // true-colour Landsat for the inset background, after-season composite
  var rgb = LS.filterDate(win(YEAR_MAP, W_POST_S), win(YEAR_MAP, W_POST_E));
  Export.image.toDrive({
    image: ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')
             .filterBounds(box)
             .filterDate(ee.Date.fromYMD(YEAR_MAP, 2, 1), ee.Date.fromYMD(YEAR_MAP, 5, 1))
             .map(function (i) {
               var qa = ee.Image(i).select('QA_PIXEL');
               var clear = qa.bitwiseAnd(1 << 3).eq(0).and(qa.bitwiseAnd(1 << 4).eq(0));
               return ee.Image(i).select(['SR_B4', 'SR_B3', 'SR_B2'], ['R', 'G', 'B'])
                        .multiply(0.0000275).add(-0.2).updateMask(clear); })
             .median().multiply(10000).toInt16().clip(box),
    description: 'FIG9_' + w.name + '_landsat_rgb',
    fileNamePrefix: 'FIG9_' + w.name + '_landsat_rgb', folder: DRIVE_DIR,
    region: box, scale: 30, crs: CRS, maxPixels: 1e10, fileFormat: 'GeoTIFF'});
});

print('Run every task in the Tasks tab. Files land in Drive folder:', DRIVE_DIR);
print('All rasters are EPSG:32646 at 30 m, so they align in QGIS with no reprojection.');
