//==============================================================
//  JHUM MAPPER — ready-to-use tool (v1.0.0)
//  Maps current-year shifting cultivation (jhum) clearing from
//  harmonized Landsat 5/7/8/9 with the validated two-branch rule.
//
//  THE RULE  (a pixel is current-year jhum if EITHER branch is true)
//    A  wide  : NBR_post <= 0.26  AND  dNBR >= 0.20
//    B  narrow: NBR_post <= 0.14  AND  dNBR >= 0.46
//
//    NBR  = (NIR - SWIR2) / (NIR + SWIR2),   dNBR = NBR_pre - NBR_post
//    wide pair  : clearing season start -> burn season end  (Mizoram: 1 Nov - 31 May)
//    narrow pair: burn season only                          (Mizoram: 1 Mar - 31 May)
//    "before" = the SAME months one year earlier; 25th percentile of
//    per-scene NBR on BOTH sides of each pair.
//
//  Validated in Mizoram (blind stratified samples):
//    2025 (310 pts): OA 0.995, UA 0.886, PA 0.901, F1 0.893
//    2010 (160 pts): OA 0.996, UA 0.983, PA 0.812, F1 0.890
//  This script reproduces the published map (01_final_map.js) when
//  run with the default settings for Mizoram.
//
//  HOW TO USE
//   1. Set AREA and YEAR in section 0.
//   2. Press Run. Layers appear; the Console prints a report.
//   3. Click a pixel to see its values. Tasks tab: export map + report.
//
//  USING IT ELSEWHERE: the rule has no region-specific term except the
//  calendar. Set CLEAR_START, BURN_START and BURN_END to the local jhum
//  calendar and check a set of known plots before trusting the output.
//  Thresholds are for Landsat surface reflectance harmonized to OLI.
//==============================================================

//--------------------------------------------------------------
// 0. SETTINGS — edit this block only
//--------------------------------------------------------------
// ---- (a) AREA: choose ONE of the three options ----
var AREA_MODE   = 'admin';      // 'admin' | 'asset' | 'draw'
var ADMIN_LEVEL = 1;            // 1 = state/province, 2 = district (FAO GAUL 2015)
var ADMIN_NAME  = 'Mizoram';
var COUNTRY     = 'India';
var AREA_ASSET  = '';           // if AREA_MODE = 'asset': 'projects/YOUR_PROJECT/assets/your_boundary'
// if AREA_MODE = 'draw': draw a polygon with the geometry tools; it is used as "geometry"

// ---- (b) YEAR to map (the year in which the burn season falls) ----
var YEAR = 2025;                // 1986 to the present (Landsat 5 from 1984)

// ---- (c) LOCAL JHUM CALENDAR, [month, day] ----
// Default = Mizoram: slashing from November, burning March-May.
var CLEAR_START = [11, 1];      // start of the clearing season (in the previous year if later than BURN_END)
var BURN_START  = [3, 1];       // start of the burn season
var BURN_END    = [6, 1];       // end of the burn season (exclusive)

// ---- (d) RULE (change only if you recalibrate) ----
var RULE = {
  wide:   {post: 0.26, dnbr: 0.20},
  narrow: {post: 0.14, dnbr: 0.46},
  useNarrowBranch: true
};
var PCT = 25;                   // same percentile on both sides: essential

// ---- (e) OUTPUT ----
var EXPORT_MAP   = true;
var EXPORT_AREA  = true;
var DRIVE_FOLDER = 'Jhum_mapper';
var EXPORT_SCALE = 30;
var EXPORT_CRS   = 'EPSG:4326'; // e.g. 'EPSG:32646' for UTM 46N

//--------------------------------------------------------------
// 1. AREA OF INTEREST
//--------------------------------------------------------------
var aoi;
if (AREA_MODE === 'asset') {
  aoi = ee.FeatureCollection(AREA_ASSET);
} else if (AREA_MODE === 'draw') {
  aoi = ee.FeatureCollection([ee.Feature(geometry)]);
} else {
  var lvl = ADMIN_LEVEL === 2 ? 'FAO/GAUL/2015/level2' : 'FAO/GAUL/2015/level1';
  var nameField = ADMIN_LEVEL === 2 ? 'ADM2_NAME' : 'ADM1_NAME';
  aoi = ee.FeatureCollection(lvl).filter(ee.Filter.eq(nameField, ADMIN_NAME))
          .filter(ee.Filter.eq('ADM0_NAME', COUNTRY));
}
Map.centerObject(aoi, 8);
var AOI_MASK = ee.Image.constant(1).clip(aoi).mask();
var WATER_FREE = ee.Image('JRC/GSW1_4/GlobalSurfaceWater').select('occurrence').unmask(0).lt(50);

//--------------------------------------------------------------
// 2. HARMONIZED LANDSAT 5 / 7 / 8 / 9  (Collection 2, Level-2)
//    TM and ETM+ put on the OLI scale with Roy et al. (2016).
//--------------------------------------------------------------
var NAMES  = ['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2'];
var SLOPES = [0.8474, 0.8483, 0.9047, 0.8462, 0.8937, 0.9071];
var INTER  = [0.0003, 0.0088, 0.0061, 0.0412, 0.0254, 0.0172];
var OLD = ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7'];
var NEW = ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7'];

function prep(img, bands, harmonize) {
  img = ee.Image(img);
  var qa = img.select('QA_PIXEL');
  var clear = qa.bitwiseAnd(1 << 1).eq(0)          // dilated cloud
              .and(qa.bitwiseAnd(1 << 3).eq(0))    // cloud
              .and(qa.bitwiseAnd(1 << 4).eq(0));   // cloud shadow
  var out = img.select(bands, NAMES).multiply(0.0000275).add(-0.2).updateMask(clear);
  if (harmonize) out = out.multiply(SLOPES).add(INTER);
  return ee.Image(out.toFloat().copyProperties(img, ['system:time_start']));
}
var LS = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2').filterBounds(aoi)
           .map(function (i) { return prep(i, OLD, true); })
  .merge(ee.ImageCollection('LANDSAT/LE07/C02/T1_L2').filterBounds(aoi)
           .map(function (i) { return prep(i, OLD, true); }))
  .merge(ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filterBounds(aoi)
           .map(function (i) { return prep(i, NEW, false); }))
  .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2').filterBounds(aoi)
           .map(function (i) { return prep(i, NEW, false); }));

function nbrOf(i) { return ee.Image(i).normalizedDifference(['NIR', 'SWIR2']).rename('NBR'); }
function pctImg(col) {
  return ee.Image(ee.Algorithms.If(col.size().gt(0),
           col.map(nbrOf).reduce(ee.Reducer.percentile([PCT])),
           ee.Image.constant(0).updateMask(0))).rename('NBR').toFloat();
}
function safeCount(col) {
  return ee.Image(ee.Algorithms.If(col.size().gt(0),
           col.select('NIR').count(), ee.Image.constant(0))).rename('n').unmask(0);
}

// windows: the clearing season starts in the previous calendar year when it begins after the burn season ends
var CROSS = (CLEAR_START[0] > BURN_END[0]) ? 1 : 0;
function wideCol(y)   { return LS.filterDate(ee.Date.fromYMD(y - CROSS, CLEAR_START[0], CLEAR_START[1]),
                                             ee.Date.fromYMD(y, BURN_END[0], BURN_END[1])); }
function narrowCol(y) { return LS.filterDate(ee.Date.fromYMD(y, BURN_START[0], BURN_START[1]),
                                             ee.Date.fromYMD(y, BURN_END[0], BURN_END[1])); }

//--------------------------------------------------------------
// 3. THE JHUM MAP
//--------------------------------------------------------------
var wPreCol = wideCol(YEAR - 1), wPostCol = wideCol(YEAR);
var nPreCol = narrowCol(YEAR - 1), nPostCol = narrowCol(YEAR);

var preW  = pctImg(wPreCol).rename('NBR_pre_wide');
var postW = pctImg(wPostCol).rename('NBR_post_wide');
var dW    = preW.subtract(postW).rename('dNBR_wide');
var preN  = pctImg(nPreCol).rename('NBR_pre_narrow');
var postN = pctImg(nPostCol).rename('NBR_post_narrow');
var dN    = preN.subtract(postN).rename('dNBR_narrow');

var A = postW.lte(RULE.wide.post).and(dW.gte(RULE.wide.dnbr));
var B = postN.lte(RULE.narrow.post).and(dN.gte(RULE.narrow.dnbr));
var valid = safeCount(wPreCol).gt(0).and(safeCount(wPostCol).gt(0));

var jhum = (RULE.useNarrowBranch ? A.or(B) : A)
             .updateMask(valid).updateMask(WATER_FREE).updateMask(AOI_MASK).rename('Jhum').clip(aoi);
var indices = ee.Image.cat([preW, postW, dW, preN, postN, dN]).updateMask(valid).clip(aoi);

//--------------------------------------------------------------
// 4. DISPLAY
//--------------------------------------------------------------
var TRUE_VIS = {bands: ['RED', 'GREEN', 'BLUE'], min: 0.02, max: 0.25, gamma: 1.3};
var BURN_VIS = {bands: ['SWIR2', 'NIR', 'RED'],  min: 0.02, max: 0.40, gamma: 1.1};
Map.addLayer(nPreCol.median().clip(aoi),  TRUE_VIS, (YEAR - 1) + ' burn season, true colour', false);
Map.addLayer(nPostCol.median().clip(aoi), TRUE_VIS, YEAR + ' burn season, true colour', true);
Map.addLayer(nPostCol.median().clip(aoi), BURN_VIS, YEAR + ' SWIR view (burns = orange-red)', false);
Map.addLayer(indices.select('dNBR_wide'),
             {min: -0.1, max: 0.8, palette: ['#1a9850', '#ffffbf', '#d73027']}, YEAR + ' dNBR (wide)', false);
Map.addLayer(jhum.selfMask(), {palette: ['#ff00ff']}, YEAR + ' JHUM MAP', true);
Map.addLayer(ee.Image().byte().paint(aoi, 1, 2), {palette: ['#ffffff']}, 'Area boundary');

//--------------------------------------------------------------
// 5. AREA AND DATA-QUALITY REPORT
//--------------------------------------------------------------
var px = ee.Image.pixelArea().divide(1e6);
var sums = ee.Image.cat([
    px.updateMask(jhum.selfMask()).rename('jhum_area_km2'),
    px.updateMask(A.and(jhum).selfMask()).rename('jhum_wide_branch_km2'),
    px.updateMask(B.and(A.not()).and(jhum).selfMask()).rename('jhum_narrow_only_km2'),
    px.updateMask(valid.and(AOI_MASK).selfMask()).rename('valid_area_km2'),
    px.updateMask(AOI_MASK).rename('total_area_km2')])
  .reduceRegion({reducer: ee.Reducer.sum(), geometry: aoi.geometry(), scale: EXPORT_SCALE,
                 maxPixels: 1e13, tileScale: 4});
var report = ee.Feature(null, sums).set({
  year: YEAR,
  images_before_wide: wPreCol.size(), images_after_wide: wPostCol.size(),
  valid_percent: ee.Number(sums.get('valid_area_km2')).divide(sums.get('total_area_km2')).multiply(100)
});
print('REPORT — ' + YEAR, report);
print('Mizoram 2025 with default settings should give about 528 km² (validated map).');
print('If images_before_wide or images_after_wide is small (< 10), treat the result with care.');

//--------------------------------------------------------------
// 6. CLICK INSPECTOR
//--------------------------------------------------------------
var out = ui.Label('Click any pixel to inspect it.');
Map.add(ui.Panel({widgets: [ui.Label('Jhum mapper — ' + YEAR, {fontWeight: 'bold'}), out],
                  style: {position: 'top-left', width: '320px', padding: '8px'}}));
Map.style().set('cursor', 'crosshair');
Map.onClick(function (c) {
  out.setValue('reading…');
  ee.Image.cat([indices, jhum.unmask(0)]).reduceRegion({
    reducer: ee.Reducer.first(), geometry: ee.Geometry.Point([c.lon, c.lat]), scale: 30
  }).evaluate(function (v) {
    if (!v || v.dNBR_wide === null) { out.setValue('No valid data here.'); return; }
    out.setValue('lon ' + c.lon.toFixed(5) + ', lat ' + c.lat.toFixed(5) +
      '\nWIDE    post ' + v.NBR_post_wide.toFixed(3) + '   dNBR ' + v.dNBR_wide.toFixed(3) +
      '\nNARROW  post ' + (v.NBR_post_narrow === null ? 'n/a' : v.NBR_post_narrow.toFixed(3)) +
      '   dNBR ' + (v.dNBR_narrow === null ? 'n/a' : v.dNBR_narrow.toFixed(3)) +
      '\nrule A: post <= ' + RULE.wide.post + ' and dNBR >= ' + RULE.wide.dnbr +
      '\nrule B: post <= ' + RULE.narrow.post + ' and dNBR >= ' + RULE.narrow.dnbr +
      '\n=> ' + (v.Jhum === 1 ? 'JHUM' : 'not jhum'));
  });
});

//--------------------------------------------------------------
// 7. EXPORTS
//--------------------------------------------------------------
var tag = (AREA_MODE === 'admin' ? ADMIN_NAME.replace(/ /g, '_') : 'AOI') + '_' + YEAR;
if (EXPORT_MAP) {
  Export.image.toDrive({
    image: ee.Image.cat([jhum.unmask(0).toByte().rename('Jhum'),
                         A.updateMask(valid).unmask(0).toByte().rename('branch_wide'),
                         B.updateMask(valid).unmask(0).toByte().rename('branch_narrow')]).clip(aoi),
    description: 'JHUM_' + tag + '_30m', fileNamePrefix: 'JHUM_' + tag + '_30m', folder: DRIVE_FOLDER,
    region: aoi.geometry(), scale: EXPORT_SCALE, crs: EXPORT_CRS, maxPixels: 1e13,
    fileFormat: 'GeoTIFF', formatOptions: {cloudOptimized: true}
  });
}
if (EXPORT_AREA) {
  Export.table.toDrive({
    collection: ee.FeatureCollection([report]),
    description: 'JHUM_' + tag + '_report', fileNamePrefix: 'JHUM_' + tag + '_report', folder: DRIVE_FOLDER,
    fileFormat: 'CSV',
    selectors: ['year', 'images_before_wide', 'images_after_wide', 'jhum_area_km2', 'jhum_wide_branch_km2',
                'jhum_narrow_only_km2', 'valid_area_km2', 'total_area_km2', 'valid_percent']
  });
}
