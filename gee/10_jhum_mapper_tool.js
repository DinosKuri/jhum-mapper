//==============================================================
//  JHUM MAPPER — ready-to-use tool
//  Maps current-year shifting cultivation (jhum) from harmonized
//  Landsat 5/7/8/9, using the equation calibrated in Mizoram, NE India.
//
//  A pixel is jhum when it is BURNED/BARE this dry season
//  and it CHANGED strongly since the same season one year earlier.
//
//      z = -3.484372 - 11.261907*NBR_post + 9.137461*dNBR
//      P = 1 / (1 + exp(-z))
//      jhum  if  P >= 0.2164  AND  dNBR >= 0.30
//
//  HOW TO USE
//   1. Set your AREA and YEAR in section 0 below.
//   2. Press Run. The map shows true colour before/after and the jhum map.
//   3. Click any pixel to see its values (panel, top-left).
//   4. Open the Tasks tab to export the map (GeoTIFF) and the area table.
//
//  NOTE ON TRANSFER: the equation was calibrated for Mizoram's jhum
//  calendar (slash Nov–Jan, burn Feb–Apr, crop Jun–Sep). In another
//  region, first shift BURN_START / BURN_END to the local burn season,
//  and check a few known plots before trusting the output.
//==============================================================

//--------------------------------------------------------------
// 0. SETTINGS — edit this block only
//--------------------------------------------------------------

// ---- (a) AREA: choose ONE of the three options ----
var AREA_MODE = 'admin';        // 'admin' | 'asset' | 'draw'

var ADMIN_LEVEL = 1;            // 1 = state/province, 2 = district
var ADMIN_NAME  = 'Mizoram';    // e.g. 'Mizoram', 'Nagaland', 'Chin' ...
var COUNTRY     = 'India';      // country the admin unit belongs to

var AREA_ASSET  = '';           // if AREA_MODE = 'asset': 'projects/xxx/assets/my_boundary'
// if AREA_MODE = 'draw': draw a polygon with the geometry tools and it will be used

// ---- (b) YEAR to map ----
var YEAR = 2025;                // any year from 1985 to the present

// ---- (c) BURN SEASON (local jhum calendar) ----
// Default = NE India: burning Mar–May. The "before" image is the SAME
// months one year earlier, which removes seasonal differences.
var BURN_START = [3, 1];        // [month, day] start of burn season
var BURN_END   = [6, 1];        // [month, day] end (exclusive)

// ---- (d) EQUATION (change only if you recalibrate) ----
var EQ = {
  b0: -3.484372,
  b_NBR_post: -11.261907,
  b_dNBR: 9.137461,
  cut: 0.2164,          // probability threshold
  min_dNBR: 0.30        // change guard: rejects towns, riverbeds, landslides
};

// ---- (e) OUTPUT ----
var EXPORT_MAP   = true;              // GeoTIFF of the jhum map
var EXPORT_AREA  = true;              // CSV with jhum area
var DRIVE_FOLDER = 'Jhum_mapper';
var EXPORT_SCALE = 30;                // metres
var EXPORT_CRS   = 'EPSG:4326';       // e.g. 'EPSG:32646' for UTM 46N

//--------------------------------------------------------------
// 1. AREA OF INTEREST
//--------------------------------------------------------------
var aoi;
if (AREA_MODE === 'asset') {
  aoi = ee.FeatureCollection(AREA_ASSET);
} else if (AREA_MODE === 'draw') {
  aoi = ee.FeatureCollection([ee.Feature(geometry)]);   // "geometry" = your drawn polygon
} else {
  var lvl = ADMIN_LEVEL === 2 ? 'FAO/GAUL/2015/level2' : 'FAO/GAUL/2015/level1';
  var nameField = ADMIN_LEVEL === 2 ? 'ADM2_NAME' : 'ADM1_NAME';
  aoi = ee.FeatureCollection(lvl)
          .filter(ee.Filter.eq(nameField, ADMIN_NAME))
          .filter(ee.Filter.eq('ADM0_NAME', COUNTRY));
}
Map.centerObject(aoi, 8);

//--------------------------------------------------------------
// 2. HARMONIZED LANDSAT 5 / 7 / 8 / 9  (Collection 2, Level-2)
//--------------------------------------------------------------
var NAMES  = ['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2'];
var SLOPES = [0.8474, 0.8483, 0.9047, 0.8462, 0.8937, 0.9071];   // Roy et al. 2016
var INTER  = [0.0003, 0.0088, 0.0061, 0.0412, 0.0254, 0.0172];   // TM/ETM+ -> OLI
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

var LS = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2').filterBounds(aoi)
           .map(function (i) { return prep(i, OLD, true); })
  .merge(ee.ImageCollection('LANDSAT/LE07/C02/T1_L2').filterBounds(aoi)
           .map(function (i) { return prep(i, OLD, true); }))
  .merge(ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filterBounds(aoi)
           .map(function (i) { return prep(i, NEW, false); }))
  .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2').filterBounds(aoi)
           .map(function (i) { return prep(i, NEW, false); }));

function season(y) {
  return LS.filterDate(ee.Date.fromYMD(y, BURN_START[0], BURN_START[1]),
                       ee.Date.fromYMD(y, BURN_END[0],   BURN_END[1]));
}
function nbr(img) { return ee.Image(img).normalizedDifference(['NIR', 'SWIR2']).rename('NBR'); }
function pct(col, p) {
  return ee.Image(ee.Algorithms.If(col.size().gt(0),
           col.map(nbr).reduce(ee.Reducer.percentile([p])),
           ee.Image.constant(0).updateMask(0)));
}

//--------------------------------------------------------------
// 3. BUILD THE JHUM MAP
//--------------------------------------------------------------
var beforeCol = season(YEAR - 1);     // same season, previous year
var afterCol  = season(YEAR);         // burn season of the mapping year

var NBR_pre  = pct(beforeCol, 75).rename('NBR_pre');    // typical fallow/vegetated state
var NBR_post = pct(afterCol, 25).rename('NBR_post');    // most burned state
var dNBR     = NBR_pre.subtract(NBR_post).rename('dNBR');

var z = ee.Image.constant(EQ.b0)
          .add(NBR_post.multiply(EQ.b_NBR_post))
          .add(dNBR.multiply(EQ.b_dNBR));
var P = z.multiply(-1).exp().add(1).pow(-1).rename('P_jhum');

var water = ee.Image('JRC/GSW1_4/GlobalSurfaceWater')
              .select('occurrence').unmask(0).lt(50);
var valid = beforeCol.select('NIR').count().unmask(0).gt(0)
            .and(afterCol.select('NIR').count().unmask(0).gt(0));

var jhum = P.gte(EQ.cut).and(dNBR.gte(EQ.min_dNBR))
             .updateMask(valid).updateMask(water).rename('Jhum').clip(aoi);

var indices = ee.Image.cat([NBR_pre, NBR_post, dNBR, P]).updateMask(valid).clip(aoi);

//--------------------------------------------------------------
// 4. DISPLAY
//--------------------------------------------------------------
var TRUE_VIS = {bands: ['RED', 'GREEN', 'BLUE'], min: 0.02, max: 0.25, gamma: 1.3};
var BURN_VIS = {bands: ['SWIR2', 'NIR', 'RED'],  min: 0.02, max: 0.40, gamma: 1.1};

Map.addLayer(beforeCol.median().clip(aoi), TRUE_VIS, (YEAR - 1) + ' true colour (before)', false);
Map.addLayer(afterCol.median().clip(aoi),  TRUE_VIS, YEAR + ' true colour (after burning)', true);
Map.addLayer(beforeCol.median().clip(aoi), BURN_VIS, (YEAR - 1) + ' SWIR view', false);
Map.addLayer(afterCol.median().clip(aoi),  BURN_VIS, YEAR + ' SWIR view (burns = orange-red)', false);
Map.addLayer(indices.select('dNBR'),
             {min: -0.1, max: 0.8, palette: ['#1a9850', '#ffffbf', '#d73027']}, YEAR + ' dNBR', false);
Map.addLayer(indices.select('P_jhum'),
             {min: 0, max: 1, palette: ['#ffffff', '#fdae61', '#d73027']}, YEAR + ' P(jhum)', false);
Map.addLayer(jhum.selfMask(), {palette: ['#ff00ff']}, YEAR + ' JHUM MAP', true);
Map.addLayer(ee.Image().byte().paint(aoi, 1, 2), {palette: ['#ffffff']}, 'Area boundary');

//--------------------------------------------------------------
// 5. AREA AND DATA-QUALITY REPORT
//--------------------------------------------------------------
var px = ee.Image.pixelArea().divide(1e6);
var jhumKm2 = px.updateMask(jhum.selfMask()).reduceRegion({
  reducer: ee.Reducer.sum(), geometry: aoi.geometry(), scale: EXPORT_SCALE,
  maxPixels: 1e13, tileScale: 4}).values().get(0);
var validKm2 = px.updateMask(valid).clip(aoi).reduceRegion({
  reducer: ee.Reducer.sum(), geometry: aoi.geometry(), scale: EXPORT_SCALE,
  maxPixels: 1e13, tileScale: 4}).values().get(0);
var totalKm2 = px.clip(aoi).reduceRegion({
  reducer: ee.Reducer.sum(), geometry: aoi.geometry(), scale: EXPORT_SCALE,
  maxPixels: 1e13, tileScale: 4}).values().get(0);

var report = ee.Feature(null, {
  year: YEAR,
  images_before: beforeCol.size(),
  images_after: afterCol.size(),
  jhum_area_km2: jhumKm2,
  valid_area_km2: validKm2,
  total_area_km2: totalKm2,
  valid_percent: ee.Number(validKm2).divide(totalKm2).multiply(100),
  jhum_percent_of_valid: ee.Number(jhumKm2).divide(validKm2).multiply(100)
});
print('REPORT — ' + YEAR, report);
print('If images_before or images_after is small (<10), results may be unreliable.');

//--------------------------------------------------------------
// 6. CLICK INSPECTOR
//--------------------------------------------------------------
var out = ui.Label('Click any pixel to inspect it.');
Map.add(ui.Panel({
  widgets: [ui.Label('Jhum mapper — ' + YEAR, {fontWeight: 'bold'}), out],
  style: {position: 'top-left', width: '300px', padding: '8px'}
}));
Map.style().set('cursor', 'crosshair');
Map.onClick(function (c) {
  out.setValue('reading…');
  ee.Image.cat([indices, jhum.unmask(0)]).reduceRegion({
    reducer: ee.Reducer.first(), geometry: ee.Geometry.Point([c.lon, c.lat]), scale: 30
  }).evaluate(function (v) {
    if (!v || v.dNBR === null) { out.setValue('No valid data here.'); return; }
    out.setValue('lon ' + c.lon.toFixed(5) + ', lat ' + c.lat.toFixed(5) +
      '\nNBR_pre  ' + v.NBR_pre.toFixed(3) +
      '\nNBR_post ' + v.NBR_post.toFixed(3) +
      '\ndNBR     ' + v.dNBR.toFixed(3) +
      '\nP(jhum)  ' + (v.P_jhum * 100).toFixed(1) + ' %' +
      '\n=> ' + (v.Jhum === 1 ? 'JHUM' : 'not jhum'));
  });
});

//--------------------------------------------------------------
// 7. EXPORTS
//--------------------------------------------------------------
var tag = (AREA_MODE === 'admin' ? ADMIN_NAME.replace(/ /g, '_') : 'AOI') + '_' + YEAR;

if (EXPORT_MAP) {
  Export.image.toDrive({
    image: ee.Image.cat([jhum.unmask(0).toByte(),
                         P.multiply(100).round().toByte().rename('P_jhum_pct')]).clip(aoi),
    description: 'JHUM_' + tag + '_30m',
    fileNamePrefix: 'JHUM_' + tag + '_30m',
    folder: DRIVE_FOLDER,
    region: aoi.geometry(),
    scale: EXPORT_SCALE,
    crs: EXPORT_CRS,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF',
    formatOptions: {cloudOptimized: true}
  });
}
if (EXPORT_AREA) {
  Export.table.toDrive({
    collection: ee.FeatureCollection([report]),
    description: 'JHUM_' + tag + '_report',
    fileNamePrefix: 'JHUM_' + tag + '_report',
    folder: DRIVE_FOLDER,
    fileFormat: 'CSV',
    selectors: ['year', 'images_before', 'images_after', 'jhum_area_km2',
                'valid_area_km2', 'total_area_km2', 'valid_percent', 'jhum_percent_of_valid']
  });
}
