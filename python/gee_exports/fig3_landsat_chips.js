//==============================================================
// FIGURE 3 — LANDSAT CHIPS FOR THE FOUR VISUAL-VALIDATION SITES
//
// For each site and each month shown in Figure 3 (same months as the
// PlanetScope and Sentinel-2 rows), this script builds a Landsat 8 + 9
// monthly composite and exports it as ONE multi-band GeoTIFF per site:
//
//   <MON>_R, <MON>_G, <MON>_B          true colour (reflectance x 10000)
//   <MON>_NIR, <MON>_SWIR1, <MON>_SWIR2 for the SWIR false colour (burns orange-red)
//   + the rule layers at the same site: NBR_post_wide, dNBR_wide,
//     NBR_post_narrow, dNBR_narrow, jhum_2025 (the validated v6 map)
//
// It also exports one CSV with, for every site x month: number of
// scenes, clear observations at the point, whether the gap-fill was
// used, the scene dates, and NBR at the point.
//
// Chips: 3 x 3 km around each point, EPSG:32646, 30 m, snapped to the
// native Landsat grid (no resampling). Crop to any smaller window later.
// Nodata = -9999. Drive folder: Jhum_fig3.
//==============================================================

//--------------------------------------------------------------
// 0. SETTINGS
//--------------------------------------------------------------
var DRIVE_DIR = 'Jhum_fig3';
var HALF_M    = 1500;                 // half-width of the chip (m): 3 x 3 km
var CRS       = 'EPSG:32646';
var GRID      = [30, 0, 15, 0, -30, 15];   // native Landsat C2 UTM grid (pixel corners at 15 m offsets)

// months as [year, month, label]; identical to the Planet / Sentinel-2 rows of Figure 3
var SITES = [
  {panel: 'a', id: 'V2025_051', lon: 93.277870, lat: 23.323005,
   months: [[2024,11,'NOV2024'],[2024,12,'DEC2024'],[2025,1,'JAN2025'],
            [2025,3,'MAR2025'],[2025,5,'MAY2025'],[2025,9,'SEP2025']]},
  {panel: 'b', id: 'V2025_030', lon: 92.909740, lat: 22.668672,
   months: [[2024,11,'NOV2024'],[2024,12,'DEC2024'],[2025,1,'JAN2025'],
            [2025,3,'MAR2025'],[2025,5,'MAY2025'],[2025,8,'AUG2025']]},
  {panel: 'c', id: 'V2025_015', lon: 92.666387, lat: 23.851484,
   months: [[2024,11,'NOV2024'],[2024,12,'DEC2024'],[2025,1,'JAN2025'],
            [2025,3,'MAR2025'],[2025,5,'MAY2025'],[2025,8,'AUG2025']]},
  {panel: 'd', id: 'V2025_100', lon: 92.586077, lat: 22.451729,
   months: [[2024,9,'SEP2024'],[2024,11,'NOV2024'],[2025,1,'JAN2025'],
            [2025,3,'MAR2025'],[2025,5,'MAY2025'],[2025,8,'AUG2025']]}
];

// v6 rule (identical to 00_MIZORAM_FINAL.js)
var RULE = {wide: {post: 0.26, dnbr: 0.20}, narrow: {post: 0.14, dnbr: 0.46}};
var W_PRE_S  = [-2, 11, 1], W_PRE_E  = [-1, 6, 1];
var W_POST_S = [-1, 11, 1], W_POST_E = [ 0, 6, 1];
var N_PRE_S  = [-1,  3, 1], N_PRE_E  = [-1, 6, 1];
var N_POST_S = [ 0,  3, 1], N_POST_E = [ 0, 6, 1];
var PCT = 25, YEAR_MAP = 2025;

//--------------------------------------------------------------
// 1. LANDSAT
//--------------------------------------------------------------
var NAMES = ['B', 'G', 'R', 'NIR', 'SWIR1', 'SWIR2'];
function prep89(img) {               // Landsat 8/9 C2 L2, clear pixels only
  img = ee.Image(img);
  var qa = img.select('QA_PIXEL');
  var clear = qa.bitwiseAnd(1 << 1).eq(0)          // dilated cloud
              .and(qa.bitwiseAnd(1 << 3).eq(0))    // cloud
              .and(qa.bitwiseAnd(1 << 4).eq(0));   // cloud shadow
  return img.select(['SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7'], NAMES)
            .multiply(0.0000275).add(-0.2).updateMask(clear)
            .copyProperties(img, ['system:time_start', 'CLOUD_COVER', 'SPACECRAFT_ID']);
}
function raw89(img) {                // same bands, no cloud mask (gap-fill only)
  img = ee.Image(img);
  return img.select(['SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7'], NAMES)
            .multiply(0.0000275).add(-0.2)
            .copyProperties(img, ['system:time_start', 'CLOUD_COVER']);
}
function l89(geom) {
  return ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
           .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2'))
           .filterBounds(geom);
}

// harmonized L5/7/8/9 for the rule layers (as in the final script)
var SLOPES = [0.8474, 0.8483, 0.9047, 0.8462, 0.8937, 0.9071];
var INTER  = [0.0003, 0.0088, 0.0061, 0.0412, 0.0254, 0.0172];
function prepAny(img, bands, harmonize) {
  img = ee.Image(img);
  var qa = img.select('QA_PIXEL');
  var clear = qa.bitwiseAnd(1 << 1).eq(0).and(qa.bitwiseAnd(1 << 3).eq(0)).and(qa.bitwiseAnd(1 << 4).eq(0));
  var out = img.select(bands, NAMES).multiply(0.0000275).add(-0.2).updateMask(clear);
  if (harmonize) out = out.multiply(SLOPES).add(INTER);
  return ee.Image(out.toFloat().copyProperties(img, ['system:time_start']));
}
var OLD = ['SR_B1','SR_B2','SR_B3','SR_B4','SR_B5','SR_B7'];
var NEW = ['SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7'];
function lsAll(geom) {
  return ee.ImageCollection('LANDSAT/LE07/C02/T1_L2').filterBounds(geom).map(function (i) { return prepAny(i, OLD, true); })
    .merge(ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filterBounds(geom).map(function (i) { return prepAny(i, NEW, false); }))
    .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2').filterBounds(geom).map(function (i) { return prepAny(i, NEW, false); }));
}
function win(y, w) { return ee.Date.fromYMD(y + w[0], w[1], w[2]); }
function pct(col) {
  return ee.Image(ee.Algorithms.If(col.size().gt(0),
    col.map(function (i) { return ee.Image(i).normalizedDifference(['NIR', 'SWIR2']); })
       .reduce(ee.Reducer.percentile([PCT])),
    ee.Image.constant(0).updateMask(0))).toFloat();
}
function ruleLayers(geom) {
  var LS = lsAll(geom), y = YEAR_MAP;
  var pW = pct(LS.filterDate(win(y, W_POST_S), win(y, W_POST_E)));
  var dW = pct(LS.filterDate(win(y, W_PRE_S),  win(y, W_PRE_E))).subtract(pW);
  var pN = pct(LS.filterDate(win(y, N_POST_S), win(y, N_POST_E)));
  var dN = pct(LS.filterDate(win(y, N_PRE_S),  win(y, N_PRE_E))).subtract(pN);
  var jhum = pW.lte(RULE.wide.post).and(dW.gte(RULE.wide.dnbr))
             .or(pN.lte(RULE.narrow.post).and(dN.gte(RULE.narrow.dnbr)));
  return ee.Image.cat([pW.rename('NBR_post_wide'), dW.rename('dNBR_wide'),
                       pN.rename('NBR_post_narrow'), dN.rename('dNBR_narrow'),
                       jhum.rename('jhum_2025').toFloat()]);
}

//--------------------------------------------------------------
// 2. PER SITE: MONTHLY CHIPS + RULE LAYERS -> ONE GEOTIFF
//--------------------------------------------------------------
var rows = [];
SITES.forEach(function (s) {
  var pt  = ee.Geometry.Point([s.lon, s.lat]);
  var box = pt.buffer(HALF_M).bounds();
  var col = l89(box);
  var bands = [];

  s.months.forEach(function (m) {
    var t0 = ee.Date.fromYMD(m[0], m[1], 1), t1 = t0.advance(1, 'month');
    var mon = col.filterDate(t0, t1);
    var clear = mon.map(prep89);
    var comp = clear.median();                                   // clear-sky monthly composite
    // gap-fill: least-cloudy scene of the month, unmasked, only where the composite has no clear pixel
    var best = ee.Image(ee.Algorithms.If(mon.size().gt(0),
                 raw89(mon.sort('CLOUD_COVER').first()), ee.Image.constant([0,0,0,0,0,0]).rename(NAMES).updateMask(0)));
    var filled = ee.Image(ee.Algorithms.If(mon.size().gt(0), comp.unmask(best), best));  // clear median, gaps from least-cloudy scene
    bands.push(filled.multiply(10000).round().toInt16()
                 .rename(NAMES.map(function (b) { return m[2] + '_' + b; })));

    var has = mon.size().gt(0);                    // guard: an empty month has no bands
    var nClear = ee.Image(ee.Algorithms.If(has, clear.select('NIR').count(), ee.Image.constant(0)))
                   .rename('n').unmask(0);
    var nbr = ee.Image(ee.Algorithms.If(has, comp.normalizedDifference(['NIR', 'SWIR2']),
                                        ee.Image.constant(0).updateMask(0))).rename('NBR_comp');
    var atPt = ee.Image.cat([nClear, nbr])
                 .reduceRegion({reducer: ee.Reducer.first(), geometry: pt, scale: 30});
    rows.push(ee.Feature(null, {
      panel: s.panel, point_id: s.id, month: m[2],
      n_scenes: mon.size(),
      scene_dates: ee.List(mon.aggregate_array('system:time_start'))
                     .map(function (d) { return ee.Date(d).format('YYYY-MM-dd'); }).join(';'),
      clear_obs_at_point: atPt.get('n'),
      gapfill_at_point: ee.Number(atPt.get('n')).eq(0),          // 1 = no clear view at the point; scene shown as-is
      NBR_at_point: atPt.get('NBR_comp')
    }));
  });

  var rule = ruleLayers(box).multiply(10000).round().toInt16();   // scaled x 10000; jhum_2025 = 0 / 10000
  var stack = ee.Image.cat(bands.concat([rule])).unmask(-9999).clip(box);

  Export.image.toDrive({
    image: stack, description: 'FIG3_landsat_' + s.panel + '_' + s.id,
    fileNamePrefix: 'FIG3_landsat_' + s.panel + '_' + s.id, folder: DRIVE_DIR,
    region: box, crs: CRS, crsTransform: GRID, maxPixels: 1e9,
    fileFormat: 'GeoTIFF', formatOptions: {noData: -9999}
  });

  // quick look in the Code Editor: May composite, SWIR false colour, with the 30 m analysis cell
  var may = ee.Image(bands[4]);
  Map.addLayer(may.select([5, 3, 2]), {min: 200, max: 4000}, s.panel + ' ' + s.months[4][2] + ' SWIR', s.panel === 'a');
  Map.addLayer(ee.Image().byte().paint(ee.FeatureCollection([ee.Feature(pt.buffer(15).bounds())]), 1, 2),
               {palette: ['yellow']}, s.panel + ' point', s.panel === 'a');
});
Map.centerObject(ee.Geometry.Point([SITES[0].lon, SITES[0].lat]), 14);

var table = ee.FeatureCollection(rows);
print('Scenes and clear observations per site x month', table);
Export.table.toDrive({
  collection: table, description: 'FIG3_landsat_scene_log', fileNamePrefix: 'FIG3_landsat_scene_log',
  folder: DRIVE_DIR, fileFormat: 'CSV',
  selectors: ['panel', 'point_id', 'month', 'n_scenes', 'clear_obs_at_point', 'gapfill_at_point',
              'NBR_at_point', 'scene_dates']
});

print('Run the 5 tasks in the Tasks tab (4 GeoTIFFs + 1 CSV). Folder:', DRIVE_DIR);
print('Band order per month: B, G, R, NIR, SWIR1, SWIR2; then the 5 rule layers. Values x 10000, nodata -9999.');
