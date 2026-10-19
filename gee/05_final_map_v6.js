//==============================================================
// MIZORAM JHUM — FINAL MAP (v6)
//
// Calibrated on 400 PlanetScope-verified points (2025):
//   user's accuracy 0.87, producer's accuracy 0.86, overall 0.994
//
// RULE (two branches, a pixel is jhum if EITHER is true)
//   A. wide  : NBR_post_wide   <= 0.26  AND  dNBR_wide   >= 0.20
//   B. narrow: NBR_post_narrow <= 0.14  AND  dNBR_narrow >= 0.46
//
//   wide   = 1 Nov – 31 May, compared with the SAME months one year earlier
//            (catches Nov–Feb clearing; same season both sides so deciduous
//             leaf-off cancels out)
//   narrow = 1 Mar – 31 May, same months one year earlier
//            (recovers plots cleared in two consecutive years, which the
//             wide comparison cannot see)
//   Both composites use the 25th percentile of per-image NBR on BOTH sides.
//==============================================================

//--------------------------------------------------------------
// 0. SETTINGS
//--------------------------------------------------------------
var YEARS = [1990, 2000, 2005, 2010, 2015, 2018, 2020, 2022, 2025];
var VERSION   = 'v6';
var MAP_DIR   = 'Jhum_maps_v6';
var ROI_NAME  = 'Mizoram';
var EXPORT_MAPS = true;
var EXPORT_CRS  = 'EPSG:32646';

// Decision rule
var RULE = {
  wide:   {post: 0.26, dnbr: 0.20},
  narrow: {post: 0.14, dnbr: 0.46},
  useNarrowBranch: true          // false = simple single-branch rule
};

// Windows [yearOffset, month, day] and percentiles
var W_PRE_S = [-2, 11, 1], W_PRE_E = [-1, 6, 1];
var W_POST_S = [-1, 11, 1], W_POST_E = [0, 6, 1];
var N_PRE_S = [-1, 3, 1], N_PRE_E = [-1, 6, 1];
var N_POST_S = [0, 3, 1], N_POST_E = [0, 6, 1];
var PCT = 25;                    // same percentile on both sides — essential

//--------------------------------------------------------------
// 1. AREA, MASKS, HARMONIZED LANDSAT 5/7/8/9
//--------------------------------------------------------------
var roi = ee.FeatureCollection('FAO/GAUL/2015/level1')
            .filter(ee.Filter.eq('ADM1_NAME', ROI_NAME));
Map.centerObject(roi, 8);

var WATER_FREE = ee.Image('JRC/GSW1_4/GlobalSurfaceWater')
                   .select('occurrence').unmask(0).lt(50);

function prepL(img, nir, swir2, harmonize) {
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
           .map(function (i) { return prepL(i, 'SR_B4', 'SR_B7', true); })
  .merge(ee.ImageCollection('LANDSAT/LE07/C02/T1_L2').filterBounds(roi)
           .map(function (i) { return prepL(i, 'SR_B4', 'SR_B7', true); }))
  .merge(ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filterBounds(roi)
           .map(function (i) { return prepL(i, 'SR_B5', 'SR_B7', false); }))
  .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2').filterBounds(roi)
           .map(function (i) { return prepL(i, 'SR_B5', 'SR_B7', false); }));

function win(y, w) { return ee.Date.fromYMD(y + w[0], w[1], w[2]); }
function nbrOf(i) { return ee.Image(i).normalizedDifference(['NIR', 'SWIR2']).rename('NBR'); }
function pctImg(col) {
  return ee.Image(ee.Algorithms.If(col.size().gt(0),
           col.map(nbrOf).reduce(ee.Reducer.percentile([PCT])),
           ee.Image.constant(0).updateMask(0)));
}

//--------------------------------------------------------------
// 2. JHUM MAP FOR ONE YEAR
//--------------------------------------------------------------
function jhumFor(y) {
  var wPre  = LS.filterDate(win(y, W_PRE_S),  win(y, W_PRE_E));
  var wPost = LS.filterDate(win(y, W_POST_S), win(y, W_POST_E));
  var nPre  = LS.filterDate(win(y, N_PRE_S),  win(y, N_PRE_E));
  var nPost = LS.filterDate(win(y, N_POST_S), win(y, N_POST_E));

  var preW  = pctImg(wPre).rename('NBR_pre_wide');
  var postW = pctImg(wPost).rename('NBR_post_wide');
  var dW    = preW.subtract(postW).rename('dNBR_wide');
  var preN  = pctImg(nPre).rename('NBR_pre_narrow');
  var postN = pctImg(nPost).rename('NBR_post_narrow');
  var dN    = preN.subtract(postN).rename('dNBR_narrow');

  var branchA = postW.lte(RULE.wide.post).and(dW.gte(RULE.wide.dnbr));
  var branchB = postN.lte(RULE.narrow.post).and(dN.gte(RULE.narrow.dnbr));
  var jhum = RULE.useNarrowBranch ? branchA.or(branchB) : branchA;

  var valid = wPre.select('NIR').count().unmask(0).gt(0)
              .and(wPost.select('NIR').count().unmask(0).gt(0));

  return {
    jhum: jhum.updateMask(valid).updateMask(WATER_FREE).rename('Jhum').clip(roi),
    branchA: branchA.updateMask(valid).clip(roi),
    branchB: branchB.updateMask(valid).clip(roi),
    idx: ee.Image.cat([preW, postW, dW, preN, postN, dN]).updateMask(valid).clip(roi),
    valid: valid.clip(roi),
    nImg: [wPre.size(), wPost.size(), nPost.size()]
  };
}

//--------------------------------------------------------------
// 3. MAP, AREAS, EXPORTS
//    NOTE: printing 9 years of area statistics at once triggers
//    "Too many concurrent aggregations" in the Code Editor.
//    So areas are computed in ONE reduceRegion per year and are
//    printed only if PRINT_TABLE = true (keep it false; use the
//    export task instead, which has no such limit).
//--------------------------------------------------------------
var PRINT_TABLE = false;      // true = print the table (may fail for many years)
var PRINT_YEAR  = 2025;       // area printed for this single year only

var px = ee.Image.pixelArea().divide(1e6);

function areaFeature(y) {
  var m = jhumFor(y);
  var stack = ee.Image.cat([
    px.updateMask(m.jhum.selfMask()).rename('jhum_km2'),
    px.updateMask(m.branchA.selfMask()).rename('jhum_km2_wide_branch'),
    px.updateMask(m.branchB.and(m.branchA.not()).selfMask()).rename('jhum_km2_narrow_only'),
    px.updateMask(m.valid.selfMask()).rename('valid_km2')
  ]);
  var sums = stack.reduceRegion({                 // ONE aggregation per year
    reducer: ee.Reducer.sum(), geometry: roi.geometry(),
    scale: 30, maxPixels: 1e13, tileScale: 4, bestEffort: false
  });
  return ee.Feature(null, sums).set('year', y);
}

var rows = YEARS.map(areaFeature);
var areaTable = ee.FeatureCollection(rows);

YEARS.forEach(function (y) {
  var m = jhumFor(y);
  Map.addLayer(m.jhum.selfMask(), {palette: ['#ff00ff']}, 'Jhum ' + y, y === 2025);
  if (EXPORT_MAPS) {
    var NAME = 'MIZ_Jhum_' + y + '_Landsat30m_' + VERSION;
    Export.image.toDrive({
      image: ee.Image.cat([m.jhum.unmask(0).toByte().rename('Jhum'),
                           m.branchA.unmask(0).toByte().rename('branch_wide'),
                           m.branchB.unmask(0).toByte().rename('branch_narrow')]).clip(roi),
      description: NAME, fileNamePrefix: NAME, folder: MAP_DIR,
      region: roi.geometry(), scale: 30, crs: EXPORT_CRS, maxPixels: 1e13,
      fileFormat: 'GeoTIFF', formatOptions: {cloudOptimized: true}
    });
  }
});

if (PRINT_TABLE) {
  print('JHUM AREA (all years) — may fail with a concurrency error',
        ui.Chart.feature.byFeature(areaTable, 'year',
          ['jhum_km2', 'jhum_km2_wide_branch', 'jhum_km2_narrow_only', 'valid_km2'])
        .setChartType('Table'));
} else {
  print('Jhum area ' + PRINT_YEAR + ' (km2)', areaFeature(PRINT_YEAR));
  print('Full table: run the export task MIZ_Jhum_Area_1990_2025_' + VERSION);
}

Export.table.toDrive({
  collection: areaTable,
  description: 'MIZ_Jhum_Area_1990_2025_' + VERSION,
  fileNamePrefix: 'MIZ_Jhum_Area_1990_2025_' + VERSION,
  folder: MAP_DIR, fileFormat: 'CSV',
  selectors: ['year', 'jhum_km2', 'jhum_km2_wide_branch', 'jhum_km2_narrow_only', 'valid_km2']
});

//--------------------------------------------------------------
// 4. CLICK INSPECTOR (2025 by default)
//--------------------------------------------------------------
var INSPECT_YEAR = 2025;
var mm = jhumFor(INSPECT_YEAR);
var out = ui.Label('Click a pixel to inspect it (' + INSPECT_YEAR + ').');
Map.add(ui.Panel({widgets: [ui.Label('Jhum map v6', {fontWeight: 'bold'}), out],
                  style: {position: 'top-left', width: '320px', padding: '8px'}}));
Map.style().set('cursor', 'crosshair');
Map.onClick(function (c) {
  out.setValue('reading…');
  ee.Image.cat([mm.idx, mm.jhum.unmask(0)]).reduceRegion({
    reducer: ee.Reducer.first(), geometry: ee.Geometry.Point([c.lon, c.lat]), scale: 30
  }).evaluate(function (v) {
    if (!v || v.dNBR_wide === null) { out.setValue('No valid data here.'); return; }
    out.setValue('lon ' + c.lon.toFixed(5) + ', lat ' + c.lat.toFixed(5) +
      '\nWIDE   post ' + v.NBR_post_wide.toFixed(3) + '  dNBR ' + v.dNBR_wide.toFixed(3) +
      '\nNARROW post ' + v.NBR_post_narrow.toFixed(3) + '  dNBR ' + v.dNBR_narrow.toFixed(3) +
      '\nrule A: post<=' + RULE.wide.post + ' & dNBR>=' + RULE.wide.dnbr +
      '\nrule B: post<=' + RULE.narrow.post + ' & dNBR>=' + RULE.narrow.dnbr +
      '\n=> ' + (v.Jhum === 1 ? 'JHUM' : 'not jhum'));
  });
});
