//==============================================================
// MIZORAM JHUM — PART E2: FINAL VALIDATION SAMPLE (v6 map)
//
// Draws an INDEPENDENT, BLIND check-point sample from the v6 map.
// Three strata (Olofsson et al. 2014 good practice):
//   1 = mapped jhum
//   2 = NOT mapped jhum but within 300 m of mapped jhum  (where
//       omissions concentrate — sampling here sharpens producer's accuracy)
//   3 = NOT mapped jhum and far from it
//
// RULES FOR THIS ROUND (this is what makes it valid):
//   - label each point ONCE, from imagery only
//   - never look at the map class or the index values while labelling
//   - never revisit a point after learning whether it agreed
//   - thresholds are already locked; nothing is tuned after this
//==============================================================

//--------------------------------------------------------------
// 0. SETTINGS
//--------------------------------------------------------------
var YEAR      = 2025;
var VERSION   = 'v6';
var DRIVE_DIR = 'Jhum_validation_v6';
var ROI_NAME  = 'Mizoram';
var SEED      = 777;

// points per stratum: [mapped jhum, near-jhum, far]
var N_POINTS = {jhum: 70, near: 70, far: 60};
var NEAR_M   = 300;        // "near" = within 300 m of mapped jhum

// exclude the calibration points so the two samples stay independent
var CALIB_ASSET  = 'projects/YOUR_PROJECT/assets/Label_2025';  // '' to skip
var CALIB_BUFFER = 100;    // metres

// v6 rule (locked — do not change)
var RULE = {wide: {post: 0.26, dnbr: 0.20}, narrow: {post: 0.14, dnbr: 0.46}, useNarrowBranch: true};
var W_PRE_S = [-2, 11, 1], W_PRE_E = [-1, 6, 1];
var W_POST_S = [-1, 11, 1], W_POST_E = [0, 6, 1];
var N_PRE_S = [-1, 3, 1], N_PRE_E = [-1, 6, 1];
var N_POST_S = [0, 3, 1], N_POST_E = [0, 6, 1];
var PCT = 25;

//--------------------------------------------------------------
// 1. LANDSAT AND THE v6 MAP (identical to the final map script)
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
function jhumFor(y) {
  var wPre = LS.filterDate(win(y, W_PRE_S), win(y, W_PRE_E));
  var wPost = LS.filterDate(win(y, W_POST_S), win(y, W_POST_E));
  var nPre = LS.filterDate(win(y, N_PRE_S), win(y, N_PRE_E));
  var nPost = LS.filterDate(win(y, N_POST_S), win(y, N_POST_E));
  var postW = pctImg(wPost), dW = pctImg(wPre).subtract(postW);
  var postN = pctImg(nPost), dN = pctImg(nPre).subtract(postN);
  var A = postW.lte(RULE.wide.post).and(dW.gte(RULE.wide.dnbr));
  var B = postN.lte(RULE.narrow.post).and(dN.gte(RULE.narrow.dnbr));
  var jhum = RULE.useNarrowBranch ? A.or(B) : A;
  var valid = wPre.select('NIR').count().unmask(0).gt(0)
              .and(wPost.select('NIR').count().unmask(0).gt(0));
  return ee.Image.cat([
      jhum.rename('map_jhum').toByte(),
      A.rename('branch_wide').toByte(), B.rename('branch_narrow').toByte(),
      postW.rename('NBR_post_wide'), dW.rename('dNBR_wide'),
      postN.rename('NBR_post_narrow'), dN.rename('dNBR_narrow')])
    .updateMask(valid).updateMask(WATER_FREE).clip(roi);
}

var m = jhumFor(YEAR);

//--------------------------------------------------------------
// 2. THREE STRATA
//--------------------------------------------------------------
var jhumMask = m.select('map_jhum').eq(1);
var near = jhumMask.unmask(0)
             .focalMax({radius: NEAR_M, units: 'meters'})
             .and(jhumMask.not());

var strata = ee.Image(3).where(near, 2).where(jhumMask, 1)
               .rename('stratum').toByte()
               .updateMask(m.select('map_jhum').mask());

if (CALIB_ASSET) {
  var calib = ee.FeatureCollection(CALIB_ASSET).map(function (f) {
    var x = ee.Number.parse(ee.Algorithms.String(f.get('lon')));
    var y = ee.Number.parse(ee.Algorithms.String(f.get('lat')));
    return ee.Feature(ee.Geometry.Point([x, y]).buffer(CALIB_BUFFER));
  });
  strata = strata.updateMask(ee.Image.constant(1).paint(calib, 0));
}

//--------------------------------------------------------------
// 3. SAMPLE
//--------------------------------------------------------------
var pts = strata.addBands(m).stratifiedSample({
  numPoints: 0, classBand: 'stratum', region: roi.geometry(), scale: 30,
  classValues: [1, 2, 3],
  classPoints: [N_POINTS.jhum, N_POINTS.near, N_POINTS.far],
  seed: SEED, tileScale: 8, geometries: true, dropNulls: true
}).randomColumn('shuffle', SEED).sort('shuffle');

var list = pts.toList(pts.size());
pts = ee.FeatureCollection(ee.List.sequence(0, list.size().subtract(1)).map(function (i) {
  var f = ee.Feature(list.get(i));
  var c = f.geometry().coordinates();
  return f.set({
    point_id: ee.String('V' + YEAR + '_').cat(ee.Number(i).add(1).format('%03d')),
    year: YEAR,
    lon: ee.Number(c.get(0)).format('%.6f'),
    lat: ee.Number(c.get(1)).format('%.6f'),
    ref_jhum: '', ref_confidence: '', ref_landcover: '', notes: ''
  });
}));

var tag = YEAR + '_' + VERSION;

// BLIND: no stratum, no map class, no index values
Export.table.toDrive({collection: pts,
  description: 'MIZ_Jhum_V1_BLIND_' + tag, fileNamePrefix: 'MIZ_Jhum_V1_BLIND_' + tag,
  folder: DRIVE_DIR, fileFormat: 'CSV',
  selectors: ['point_id', 'year', 'lon', 'lat', 'ref_jhum', 'ref_confidence', 'ref_landcover', 'notes']});
Export.table.toDrive({collection: pts.select(['point_id']),
  description: 'MIZ_Jhum_V1_BLIND_KML_' + tag, fileNamePrefix: 'MIZ_Jhum_V1_BLIND_' + tag,
  folder: DRIVE_DIR, fileFormat: 'KML'});
Export.table.toDrive({
  collection: pts.map(function (f) {
    return ee.Feature(f.geometry().buffer(15).bounds(), {point_id: f.get('point_id')}); }),
  description: 'MIZ_Jhum_V1_PixelSquares_KML_' + tag,
  fileNamePrefix: 'MIZ_Jhum_V1_PixelSquares_' + tag,
  folder: DRIVE_DIR, fileFormat: 'KML'});

// KEY: open only AFTER all labelling is finished
Export.table.toDrive({collection: pts,
  description: 'MIZ_Jhum_V2_KEY_' + tag, fileNamePrefix: 'MIZ_Jhum_V2_KEY_' + tag,
  folder: DRIVE_DIR, fileFormat: 'CSV',
  selectors: ['point_id', 'stratum', 'map_jhum', 'branch_wide', 'branch_narrow',
              'NBR_post_wide', 'dNBR_wide', 'NBR_post_narrow', 'dNBR_narrow']});

//--------------------------------------------------------------
// 4. STRATUM AREAS (needed for the weighted accuracy)
//--------------------------------------------------------------
var px = ee.Image.pixelArea().divide(1e6);
var areaImg = ee.Image.cat([
  px.updateMask(strata.eq(1)).rename('area_stratum1_jhum_km2'),
  px.updateMask(strata.eq(2)).rename('area_stratum2_near_km2'),
  px.updateMask(strata.eq(3)).rename('area_stratum3_far_km2')]);
var areaFeat = ee.Feature(null, areaImg.reduceRegion({
  reducer: ee.Reducer.sum(), geometry: roi.geometry(),
  scale: 30, maxPixels: 1e13, tileScale: 4})).set('year', YEAR);
print('STRATUM AREAS (km2)', areaFeat);

Export.table.toDrive({collection: ee.FeatureCollection([areaFeat]),
  description: 'MIZ_Jhum_V3_StratumAreas_' + tag,
  fileNamePrefix: 'MIZ_Jhum_V3_StratumAreas_' + tag,
  folder: DRIVE_DIR, fileFormat: 'CSV',
  selectors: ['year', 'area_stratum1_jhum_km2', 'area_stratum2_near_km2', 'area_stratum3_far_km2']});

//--------------------------------------------------------------
// 5. QUICK VIEW
//--------------------------------------------------------------
Map.addLayer(strata.eq(1).selfMask(), {palette: ['#ff00ff']}, 'Stratum 1: mapped jhum');
Map.addLayer(strata.eq(2).selfMask(), {palette: ['#ffaa00']}, 'Stratum 2: near jhum', false);
Map.addLayer(pts, {color: 'yellow'}, 'Validation points (blind)');
print('Points drawn', pts.size());
