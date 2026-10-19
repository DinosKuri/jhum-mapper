//==============================================================
// MIZORAM JHUM — DUAL-MAP VALIDATION SAMPLE (Landsat and Sentinel-2)
//
// WHY THIS EXISTS
// The published Sentinel-2 comparison used a sample stratified on
// the LANDSAT map. That design can measure how often Sentinel-2
// agrees with Landsat where Landsat looked, but it cannot measure
// what Sentinel-2 finds on its own: the 13,473 km2 "beyond 300 m"
// stratum carries one point per 135 km2, so Sentinel-2-only
// detections there are effectively unsampled. No defensible
// accuracy statement about the Sentinel-2 MAP can come from it.
//
// WHAT THIS DOES
// Draws ONE stratified sample whose strata are built from BOTH
// maps at once, so a single set of interpreted points yields:
//   * area-weighted accuracy for the Landsat map
//   * area-weighted accuracy for the Sentinel-2 map
//   * error-adjusted area for each
//   * the area of disagreement, with a variance
//   * a paired (McNemar) test of whether the sensors differ
// Strata different from map classes are admissible under
// Stehman (2014); estimators follow Olofsson et al. (2014).
//
// STRATA
//   A  both maps say jhum                     (agreement, positive)
//   B  Landsat only                           (disagreement)
//   C  Sentinel-2 only                        (disagreement)
//   D  neither, but within 300 m of either    (near miss)
//   E  neither, beyond 300 m                  (the rest of the state)
// B and C are small in area and carry the scientific question, so
// they are deliberately oversampled. Weighting restores the areas.
//
// OUTPUT
//   ..._BLIND  shuffled points, no map class  -> interpret this one
//   ..._KEY    the same points with map classes and features
//   ..._STRATA stratum areas in km2, for the weights
// Interpret the BLIND file only. Opening the KEY first destroys
// the independence of the assessment.
//
// L. Khiangte, Mizoram University.
//==============================================================

//--------------------------------------------------------------
// 0. SETTINGS
//--------------------------------------------------------------
var YEAR      = 2025;
var ROI_NAME  = 'Mizoram';
var DRIVE_DIR = 'Jhum_dualmap_validation';
var SEED      = 20261013;   // new seed: the first draw's KEY has been seen,
                            // so those points can no longer be interpreted blind
var BUFFER_M  = 300;

// Points per stratum, set from the areas measured on 12 October 2026:
//   A both        425.6 km2      B Landsat only   102.5 km2
//   C S2 only     208.5 km2      D near, neither 8388.3 km2
//   E rest      12059.2 km2
// D carries most of the omission risk for both maps and had only
// 60 points in the first draw, one per 140 km2, worse than the 51
// km2 per point of the published 2025 design. It is raised to 100.
// Total 320.
var N = {A: 40, B: 45, C: 45, D: 100, E: 90};

// Each sensor carries its OWN thresholds, fitted the same way on the
// same 400 calibration points. Fitting Landsat thresholds and then
// applying them to Sentinel-2 would test threshold transfer, not the
// sensor. A bootstrap over those points (2,000 resamples) puts the
// Sentinel-2 minus Landsat difference at -0.04 for the state
// threshold (95 % CI -0.15 to +0.04) and +0.01 for the change
// threshold (-0.08 to +0.12), so the two coincide within sampling
// error; they are kept separate anyway, because symmetry of
// treatment is what makes the comparison fair.
var RULE_L  = {wide: {post: 0.26, dnbr: 0.20}, narrow: {post: 0.14, dnbr: 0.46}};
var RULE_S2 = {wide: {post: 0.27, dnbr: 0.25}, narrow: {post: 0.14, dnbr: 0.46}};
var W_PRE_S  = [-2, 11, 1], W_PRE_E  = [-1, 6, 1];
var W_POST_S = [-1, 11, 1], W_POST_E = [ 0, 6, 1];
var N_PRE_S  = [-1,  3, 1], N_PRE_E  = [-1, 6, 1];
var N_POST_S = [ 0,  3, 1], N_POST_E = [ 0, 6, 1];
var PCT = 25;
var CS_THRESH = 0.60;

//--------------------------------------------------------------
// 1. AREA AND MASKS
//--------------------------------------------------------------
var roi = ee.FeatureCollection('FAO/GAUL/2015/level1')
            .filter(ee.Filter.eq('ADM1_NAME', ROI_NAME));
var ROI_MASK = ee.Image.constant(1).clip(roi).mask();
var WATER_FREE = ee.Image('JRC/GSW1_4/GlobalSurfaceWater')
                   .select('occurrence').unmask(0).lt(50);
Map.centerObject(roi, 8);

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

//--------------------------------------------------------------
// 2. LANDSAT MAP (the published rule)
//--------------------------------------------------------------
var OLD = ['SR_B4', 'SR_B7'], NEW = ['SR_B5', 'SR_B7'];
function prepL(img, bands, harmonize) {
  img = ee.Image(img);
  var qa = img.select('QA_PIXEL');
  var clear = qa.bitwiseAnd(1 << 1).eq(0)
              .and(qa.bitwiseAnd(1 << 3).eq(0))
              .and(qa.bitwiseAnd(1 << 4).eq(0));
  var out = img.select(bands, ['NIR', 'SWIR2'])
               .multiply(0.0000275).add(-0.2).updateMask(clear);
  if (harmonize) out = out.multiply([0.8462, 0.9071]).add([0.0412, 0.0172]);
  return ee.Image(out.toFloat().copyProperties(img, ['system:time_start']));
}
var LS = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2').filterBounds(roi)
           .map(function (i) { return prepL(i, OLD, true); })
  .merge(ee.ImageCollection('LANDSAT/LE07/C02/T1_L2').filterBounds(roi)
           .map(function (i) { return prepL(i, OLD, true); }))
  .merge(ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filterBounds(roi)
           .map(function (i) { return prepL(i, NEW, false); }))
  .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2').filterBounds(roi)
           .map(function (i) { return prepL(i, NEW, false); }));

//--------------------------------------------------------------
// 3. SENTINEL-2 MAP (same rule, Cloud Score+, B8A / B12)
//--------------------------------------------------------------
var csPlus = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED');
var S2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
           .filterBounds(roi).linkCollection(csPlus, ['cs_cdf'])
           .map(function (img) {
             img = ee.Image(img);
             var clear = img.select('cs_cdf').gte(CS_THRESH);
             var out = img.select(['B8A', 'B12'], ['NIR', 'SWIR2'])
                          .divide(10000).updateMask(clear).toFloat();
             return ee.Image(out.copyProperties(img, ['system:time_start']));
           });

// MIN_SCENES stays at 1, which is the definition used by the
// published rule. The first draw threw up four points where a
// Sentinel-2 percentile was computed from so few clear scenes that
// NBR reached -1.000, a degenerate value. Rather than change the
// rule, the clear-scene counts are carried into the KEY so the
// effect can be measured after interpretation: raise MIN_SCENES
// then and report how much area it removes.
var MIN_SCENES_L  = 1;
var MIN_SCENES_S2 = 1;

function mapFor(col, y, tag, RULE, minScenes) {
  var wPre = col.filterDate(win(y, W_PRE_S),  win(y, W_PRE_E));
  var wPost= col.filterDate(win(y, W_POST_S), win(y, W_POST_E));
  var nPre = col.filterDate(win(y, N_PRE_S),  win(y, N_PRE_E));
  var nPost= col.filterDate(win(y, N_POST_S), win(y, N_POST_E));
  var postW = pctImg(wPost), dW = pctImg(wPre).subtract(postW);
  var postN = pctImg(nPost), dN = pctImg(nPre).subtract(postN);
  var A = postW.lte(RULE.wide.post).and(dW.gte(RULE.wide.dnbr));
  var B = postN.lte(RULE.narrow.post).and(dN.gte(RULE.narrow.dnbr));
  var nPreI  = safeCount(wPre).rename('n_pre_'  + tag);
  var nPostI = safeCount(wPost).rename('n_post_' + tag);
  var valid = nPreI.gte(minScenes).and(nPostI.gte(minScenes));
  return {
    jhum: A.or(B).updateMask(valid).updateMask(WATER_FREE)
            .updateMask(ROI_MASK).rename('map_' + tag).clip(roi),
    idx: ee.Image.cat([postW.rename('post_wide_' + tag),
                       dW.rename('dnbr_wide_' + tag),
                       postN.rename('post_narrow_' + tag),
                       dN.rename('dnbr_narrow_' + tag),
                       nPreI, nPostI])
           .updateMask(valid).updateMask(ROI_MASK).clip(roi)
  };
}
var L = mapFor(LS, YEAR, 'landsat', RULE_L,  MIN_SCENES_L);
var S = mapFor(S2, YEAR, 's2',      RULE_S2, MIN_SCENES_S2);

// Both maps are put on the same 30 m grid so the strata are well
// defined, and Sentinel-2 is aggregated by majority within each
// Landsat cell.
//
// setDefaultProjection() is required here. Images built by reducing
// an ImageCollection are computed objects with no inherent
// projection, and reduceResolution refuses to run on one: it has to
// know the input cell size before it can average cells together.
// Declaring the native grid of each sensor first makes the
// aggregation meaningful rather than arbitrary.
var PROJ   = ee.Projection('EPSG:32646').atScale(30);   // Landsat grid
var S2PROJ = ee.Projection('EPSG:32646').atScale(20);   // B8A / B12 grid

var sOn30 = S.jhum.unmask(0)
             .setDefaultProjection(S2PROJ)
             .reduceResolution({reducer: ee.Reducer.mean(), maxPixels: 36})
             .reproject(PROJ)
             .gte(0.5).rename('map_s2_30m');

var lJ = L.jhum.unmask(0)
          .setDefaultProjection(PROJ)
          .reproject(PROJ).rename('map_landsat_30m');

//--------------------------------------------------------------
// 4. THE FIVE STRATA
//--------------------------------------------------------------
var either = lJ.or(sOn30).rename('either');

// distance in metres to the nearest mapped pixel, from either map.
// A Euclidean kernel states the radius in metres directly, which is
// safer than converting a pixel-unit distance transform.
var near = either.selfMask()
             .distance(ee.Kernel.euclidean(BUFFER_M, 'meters'), false)
             .lte(BUFFER_M).unmask(0).reproject(PROJ);

var strata = ee.Image(0)
  .where(lJ.eq(1).and(sOn30.eq(1)), 1)                 // A both
  .where(lJ.eq(1).and(sOn30.eq(0)), 2)                 // B Landsat only
  .where(lJ.eq(0).and(sOn30.eq(1)), 3)                 // C Sentinel-2 only
  .where(either.eq(0).and(near), 4)                    // D near miss
  .where(either.eq(0).and(near.not()), 5)              // E rest
  .rename('stratum').updateMask(ROI_MASK).clip(roi).toByte()
  .setDefaultProjection(PROJ);   // ee.Image(0) carries no projection of its own

// The strata layer is DISPLAY ONLY and is off by default.
// Turning it on at state-wide zoom asks Earth Engine to compute a
// 300 m Euclidean distance kernel on a reprojected 30 m grid across
// 21 000 km2 for every map tile, which overruns the interactive
// memory limit ("Object too large", "User memory limit exceeded").
// Those are tile errors in the map pane. They do NOT affect the
// export tasks, which run on the batch backend. If you want to see
// the strata, tick the layer on after zooming in to a district.
var SHOW_STRATA = false;
Map.addLayer(strata.selfMask(), {min: 1, max: 5,
  palette: ['#b2182b', '#ef8a62', '#67a9cf', '#d9d9d9', '#f7f7f7']},
  'strata: A both, B Landsat only, C S2 only, D near, E rest', SHOW_STRATA);

var areaKm2 = ee.Image.pixelArea().divide(1e6).addBands(strata)
  .reduceRegion({reducer: ee.Reducer.sum().group(1, 'stratum'),
                 geometry: roi.geometry(), scale: 30,
                 maxPixels: 1e13, tileScale: 8});
print('STRATUM AREAS (km2) — read these before sampling', areaKm2);

//--------------------------------------------------------------
// 5. DRAW THE SAMPLE
//--------------------------------------------------------------
var classValues = [1, 2, 3, 4, 5];
var classPoints = [N.A, N.B, N.C, N.D, N.E];

var sample = strata.stratifiedSample({
  numPoints: 0, classBand: 'stratum', region: roi.geometry(),
  scale: 30, seed: SEED, classValues: classValues,
  classPoints: classPoints, geometries: true, tileScale: 8
});

// attach map classes and both sets of features, then shuffle
var feats = ee.Image.cat([lJ.rename('map_landsat'), sOn30.rename('map_s2'),
                          L.idx, S.idx]);
var keyed = feats.sampleRegions({collection: sample, properties: ['stratum'],
                                 scale: 30, tileScale: 8, geometries: true})
  .randomColumn('shuffle', SEED)
  .sort('shuffle');

var ids = ee.List.sequence(1, keyed.size());
var keyList = keyed.toList(keyed.size());
var tagged = ee.FeatureCollection(ids.map(function (i) {
  i = ee.Number(i);
  var f = ee.Feature(keyList.get(i.subtract(1)));
  var pid = ee.String('D2025_').cat(ee.Number(i).format('%03d'));
  var g = f.geometry().coordinates();
  return f.set('point_id', pid)
          .set('lon', ee.Number(ee.List(g).get(0)))
          .set('lat', ee.Number(ee.List(g).get(1)));
}));

// BLIND: coordinates and empty label columns only
var blind = tagged.map(function (f) {
  return ee.Feature(null, {
    point_id: f.get('point_id'), lon: f.get('lon'), lat: f.get('lat'),
    ref_jhum: '', ref_confidence: '', ref_landcover: '', notes: ''});
});

Export.table.toDrive({
  collection: blind,
  description: 'MIZ_Jhum_D1_BLIND_dualmap_' + YEAR,
  fileNamePrefix: 'MIZ_Jhum_D1_BLIND_dualmap_' + YEAR,
  folder: DRIVE_DIR, fileFormat: 'CSV',
  selectors: ['point_id', 'lon', 'lat', 'ref_jhum', 'ref_confidence',
              'ref_landcover', 'notes']});

Export.table.toDrive({
  collection: tagged,
  description: 'MIZ_Jhum_D2_KEY_dualmap_' + YEAR,
  fileNamePrefix: 'MIZ_Jhum_D2_KEY_dualmap_' + YEAR,
  folder: DRIVE_DIR, fileFormat: 'CSV',
  selectors: ['point_id', 'lon', 'lat', 'stratum', 'map_landsat', 'map_s2',
              'post_wide_landsat', 'dnbr_wide_landsat',
              'post_narrow_landsat', 'dnbr_narrow_landsat',
              'n_pre_landsat', 'n_post_landsat',
              'post_wide_s2', 'dnbr_wide_s2',
              'post_narrow_s2', 'dnbr_narrow_s2',
              'n_pre_s2', 'n_post_s2']});

// stratum areas as a table, for the weights
var areaRows = ee.List(ee.Dictionary(areaKm2).get('groups')).map(function (g) {
  g = ee.Dictionary(g);
  return ee.Feature(null, {stratum: g.get('stratum'), area_km2: g.get('sum')});
});
Export.table.toDrive({
  collection: ee.FeatureCollection(areaRows),
  description: 'MIZ_Jhum_D3_StratumAreas_dualmap_' + YEAR,
  fileNamePrefix: 'MIZ_Jhum_D3_StratumAreas_dualmap_' + YEAR,
  folder: DRIVE_DIR, fileFormat: 'CSV',
  selectors: ['stratum', 'area_km2']});

//--------------------------------------------------------------
// 6. THE AGREEMENT RASTER  (Figure 7c and Figure 9c)
//--------------------------------------------------------------
// One byte band on the Landsat 30 m grid, masked everywhere
// neither map claims a clearing, so the GeoTIFF carries only the
// three classes a reader has to tell apart:
//
//    1  Sentinel-2 only     #1b9e77   green
//    2  Landsat only        #e08214   orange
//    3  both maps           #2c4a63   dark blue
//
// The arithmetic is lJ*2 + sOn30, which is why the order runs
// S2, Landsat, both rather than alphabetically.
var agree = lJ.multiply(2).add(sOn30)
             .updateMask(lJ.or(sOn30))
             .rename('agreement').toByte()
             .clip(roi);

Map.addLayer(agree, {min: 1, max: 3,
  palette: ['#1b9e77', '#e08214', '#2c4a63']},
  'agreement: 1 S2 only, 2 Landsat only, 3 both', false);

// Areas of the three classes, printed so the figure caption can
// quote them without a second run.
print('AGREEMENT CLASS AREAS (km2)',
  ee.Image.pixelArea().divide(1e6).addBands(agree)
    .reduceRegion({reducer: ee.Reducer.sum().group(1, 'agreement'),
                   geometry: roi.geometry(), scale: 30,
                   maxPixels: 1e13, tileScale: 8}));

// GeoTIFF for QGIS. 30 m, UTM 46N, so it drops straight on top of
// the Figure 9 layers with no reprojection on the fly.
Export.image.toDrive({
  image: agree.unmask(0).toByte(),
  description: 'MIZ_Jhum_Agreement_L8S2_' + YEAR,
  fileNamePrefix: 'MIZ_Jhum_Agreement_L8S2_' + YEAR,
  folder: DRIVE_DIR,
  region: roi.geometry(),
  scale: 30,
  crs: 'EPSG:32646',
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF',
  formatOptions: {cloudOptimized: true, noData: 0}});

//--------------------------------------------------------------
// 7. HOW TO USE THE RESULT
//--------------------------------------------------------------
// 1. Interpret the BLIND file on PlanetScope with the existing
//    protocol (Text S1). Do not open the KEY until every point is
//    labelled. Use the blind viewer, not this script.
// 2. Join BLIND to KEY on point_id.
// 3. For each map separately, build the stratified estimator of
//    Olofsson et al. (2014) with the stratum areas from D3 as the
//    weights, and report user's, producer's and error-adjusted area
//    with 95 % intervals.
// 4. For the sensor comparison, count the discordant pairs:
//       b = points where Landsat is right and Sentinel-2 wrong
//       c = points where Sentinel-2 is right and Landsat wrong
//    and apply an exact McNemar test. That is the paired test the
//    question needs; two independent accuracy figures with
//    overlapping intervals do not answer it.
// 5. Report the area of strata B and C, which is the mapped
//    disagreement, with its own confidence interval.
