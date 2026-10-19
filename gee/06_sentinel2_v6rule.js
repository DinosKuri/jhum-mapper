//==============================================================
// MIZORAM JHUM — SENTINEL-2 IMPLEMENTATION OF THE VALIDATED v6 RULE
//
// This is the Landsat v6 rule of Khiangte (2026) ported to
// Sentinel-2. It keeps every element of the validated design:
//
//   * symmetric ANNIVERSARY compositing: the same calendar months
//     on both sides of the comparison, one year apart
//   * the SAME percentile (25th) on both sides, taken over
//     per-scene NBR rather than over reflectance
//   * two branches: a wide dry-season pair and a narrow
//     burning-season pair, OR'd together
//   * a state condition (dark after clearing) AND a change
//     condition (large dNBR), not change alone
//   * explicit water, valid-observation and boundary masks
//
// WHAT IS DIFFERENT FROM THE LANDSAT VERSION, AND WHY
//   Bands   B8A (narrow NIR, 20 m) and B12 (SWIR2, 20 m).
//           B8A, not B8: B8 is a broad 10 m band whose response
//           does not match the Landsat NIR, so an index built on
//           it is not comparable with the published series.
//   Cloud   Cloud Score+ (cs_cdf >= 0.60). QA60 is NOT used:
//           ESA stopped populating it with processing baseline
//           04.00 in January 2022, so a QA60 mask silently passes
//           clouds through for most of 2022 onward.
//   Scale   20 m. B12 is a 20 m band; exporting at 10 m only
//           resamples and quadruples the file.
//   Period  2020-2025. For mapped year Y the rule needs a pre
//           window of Nov(Y-2) to May(Y-1); Level-2A is not
//           dependable over this region before late 2018, so 2020
//           is the first year with a sound pre window.
//
// ACCURACY — READ THIS BEFORE PUBLISHING ANY NUMBER FROM IT
//   On the published 2025 sample (520 points, five strata drawn
//   jointly on this map and the Landsat map; gee/07 and gee/20),
//   with area weighting, this script reached
//       OA 0.987, user's 0.739, producer's 0.811, F1 0.774
//   against
//       OA 0.981, user's 0.672, producer's 0.614, F1 0.642
//   for Landsat at 30 m. Sentinel-2 has the higher point estimate
//   on every measure; exact McNemar on the 90 discordant points
//   gives p = 0.00008, and the design-based paired differences
//   in producer's, user's, overall accuracy and F1 all exclude
//   zero (user's +0.067, 95 % CI 0.015 to 0.119).
//
//   An earlier 240-point test reported OA 0.985, user's 0.737,
//   producer's 0.688 for Sentinel-2 against user's 0.886,
//   producer's 0.820 for Landsat, and concluded that Sentinel-2
//   did not outperform Landsat. That test used the Landsat wide
//   branch on Sentinel-2 and a sample stratified on the Landsat
//   map alone. Both are superseded. Do not quote those numbers.
//
//   Sentinel-2 covers 2020 onward only. This script is the 20 m
//   companion for recent years, not a replacement for the
//   1989-2025 Landsat series.
//
// THRESHOLDS
//   The wide branch below is refitted on the 400 Sentinel-2
//   calibration points: post <= 0.27, dNBR >= 0.25. The narrow
//   branch is carried over from Landsat unchanged. These are the
//   published settings and they reproduce the accuracies above.
//   v1.0.0 of this repository shipped the Landsat wide branch
//   (0.26 / 0.20) here; it reproduces none of them. The refit
//   block in Section 7 lets you repeat the exercise with your
//   own points.
//
// OUTPUT
//   jhum_YYYY           one band per year, 1 = clearing that year
//   last_clearing       calendar year of the most recent clearing
//   years_since         years since that clearing, RIGHT-CENSORED
//   jhum_count          number of years each pixel was cleared
//   valid_years         number of years with a usable image pair
//
//   "years_since" is censored at the length of the record. A pixel
//   with no detection in 2020-2025 is NOT mature forest; it is a
//   pixel that was not cleared during six years. The jhum cycle in
//   Mizoram runs roughly 5 to 10 years, so most fallow land falls
//   in that class. Fallow age needs the 37-year Landsat stack.
//
// L. Khiangte, Mizoram University. Companion to the Landsat v6
// script MIZORAM_JHUM_FINAL_GEE.js.
//==============================================================

//--------------------------------------------------------------
// 0. SETTINGS
//--------------------------------------------------------------
var ROI_NAME    = 'Mizoram';
var ROI_ASSET   = '';          // '' = FAO GAUL level 1, as published
var START_YEAR  = 2020;
var END_YEAR    = 2025;
var DRIVE_DIR   = 'Jhum_S2_v6rule';
var EXPORT_CRS  = 'EPSG:32646';
var SCALE       = 20;          // native resolution of B8A and B12

// Cloud Score+ threshold. Higher is stricter. 0.60 is the value
// used in the published comparison.
var CS_BAND     = 'cs_cdf';
var CS_THRESH   = 0.60;

// The validated two-branch rule. Keep these unless you have
// refitted on your own interpreted points (Section 7).
// PUBLISHED SENTINEL-2 RULE. The wide branch is refitted on the same 400
// calibration points; the narrow branch is carried over from Landsat
// unchanged, because refitting it on Sentinel-2 gave no separable optimum.
// These are the thresholds behind the 634.1 km2 map and every Sentinel-2
// figure in the paper. Do not substitute the Landsat wide branch
// ({post: 0.26, dnbr: 0.20}): it was the setting of an earlier draft and
// reproduces none of the published numbers.
var RULE = {
  wide:   {post: 0.27, dnbr: 0.25},
  narrow: {post: 0.14, dnbr: 0.46},
  useNarrowBranch: true
};
// Sentinel-2 refit on 400 calibration points, for reference only:
// An alternative in which the narrow branch was also refitted
// ({post: -0.05, dnbr: 0.53}) was tested and not adopted; it is recorded
// here only so the choice is on the record.

// Composite windows as [yearOffset, month, day]
var W_PRE_S  = [-2, 11, 1], W_PRE_E  = [-1, 6, 1];   // wide, previous year
var W_POST_S = [-1, 11, 1], W_POST_E = [ 0, 6, 1];   // wide, this year
var N_PRE_S  = [-1,  3, 1], N_PRE_E  = [-1, 6, 1];   // narrow, previous year
var N_POST_S = [ 0,  3, 1], N_POST_E = [ 0, 6, 1];   // narrow, this year
var PCT = 25;                                        // same on both sides

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
// 2. SENTINEL-2 WITH CLOUD SCORE+
//--------------------------------------------------------------
var csPlus = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED');

var S2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
           .filterBounds(roi)
           .linkCollection(csPlus, [CS_BAND])
           .map(function (img) {
             img = ee.Image(img);
             var clear = img.select(CS_BAND).gte(CS_THRESH);
             var out = img.select(['B8A', 'B12'], ['NIR', 'SWIR2'])
                          .divide(10000)
                          .updateMask(clear)
                          .toFloat();
             // copyProperties returns an Element; cast back to Image
             return ee.Image(out.copyProperties(img, ['system:time_start']));
           });

function win(y, w) { return ee.Date.fromYMD(y + w[0], w[1], w[2]); }
function nbrOf(i)  { return ee.Image(i).normalizedDifference(['NIR', 'SWIR2']).rename('NBR'); }

// 25th percentile of per-scene NBR. An empty window returns a
// fully masked image rather than breaking the chain.
function pctImg(col) {
  return ee.Image(ee.Algorithms.If(col.size().gt(0),
           col.map(nbrOf).reduce(ee.Reducer.percentile([PCT])),
           ee.Image.constant(0).updateMask(0))).rename('NBR').toFloat();
}
// An empty collection has no bands, which breaks unmask(); guard it.
function safeCount(col) {
  return ee.Image(ee.Algorithms.If(col.size().gt(0),
           col.select('NIR').count(), ee.Image.constant(0)))
           .rename('n').unmask(0);
}

//--------------------------------------------------------------
// 3. THE RULE FOR ONE YEAR
//--------------------------------------------------------------
function jhumFor(y) {
  var wPreCol  = S2.filterDate(win(y, W_PRE_S),  win(y, W_PRE_E));
  var wPostCol = S2.filterDate(win(y, W_POST_S), win(y, W_POST_E));
  var nPreCol  = S2.filterDate(win(y, N_PRE_S),  win(y, N_PRE_E));
  var nPostCol = S2.filterDate(win(y, N_POST_S), win(y, N_POST_E));

  var postW = pctImg(wPostCol).rename('NBR_post_wide');
  var dW    = pctImg(wPreCol).subtract(postW).rename('dNBR_wide');
  var postN = pctImg(nPostCol).rename('NBR_post_narrow');
  var dN    = pctImg(nPreCol).subtract(postN).rename('dNBR_narrow');

  var A = postW.lte(RULE.wide.post).and(dW.gte(RULE.wide.dnbr));
  var B = postN.lte(RULE.narrow.post).and(dN.gte(RULE.narrow.dnbr));
  var jhum = RULE.useNarrowBranch ? A.or(B) : A;

  // a pixel is only assessable where BOTH sides of the wide pair
  // carry at least one clear observation
  var valid = safeCount(wPreCol).gt(0).and(safeCount(wPostCol).gt(0));

  return {
    jhum:  jhum.updateMask(valid).updateMask(WATER_FREE)
               .updateMask(ROI_MASK).rename('jhum').clip(roi),
    valid: valid.updateMask(ROI_MASK).rename('valid').clip(roi),
    idx:   ee.Image.cat([postW, dW, postN, dN])
             .updateMask(valid).updateMask(ROI_MASK).clip(roi),
    nPre:  wPreCol.size(),
    nPost: wPostCol.size()
  };
}

//--------------------------------------------------------------
// 4. ANNUAL STACK, YEARS SINCE CLEARING, CLEARING FREQUENCY
//--------------------------------------------------------------
var YEARS = [];
for (var y = START_YEAR; y <= END_YEAR; y++) YEARS.push(y);

// build each year once and reuse it: jhumFor() is expensive and the
// same year is needed by the stack, the counts and the area table
var MAPS = {};
YEARS.forEach(function (yy) { MAPS[yy] = jhumFor(yy); });

var annual = YEARS.map(function (yy) {
  return MAPS[yy].jhum.unmask(0).updateMask(ROI_MASK)
           .rename('jhum_' + yy).toByte().clip(roi);
});
var validBands = YEARS.map(function (yy) {
  return MAPS[yy].valid.unmask(0).updateMask(ROI_MASK)
           .rename('valid_' + yy).toByte().clip(roi);
});

var stack      = ee.Image.cat(annual);
var jhumCount  = stack.reduce(ee.Reducer.sum()).rename('jhum_count').toByte();
var validYears = ee.Image.cat(validBands).reduce(ee.Reducer.sum())
                   .rename('valid_years').toByte();

// most recent clearing year: the maximum year among detections
var yearImages = YEARS.map(function (yy) {
  return ee.Image.constant(yy).toInt16()
           .updateMask(MAPS[yy].jhum.eq(1)).rename('yr');
});
var lastClearing = ee.ImageCollection(yearImages).max()
                     .rename('last_clearing').toInt16();

// years since that clearing, measured from END_YEAR.
// Masked where no clearing was detected: that is censoring, not
// mature forest, and it must not be filled with a number.
var yearsSince = ee.Image.constant(END_YEAR).toInt16()
                   .subtract(lastClearing)
                   .updateMask(lastClearing.mask())
                   .rename('years_since').toInt16();

var summary = ee.Image.cat([
      lastClearing, yearsSince, jhumCount, validYears
    ]).toFloat().updateMask(ROI_MASK).clip(roi);

//--------------------------------------------------------------
// 5. DISPLAY
//--------------------------------------------------------------
var lastYear = MAPS[END_YEAR];
Map.addLayer(lastYear.jhum.selfMask(), {palette: ['#d7191c']},
             END_YEAR + ' jhum clearing (Sentinel-2, 20 m)', true);
Map.addLayer(yearsSince.clip(roi),
  {min: 0, max: (END_YEAR - START_YEAR),
   palette: ['#67000d', '#d7301f', '#fc8d59', '#fdcc8a', '#addd8e', '#31a354']},
  'Years since last clearing (censored at ' + (END_YEAR - START_YEAR) + ')', false);
Map.addLayer(jhumCount.selfMask().clip(roi),
  {min: 1, max: 3, palette: ['#ffffb2', '#fd8d3c', '#bd0026']},
  'Times cleared, ' + START_YEAR + '-' + END_YEAR, false);
Map.addLayer(validYears.clip(roi), {min: 0, max: YEARS.length,
  palette: ['#000000', '#ffffff']}, 'Years with a usable image pair', false);

// scene counts per window: observation density drives percentile depth,
// so print it before comparing one year with another
YEARS.forEach(function (yy) {
  print(yy + '  clear scenes  wide-pre / wide-post', MAPS[yy].nPre, MAPS[yy].nPost);
});

//--------------------------------------------------------------
// 6. AREA TABLE AND EXPORTS
//--------------------------------------------------------------
var px = ee.Image.pixelArea().divide(1e6);
var rows = YEARS.map(function (yy) {
  var m = MAPS[yy];
  var img = ee.Image.cat([
    px.updateMask(m.jhum.selfMask()).rename('jhum_km2'),
    px.updateMask(m.valid.selfMask()).rename('valid_km2')]);
  return ee.Feature(null, img.reduceRegion({
    reducer: ee.Reducer.sum(), geometry: roi.geometry(),
    scale: SCALE, maxPixels: 1e13, tileScale: 8}))
    .set('year', yy).set('n_pre', m.nPre).set('n_post', m.nPost);
});
var areaTable = ee.FeatureCollection(rows);
print('Annual area, Sentinel-2 (check before exporting)', areaTable.first());

Export.table.toDrive({
  collection: areaTable,
  description: 'MIZ_Jhum_S2_AnnualArea_' + START_YEAR + '_' + END_YEAR,
  fileNamePrefix: 'MIZ_Jhum_S2_AnnualArea_' + START_YEAR + '_' + END_YEAR,
  folder: DRIVE_DIR, fileFormat: 'CSV',
  selectors: ['year', 'jhum_km2', 'valid_km2', 'n_pre', 'n_post']});

Export.image.toDrive({
  image: stack.updateMask(ROI_MASK).clip(roi).toByte(),
  description: 'MIZ_Jhum_S2_AnnualStack_' + START_YEAR + '_' + END_YEAR,
  fileNamePrefix: 'MIZ_Jhum_S2_AnnualStack_' + START_YEAR + '_' + END_YEAR,
  folder: DRIVE_DIR, region: roi.geometry(), scale: SCALE,
  crs: EXPORT_CRS, maxPixels: 1e13,
  formatOptions: {cloudOptimized: true}});

Export.image.toDrive({
  image: summary,
  description: 'MIZ_Jhum_S2_Summary_' + START_YEAR + '_' + END_YEAR,
  fileNamePrefix: 'MIZ_Jhum_S2_Summary_' + START_YEAR + '_' + END_YEAR,
  folder: DRIVE_DIR, region: roi.geometry(), scale: SCALE,
  crs: EXPORT_CRS, maxPixels: 1e13,
  formatOptions: {cloudOptimized: true}});

//--------------------------------------------------------------
// 7. OPTIONAL — REFIT THE THRESHOLDS ON YOUR OWN POINTS
//    Set POINTS_ASSET to a table with columns lon, lat, ref_jhum
//    (1 = jhum for TARGET_YEAR, 0 = not). The export gives one row
//    per point with the four features; fit the thresholds outside
//    Earth Engine by grid search on a CALIBRATION set, then test
//    them once on a blind set. Never fit and report on the same
//    points.
//--------------------------------------------------------------
var POINTS_ASSET = '';          // e.g. 'projects/xxx/assets/Label_2025'
var TARGET_YEAR  = END_YEAR;

if (POINTS_ASSET) {
  var pts = ee.FeatureCollection(POINTS_ASSET).map(function (f) {
    var x = ee.Number.parse(ee.Algorithms.String(f.get('lon')));
    var yy = ee.Number.parse(ee.Algorithms.String(f.get('lat')));
    return ee.Feature(ee.Geometry.Point([x, yy]), f.toDictionary());
  });
  var feats = MAPS[TARGET_YEAR].idx.sampleRegions({
    collection: pts,
    properties: ['point_id', 'ref_jhum', 'ref_confidence', 'ref_landcover', 'notes'],
    scale: SCALE, tileScale: 8, geometries: false});
  Export.table.toDrive({
    collection: feats,
    description: 'MIZ_Jhum_S2_Features_' + TARGET_YEAR,
    fileNamePrefix: 'MIZ_Jhum_S2_Features_' + TARGET_YEAR,
    folder: DRIVE_DIR, fileFormat: 'CSV',
    selectors: ['point_id', 'ref_jhum', 'ref_confidence', 'ref_landcover',
                'notes', 'NBR_post_wide', 'dNBR_wide',
                'NBR_post_narrow', 'dNBR_narrow']});
}

//--------------------------------------------------------------
// 8. CLICK INSPECTOR — use only after labelling is finished
//--------------------------------------------------------------
var panel = ui.Panel({style: {position: 'bottom-left', width: '330px', padding: '8px'}});
panel.add(ui.Label('Click a pixel', {fontWeight: 'bold'}));
Map.add(panel);
Map.style().set('cursor', 'crosshair');
Map.onClick(function (c) {
  panel.clear();
  panel.add(ui.Label('Click a pixel', {fontWeight: 'bold'}));
  var pt = ee.Geometry.Point([c.lon, c.lat]);
  ee.Image.cat([lastYear.idx, lastYear.jhum.unmask(0), yearsSince.unmask(-1),
                jhumCount.unmask(0)])
    .reduceRegion({reducer: ee.Reducer.first(), geometry: pt, scale: SCALE})
    .evaluate(function (v) {
      if (!v) { panel.add(ui.Label('no data here')); return; }
      panel.add(ui.Label(END_YEAR + ' NBR_post wide: ' +
        (v.NBR_post_wide === null ? 'masked' : v.NBR_post_wide.toFixed(3))));
      panel.add(ui.Label(END_YEAR + ' dNBR wide:     ' +
        (v.dNBR_wide === null ? 'masked' : v.dNBR_wide.toFixed(3))));
      panel.add(ui.Label(END_YEAR + ' jhum: ' + v.jhum));
      panel.add(ui.Label('years since clearing: ' +
        (v.years_since < 0 ? 'none detected in record' : v.years_since)));
      panel.add(ui.Label('times cleared: ' + v.jhum_count));
    });
});
