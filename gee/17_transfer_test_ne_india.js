//==============================================================
// JHUM RULE — TRANSFER TEST ACROSS NORTH-EAST INDIA
//
// Applies the Mizoram-calibrated rule, unchanged, to neighbouring
// states, with the burning season shifted to the local calendar.
//
// WHAT THIS CAN AND CANNOT SHOW
//   It CAN show whether the rule produces plausible, spatially
//   coherent results elsewhere, and whether detections coincide
//   with independent fire observations.
//   It CANNOT give accuracy. User's and producer's accuracy need
//   reference points labelled in each state.
//
// Independent checks included (no labelling required):
//   1. mapped area and share of each state
//   2. share of mapped pixels within 1 km of a dry-season FIRMS
//      fire detection, against the same share for random forest
//      pixels (the control)
//   3. elevation and slope distribution of mapped pixels
//
// Rule (unchanged from Mizoram):
//   A: NBR_post_wide   <= 0.26  AND  dNBR_wide   >= 0.20
//   B: NBR_post_narrow <= 0.14  AND  dNBR_narrow >= 0.46
//==============================================================

//--------------------------------------------------------------
// 0. SETTINGS
//--------------------------------------------------------------
var YEAR      = 2025;
var DRIVE_DIR = 'Jhum_transfer_NE_v2';

// Burn calendar per state. [startMonth, endMonth] of the BURNING
// season; the wide window runs from CLEAR_START to BURN_END+1.
// These are starting values from the general literature: verify
// each one locally before drawing conclusions from that state.
var STATES = [
  // clearStart / burnEnd are the months bounding the WIDE window.
  // Sources: Manipur, slashing Nov-Dec, burning March, sowing Apr-May
  //   (Kurien et al., Remote Sensing Applications 2019);
  // Meghalaya (Khasi), jungle cutting Dec-Jan, burning Feb-Mar
  //   (NEHU Journal XVII); East Garo Hills, burning Feb-Mar, cropping Apr-May;
  // Region-wide, plot selection Oct-Dec and firing Feb-Mar
  //   (NIRD-NERC, Shifting Cultivation: Towards Transformation Approach);
  // Arunachal Pradesh, burning Feb-Mar (Shifting Cultivation in NE India, 2019);
  // Mizoram, felling from Nov, drying by Chapchar Kut (late Feb to early Mar),
  //   burning Feb-Apr: this study's own interpretation found 44 of 63 plots
  //   cleared in November.
  {name: 'Mizoram',           clearStart: 11, burnEnd: 6, note: 'reference case'},
  {name: 'Manipur',           clearStart: 11, burnEnd: 5, note: 'slash Nov-Dec, burn Mar, sow Apr-May'},
  {name: 'Nagaland',          clearStart: 11, burnEnd: 5, note: 'clearing in winter, burn Feb-Mar'},
  {name: 'Meghalaya',         clearStart: 12, burnEnd: 4, note: 'cut Dec-Jan, burn Feb-Mar, crop Apr-May'},
  {name: 'Tripura',           clearStart: 12, burnEnd: 5, note: 'burn Feb-Apr'},
  {name: 'Arunachal Pradesh', clearStart: 11, burnEnd: 4, note: 'burn Feb-Mar'},
  {name: 'Assam',             clearStart: 11, burnEnd: 5, note: 'mostly plains: expect very little'}
];

var RULE = {wide: {post: 0.26, dnbr: 0.20}, narrow: {post: 0.14, dnbr: 0.46}};
var PCT = 25;
var FIRE_RADIUS = 500;    // metres; 1 km saturated the control in round 1
var SCALE       = 30;     // round 1 used 120 m, which under-counted small plots
var MIN_SLOPE   = 10;     // degrees; jhum is a hill-slope practice
var MIN_TREE    = 40;     // per cent tree cover in 2000
var USE_TERRAIN_MASK = true;   // the round-1 result showed this is needed

//--------------------------------------------------------------
// 1. HARMONIZED LANDSAT
//--------------------------------------------------------------
var ADM = ee.FeatureCollection('FAO/GAUL/2015/level1')
            .filter(ee.Filter.eq('ADM0_NAME', 'India'));

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
function landsat(region) {
  return ee.ImageCollection('LANDSAT/LT05/C02/T1_L2').filterBounds(region)
           .map(function (i) { return prep(i, 'SR_B4', 'SR_B7', true); })
    .merge(ee.ImageCollection('LANDSAT/LE07/C02/T1_L2').filterBounds(region)
           .map(function (i) { return prep(i, 'SR_B4', 'SR_B7', true); }))
    .merge(ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filterBounds(region)
           .map(function (i) { return prep(i, 'SR_B5', 'SR_B7', false); }))
    .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2').filterBounds(region)
           .map(function (i) { return prep(i, 'SR_B5', 'SR_B7', false); }));
}
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

var WATER_FREE = ee.Image('JRC/GSW1_4/GlobalSurfaceWater')
                   .select('occurrence').unmask(0).lt(50);
var SRTM = ee.Image('USGS/SRTMGL1_003');
var SLOPE = ee.Terrain.slope(SRTM);
var FOREST = ee.Image('UMD/hansen/global_forest_change_2023_v1_11')
               .select('treecover2000').gte(40);

//--------------------------------------------------------------
// 2. THE RULE, WITH A LOCAL CALENDAR
//--------------------------------------------------------------
function jhumFor(st, region) {
  var LS = landsat(region);
  var cs = st.clearStart, be = st.burnEnd;
  function d(yOff, m) { return ee.Date.fromYMD(YEAR + yOff, m, 1); }
  var wPreC  = LS.filterDate(d(-2, cs), d(-1, be));
  var wPostC = LS.filterDate(d(-1, cs), d(0, be));
  var nPreC  = LS.filterDate(d(-1, be - 3), d(-1, be));
  var nPostC = LS.filterDate(d(0, be - 3), d(0, be));

  var pW = pctImg(wPostC), dW = pctImg(wPreC).subtract(pW);
  var pN = pctImg(nPostC), dN = pctImg(nPreC).subtract(pN);
  var A = pW.lte(RULE.wide.post).and(dW.gte(RULE.wide.dnbr));
  var B = pN.lte(RULE.narrow.post).and(dN.gte(RULE.narrow.dnbr));
  var valid = safeCount(wPreC).gt(0).and(safeCount(wPostC).gt(0));
  var jh = A.or(B).updateMask(valid).updateMask(WATER_FREE);
  if (USE_TERRAIN_MASK) jh = jh.updateMask(SLOPE.gte(MIN_SLOPE)).updateMask(FOREST);
  return {jhum: jh.rename('Jhum').clip(region),
          valid: valid.clip(region),
          nPre: wPreC.size(), nPost: wPostC.size()};
}

//--------------------------------------------------------------
// 3. INDEPENDENT CHECK: FIRMS DRY-SEASON FIRE DETECTIONS
//--------------------------------------------------------------
function fireProximity(st, region) {
  // restricted to the BURNING months only, which is sharper than the
  // whole dry season used in round 1
  var fires = ee.ImageCollection('FIRMS')
    .filterDate(ee.Date.fromYMD(YEAR, Math.max(1, st.burnEnd - 4), 1),
                ee.Date.fromYMD(YEAR, st.burnEnd, 1))
    .select('T21').max().gt(0).unmask(0).clip(region);
  return fires.focalMax({radius: FIRE_RADIUS, units: 'meters'}).rename('near_fire');
}

// second independent comparator: MODIS burned area for the same months
function burnedArea(st, region) {
  return ee.ImageCollection('MODIS/061/MCD64A1')
    .filterDate(ee.Date.fromYMD(YEAR, Math.max(1, st.burnEnd - 4), 1),
                ee.Date.fromYMD(YEAR, st.burnEnd, 1))
    .select('BurnDate').max().gt(0).unmask(0).clip(region).rename('burned');
}

//--------------------------------------------------------------
// 4. RUN EACH STATE
//--------------------------------------------------------------
var px = ee.Image.pixelArea().divide(1e6);
var rows = [];

STATES.forEach(function (st) {
  var fc = ADM.filter(ee.Filter.eq('ADM1_NAME', st.name));
  var region = fc.geometry();
  var m = jhumFor(st, region);
  var nearFire = fireProximity(st, region);
  var burned = burnedArea(st, region);

  // area and share
  var areaImg = ee.Image.cat([
    px.updateMask(m.jhum.selfMask()).rename('jhum_km2'),
    px.updateMask(m.valid.selfMask()).rename('valid_km2'),
    px.rename('state_km2'),
    px.updateMask(m.jhum.selfMask().and(nearFire.eq(1))).rename('jhum_near_fire_km2'),
    px.updateMask(FOREST.selfMask().and(nearFire.eq(1))).rename('forest_near_fire_km2'),
    px.updateMask(FOREST.selfMask()).rename('forest_km2'),
    px.updateMask(m.jhum.selfMask().and(burned.eq(1))).rename('jhum_burned_km2'),
    px.updateMask(FOREST.selfMask().and(burned.eq(1))).rename('forest_burned_km2'),
    SRTM.updateMask(m.jhum.selfMask()).rename('elev_mean'),
    SLOPE.updateMask(m.jhum.selfMask()).rename('slope_mean')
  ]);

  var sums = areaImg.select(['jhum_km2','valid_km2','state_km2','jhum_near_fire_km2',
                             'forest_near_fire_km2','forest_km2','jhum_burned_km2',
                             'forest_burned_km2'])
    .reduceRegion({reducer: ee.Reducer.sum(), geometry: region, scale: SCALE,
                   maxPixels: 1e13, tileScale: 8, bestEffort: true});
  var means = areaImg.select(['elev_mean','slope_mean'])
    .reduceRegion({reducer: ee.Reducer.mean(), geometry: region, scale: SCALE,
                   maxPixels: 1e13, tileScale: 8, bestEffort: true});

  rows.push(ee.Feature(null, sums.combine(means))
    .set('state', st.name).set('year', YEAR)
    .set('clear_start_month', st.clearStart).set('burn_end_month', st.burnEnd)
    .set('n_scenes_pre', m.nPre).set('n_scenes_post', m.nPost));

  Map.addLayer(m.jhum.selfMask(), {palette: ['#ff00ff']}, st.name + ' jhum ' + YEAR, false);
});

var table = ee.FeatureCollection(rows);
print('Mizoram row (the reference case: check this first)', rows[0]);
print('Interpretation: jhum_near_fire_km2 / jhum_km2 is the share of mapped clearing',
      'near a dry-season fire detection. Compare it with forest_near_fire_km2 / forest_km2,',
      'which is the background rate. A much higher share for mapped jhum supports the result;',
      'a similar share does not.');

Export.table.toDrive({
  collection: table,
  description: 'JHUM_TransferTest_NEIndia_v2_' + YEAR,
  fileNamePrefix: 'JHUM_TransferTest_NEIndia_v2_' + YEAR,
  folder: DRIVE_DIR, fileFormat: 'CSV',
  selectors: ['state','year','clear_start_month','burn_end_month',
              'jhum_km2','valid_km2','state_km2',
              'jhum_near_fire_km2','forest_near_fire_km2','forest_km2',
              'jhum_burned_km2','forest_burned_km2',
              'elev_mean','slope_mean','n_scenes_pre','n_scenes_post']
});

Map.centerObject(ADM.filter(ee.Filter.eq('ADM1_NAME', 'Nagaland')), 7);
