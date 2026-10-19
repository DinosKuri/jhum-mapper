//==============================================================
// MIZORAM JHUM DETECTION — PART A (+ Part B)
//
// PART A : turns your jhum POLYGONS into a SAMPLE TABLE
//          (one row per pixel: year, jhum 1/0, dNBR, RBR ...)
// PART B : (later) applies the fitted equation and exports maps
//
// INPUT  : your polygon asset, class_id 1–7 = 1990 … 2025
// OUTPUT : 1. GEE asset  MIZ_Jhum_Samples_Landsat_1990_2025_v1  (for Part C)
//          2. Drive CSV  MIZ_Jhum_Samples_Landsat_1990_2025_v1  (for Colab)
//==============================================================

//--------------------------------------------------------------
// 0. SETTINGS  — edit only this block
//--------------------------------------------------------------
var POLYGON_ASSET = 'projects/phdmrsei89448/assets/Jhum_1990toJhum_2025'; // <-- your polygons   // PRIVATE ASSET — see README, Private Earth Engine assets
var OUT_FOLDER    = 'projects/phdmrsei89448/assets/';                     // where the sample table is saved   // PRIVATE ASSET — see README, Private Earth Engine assets
// GEE writes to a Drive FOLDER NAME, not a path (no nested folders).
var DRIVE_DIR     = 'Jhum_samples_v3';   // samples CSV
var MAP_DIR       = 'Jhum_maps_v3';      // yearly jhum maps + area table
var VERSION       = 'v3';   // v3 = anniversary (same-season) 'before' image

var BG_PER_YEAR = 3000;   // background (non-jhum) pixels per year
var BG_BUFFER_M = 500;    // background kept this far away from that year's jhum polygons
var SEED        = 42;     // fixed seed = same samples every time you run

// class_id -> year (matches your polygons)
var CLASS_YEAR = {1: 1990, 2: 2000, 3: 2005, 4: 2010, 5: 2015, 6: 2020, 7: 2025};

// Windows: [yearOffset, month, day]  (yearOffset -1 = previous year)
// v3 uses an ANNIVERSARY 'before' image: the SAME dry season (Mar–May),
// one year earlier, when this year's jhum plot was still fallow vegetation.
// Same season in both images = no false change from seasonal drying.
var PRE_START  = [-1, 3, 1];   // 1 Mar, previous year
var PRE_END    = [-1, 6, 1];   // 31 May, previous year (end exclusive)
var POST_START = [0, 3, 1];    // 1 Mar, this year (after burning)
var POST_END   = [0, 6, 1];    // 31 May, this year

// Percentile of per-image NBR in each window
//   PRE_PCT  = 75 -> typical-to-green fallow state (ignores haze/smoke dates)
//   POST_PCT = 25 -> the most burned state (works even with few images)
var PRE_PCT  = 75;
var POST_PCT = 25;

var SAMPLE_NAME = 'MIZ_Jhum_Samples_Landsat_1990_2025_' + VERSION;

//--------------------------------------------------------------
// 1. STUDY AREA, POLYGONS, MASKS
//--------------------------------------------------------------
var roi = ee.FeatureCollection('FAO/GAUL/2015/level1')
            .filter(ee.Filter.eq('ADM1_NAME', 'Mizoram'));
Map.centerObject(roi, 8);

var ref = ee.FeatureCollection(POLYGON_ASSET);

print('CHECK 1 – Polygons per class_id (expect 7 classes)', ref.aggregate_histogram('class_id'));
print('CHECK 2 – Polygons per jhum_year', ref.aggregate_histogram('jhum_year'));

var WATER_FREE = ee.Image('JRC/GSW1_4/GlobalSurfaceWater')
                   .select('occurrence').unmask(0).lt(50);
var TERRAIN = ee.Terrain.products(ee.Image('USGS/SRTMGL1_003'))
                .select(['elevation', 'slope', 'aspect']);

//--------------------------------------------------------------
// 2. HARMONIZED LANDSAT 5/7/8/9 (Collection 2, Level-2)
//--------------------------------------------------------------
function prepL(img, nir, swir2) {
  img = ee.Image(img);
  var qa = img.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(1 << 1).eq(0)        // dilated cloud
             .and(qa.bitwiseAnd(1 << 3).eq(0))  // cloud
             .and(qa.bitwiseAnd(1 << 4).eq(0)); // cloud shadow
  return ee.Image(img.select([nir, swir2], ['NIR', 'SWIR2'])
            .multiply(0.0000275).add(-0.2)
            .updateMask(mask)
            .toFloat()
            .copyProperties(img, ['system:time_start']));
}

// Roy et al. (2016): TM/ETM+ -> OLI (NIR, SWIR2)
function toOLI(img) {
  img = ee.Image(img);
  return ee.Image(img.multiply([0.8462, 0.9071]).add([0.0412, 0.0172])
            .toFloat()
            .copyProperties(img, ['system:time_start']));
}

var L5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2').filterBounds(roi)
           .map(function (i) { return toOLI(prepL(i, 'SR_B4', 'SR_B7')); });
var L7 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2').filterBounds(roi)
           .map(function (i) { return toOLI(prepL(i, 'SR_B4', 'SR_B7')); });
var L8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filterBounds(roi)
           .map(function (i) { return prepL(i, 'SR_B5', 'SR_B7'); });
var L9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2').filterBounds(roi)
           .map(function (i) { return prepL(i, 'SR_B5', 'SR_B7'); });

var LS = L5.merge(L7).merge(L8).merge(L9);

//--------------------------------------------------------------
// 3. YEARLY FEATURES (NBR before / after, dNBR, RBR, terrain)
//--------------------------------------------------------------
var EMPTY = ee.Image.constant(0).rename('NBR').toFloat().updateMask(0);
function nbrPercentile(col, pct) {
  var nbr = col.map(function (i) {
    return ee.Image(i).normalizedDifference(['NIR', 'SWIR2']).rename('NBR').toFloat();
  });
  return ee.Image(ee.Algorithms.If(col.size().gt(0),
           nbr.reduce(ee.Reducer.percentile([pct])), EMPTY));
}
function win(y, w) { return ee.Date.fromYMD(y + w[0], w[1], w[2]); }
function preCol(y)  { return LS.filterDate(win(y, PRE_START),  win(y, PRE_END)); }
function postCol(y) { return LS.filterDate(win(y, POST_START), win(y, POST_END)); }

function yearFeatures(y) {
  var pc = preCol(y), qc = postCol(y);
  var pre  = nbrPercentile(pc, PRE_PCT).rename('NBR_pre');
  var post = nbrPercentile(qc, POST_PCT).rename('NBR_post');
  var dnbr = pre.subtract(post).rename('dNBR');           // positive = cleared / burned
  var rbr  = dnbr.divide(pre.add(1.001)).rename('RBR');   // Parks et al. 2014
  var nPre  = pc.select('NIR').count().unmask(0).rename('n_pre');
  var nPost = qc.select('NIR').count().unmask(0).rename('n_post');
  return ee.Image.cat([pre, post, dnbr, rbr, nPre, nPost, TERRAIN]).toFloat();
}

Object.keys(CLASS_YEAR).forEach(function (cid) {
  var y = CLASS_YEAR[cid];
  print('CHECK 3 – Year ' + y + ': before (previous-year Mar–May) scenes / after-burn scenes',
        preCol(y).size(), postCol(y).size());
});

//--------------------------------------------------------------
// 4. VISUAL CHECK (2025): jhum polygons should sit on red dNBR
//--------------------------------------------------------------
var VIS_YEAR = 2025;
Map.addLayer(yearFeatures(VIS_YEAR).select('dNBR').clip(roi),
             {min: -0.1, max: 0.5, palette: ['#1a9850', '#ffffbf', '#d73027']},
             'dNBR ' + VIS_YEAR + ' (red = cleared/burned)');
Map.addLayer(ee.Image().byte().paint(ref.filter(ee.Filter.eq('class_id', 7)), 1, 2),
             {palette: ['#00ffff']}, 'Jhum polygons ' + VIS_YEAR);

//--------------------------------------------------------------
// 5. SAMPLES: jhum pixels (inside polygons) + background pixels
//--------------------------------------------------------------
function samplesForYear(cid) {
  var y = CLASS_YEAR[cid];
  var feat = yearFeatures(y).updateMask(WATER_FREE);

  var polys = ref.filter(ee.Filter.eq('class_id', Number(cid)))
    .map(function (f) { return f.set('poly_id', ee.String('Y' + y + '_').cat(f.id())); });

  // Jhum = 1 : every 30 m pixel inside that year's polygons
  var pos = feat.sampleRegions({
      collection: polys, properties: ['poly_id', 'class_id'],
      scale: 30, tileScale: 8, geometries: true   // point geometry needed for asset export
    })
    .map(function (f) { return f.set({year: y, jhum: 1, sample_type: 'jhum'}); });

  // Jhum = 0 : random pixels in Mizoram, away from that year's polygons
  var exclude = ee.Image.constant(1).paint(
      polys.map(function (f) { return f.buffer(BG_BUFFER_M); }), 0);

  var bg = feat.updateMask(exclude).sample({
      region: roi.geometry(), scale: 30, numPixels: BG_PER_YEAR * 2,
      seed: SEED + y, tileScale: 8, geometries: true   // point geometry needed for asset export
    })
    .limit(BG_PER_YEAR)
    .map(function (f) {
      return f.set({year: y, jhum: 0, sample_type: 'background', class_id: 0,
                    poly_id: ee.String('BG' + y + '_').cat(f.id())});
    });

  return pos.merge(bg);
}

var allSamples = ee.FeatureCollection([]);
Object.keys(CLASS_YEAR).forEach(function (cid) {
  allSamples = allSamples.merge(samplesForYear(cid));
});

var COLUMNS = ['year', 'jhum', 'sample_type', 'class_id', 'poly_id',
               'NBR_pre', 'NBR_post', 'dNBR', 'RBR',
               'n_pre', 'n_post', 'elevation', 'slope', 'aspect'];

//--------------------------------------------------------------
// 6. EXPORTS
//--------------------------------------------------------------
// (1) Asset — needed by Part C (GEE can only read assets, not Drive)
Export.table.toAsset({
  collection: allSamples.select(COLUMNS),
  description: SAMPLE_NAME + '_ASSET',
  assetId: OUT_FOLDER + SAMPLE_NAME
});

// (2) CSV in Drive — used later by the Colab notebook
Export.table.toDrive({
  collection: allSamples,
  description: SAMPLE_NAME + '_CSV',
  fileNamePrefix: SAMPLE_NAME,
  folder: DRIVE_DIR,
  fileFormat: 'CSV',
  selectors: COLUMNS
});

//==============================================================
// PART B — APPLY THE EQUATION (run later, after Part C / Colab)
// Paste the EQ block from Part C or Colab, set RUN_MAPPING = true.
//==============================================================
var RUN_MAPPING = true;    // true = make the yearly jhum maps

// FINAL v3 EQUATION (Colab, all years + confusers)
var EQ = {
  model: 'M3_NBRpost_dNBR',
  b0: -3.484372,
  coef: {NBR_post: -11.261907, dNBR: 9.137461},
  cut: 0.2164,
  min_dNBR: 0.30   // change guard: pixel must also have dNBR >= 0.30
};

function jhumProbability(feat) {
  var z = ee.Image.constant(EQ.b0);
  Object.keys(EQ.coef).forEach(function (band) {
    z = z.add(feat.select(band).multiply(EQ.coef[band]));
  });
  return z.multiply(-1).exp().add(1).pow(-1).rename('P_jhum');
}

if (RUN_MAPPING) {
  var areaRows = [];
  Object.keys(CLASS_YEAR).forEach(function (cid) {
    var y = CLASS_YEAR[cid];
    var feat = yearFeatures(y);
    var valid = feat.select('n_pre').gt(0).and(feat.select('n_post').gt(0));
    var P = jhumProbability(feat).updateMask(valid).updateMask(WATER_FREE);
    var jhum = P.gte(EQ.cut).and(feat.select('dNBR').gte(EQ.min_dNBR || 0)).rename('Jhum');

    Map.addLayer(jhum.selfMask().clip(roi), {palette: ['red']}, 'Jhum map ' + y, y === 2025);

    // Jhum area (km2) for this year
    var km2 = jhum.selfMask().multiply(ee.Image.pixelArea()).divide(1e6)
      .reduceRegion({reducer: ee.Reducer.sum(), geometry: roi.geometry(),
                     scale: 30, maxPixels: 1e13, tileScale: 4}).get('Jhum');
    areaRows.push(ee.Feature(null, {year: y, jhum_area_km2: km2}));

    var MAP_NAME = 'MIZ_Jhum_' + y + '_Landsat30m_' + EQ.model + '_' + VERSION;
    Export.image.toDrive({
      image: ee.Image.cat([jhum.toByte(),
                           P.multiply(100).round().toByte().rename('P_jhum_pct')]).clip(roi),
      description: MAP_NAME, fileNamePrefix: MAP_NAME, folder: MAP_DIR,
      region: roi.geometry(), scale: 30, crs: 'EPSG:32646', maxPixels: 1e13,
      fileFormat: 'GeoTIFF', formatOptions: {cloudOptimized: true}
    });
  });

  var areaTable = ee.FeatureCollection(areaRows);
  print('JHUM AREA per year (km2) – if this times out, use the export task',
        ui.Chart.feature.byFeature(areaTable, 'year', ['jhum_area_km2']).setChartType('Table'));
  Export.table.toDrive({collection: areaTable,
    description: 'MIZ_Jhum_Area_1990_2025_' + EQ.model + '_' + VERSION,
    fileNamePrefix: 'MIZ_Jhum_Area_1990_2025_' + EQ.model + '_' + VERSION,
    folder: MAP_DIR, fileFormat: 'CSV', selectors: ['year', 'jhum_area_km2']});
}
