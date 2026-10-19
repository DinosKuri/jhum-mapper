//==============================================================
// MIZORAM JHUM DETECTION — PART C: THRESHOLD TEST
//
// INPUT : the SAMPLE TABLE made by Part A (NOT the polygons)
// ANSWERS:
//   1. How different are jhum and background pixels each year?
//   2. Which threshold separates them best each year?
//   3. Does ONE threshold learned from 6 years work on the 7th?
//      (leave-one-year-out = proof of transferability)
//==============================================================

//--------------------------------------------------------------
// 0. SETTINGS — edit only this block
//--------------------------------------------------------------
var SAMPLE_ASSET = 'projects/phdmrsei89448/assets/MIZ_Jhum_Samples_Landsat_1990_2025_v3';   // PRIVATE ASSET — see README, Private Earth Engine assets
// Which variable to test, and which side of the threshold is jhum:
//   'dNBR'     , 'above'  -> jhum if dNBR >= threshold   (change-based)
//   'RBR'      , 'above'  -> jhum if RBR  >= threshold   (change-based)
//   'NBR_post' , 'below'  -> jhum if NBR_post <= threshold (post-burn state)
var VAR       = 'NBR_post';
var DIRECTION = 'below';
var DRIVE_DIR = 'Jhum_tests_v3';   // GEE uses a folder NAME, not a path
var VERSION   = 'v3';
var YEARS     = [1990, 2000, 2005, 2010, 2015, 2020, 2025];
var TH = (DIRECTION === 'above')
         ? ee.List.sequence(0, 0.6, 0.02)       // for dNBR / RBR
         : ee.List.sequence(-0.4, 0.6, 0.02);   // for NBR_post

// Rule used everywhere below: which pixels are called 'jhum' at threshold t
function jhumRule(t) {
  return (DIRECTION === 'above') ? ee.Filter.gte(VAR, t) : ee.Filter.lte(VAR, t);
}

//--------------------------------------------------------------
// 1. SAFETY CHECKS — is this really the sample table?
//--------------------------------------------------------------
var RAW = ee.FeatureCollection(SAMPLE_ASSET);
var NEED = ['year', 'jhum', 'dNBR', 'RBR', 'NBR_pre', 'NBR_post', 'n_pre', 'n_post'];
var isSampleTable = RAW.first().propertyNames().containsAll(NEED);

print('CHECK 1 – Is this the sample table?',
      ee.Algorithms.If(isSampleTable,
        'YES – sample table found.',
        'NO – this is not the sample table (maybe the polygon file). ' +
        'Set SAMPLE_ASSET to the table made by Part A.'));
print('CHECK 2 – Total rows', RAW.size());
print('CHECK 3 – Rows per year', RAW.aggregate_histogram('year'));

var S = RAW.filter(ee.Filter.gt('n_pre', 0))
           .filter(ee.Filter.gt('n_post', 0))
           .filter(ee.Filter.notNull([VAR]));
print('CHECK 4 – Rows with valid before/after images', S.size());

//--------------------------------------------------------------
// 2. SUMMARY PER YEAR (medians of jhum vs background)
//--------------------------------------------------------------
function medians(c) {
  return ee.List(ee.Dictionary(c.reduceColumns(ee.Reducer.median().repeat(4),
           ['dNBR', 'NBR_pre', 'RBR', 'NBR_post'])).get('median'));
}

var summary = ee.FeatureCollection(YEARS.map(function (y) {
  var yr = S.filter(ee.Filter.eq('year', y));
  var j  = yr.filter(ee.Filter.eq('jhum', 1));
  var b  = yr.filter(ee.Filter.eq('jhum', 0));
  var mj = medians(j), mb = medians(b);
  return ee.Feature(null, {
    year_label: 'Y' + y, year: y,
    n_jhum: j.size(), n_background: b.size(),
    median_clear_obs_post: yr.aggregate_array('n_post').reduce(ee.Reducer.median()),
    jhum_dNBR: mj.get(0), jhum_NBRpre: mj.get(1), jhum_RBR: mj.get(2), jhum_NBR_post: mj.get(3),
    bg_dNBR: mb.get(0),   bg_NBRpre: mb.get(1),   bg_RBR: mb.get(2),   bg_NBR_post: mb.get(3)
  });
}));

print('OUTPUT 1 – Summary per year', ui.Chart.feature.byFeature(summary, 'year_label',
      ['n_jhum', 'n_background', 'median_clear_obs_post',
       'jhum_dNBR', 'bg_dNBR', 'jhum_RBR', 'bg_RBR', 'jhum_NBRpre', 'bg_NBRpre',
       'jhum_NBR_post', 'bg_NBR_post'])
      .setChartType('Table'));

print(ui.Chart.feature.byFeature(summary, 'year_label', ['jhum_' + VAR.replace('NBR_pre','NBRpre'), 'bg_' + VAR.replace('NBR_pre','NBRpre')])
      .setChartType('ColumnChart')
      .setOptions({title: 'Median ' + VAR + ': jhum vs background',
                   colors: ['#d73027', '#1a9850'], vAxis: {title: VAR}}));

//--------------------------------------------------------------
// 3. ACCURACY AT EVERY THRESHOLD
//--------------------------------------------------------------
function countsFor(coll, label) {
  var pos = coll.filter(ee.Filter.eq('jhum', 1));
  var neg = coll.filter(ee.Filter.eq('jhum', 0));
  var P = pos.size(), N = neg.size();
  return ee.FeatureCollection(TH.map(function (t) {
    t = ee.Number(t).multiply(1000).round().divide(1000);   // clean values (fixes chart glitch at 0)
    return ee.Feature(null, {
      group: label, threshold: t, P: P, N: N,
      TP: pos.filter(jhumRule(t)).size(),
      FP: neg.filter(jhumRule(t)).size()
    });
  }));
}

function addMetrics(f) {
  var tp = ee.Number(f.get('TP')), fp = ee.Number(f.get('FP'));
  var P  = ee.Number(f.get('P')),  N  = ee.Number(f.get('N'));
  var tn = N.subtract(fp);
  var pa = tp.divide(P.max(1));                    // producer's accuracy: jhum found
  var ua = tp.divide(tp.add(fp).max(1));           // user's accuracy: jhum correct
  var f1 = pa.multiply(ua).multiply(2).divide(pa.add(ua).max(1e-9));
  var oa = tp.add(tn).divide(P.add(N).max(1));
  var youden = pa.subtract(fp.divide(N.max(1)));   // choose threshold by max Youden J
  return f.set({PA_jhum: pa, UA_jhum: ua, F1: f1, OA: oa, Youden: youden});
}

var perYear = ee.FeatureCollection([]);
YEARS.forEach(function (y) {
  perYear = perYear.merge(countsFor(S.filter(ee.Filter.eq('year', y)), 'Y' + y));
});
perYear = perYear.merge(countsFor(S, 'ALL')).map(addMetrics);

var groups = YEARS.map(function (y) { return 'Y' + y; }).concat(['ALL']);
var best = ee.FeatureCollection(groups.map(function (g) {
  return ee.Feature(perYear.filter(ee.Filter.eq('group', g))
                           .sort('Youden', false).first());
}));

print('OUTPUT 2 – Best ' + VAR + ' threshold per year',
      ui.Chart.feature.byFeature(best, 'group',
        ['threshold', 'Youden', 'OA', 'PA_jhum', 'UA_jhum', 'F1', 'P', 'N'])
      .setChartType('Table'));

print(ui.Chart.feature.groups(perYear.filter(ee.Filter.neq('group', 'ALL')),
        'threshold', 'Youden', 'group')
      .setOptions({title: 'Accuracy (Youden J) vs ' + VAR + ' threshold — one line per year',
                   hAxis: {title: VAR + ' threshold'}, vAxis: {title: 'Youden J'},
                   lineWidth: 2, pointSize: 3}));

//--------------------------------------------------------------
// 4. LEAVE-ONE-YEAR-OUT: learn on 6 years, test on the 7th
//--------------------------------------------------------------
var loyo = ee.FeatureCollection(YEARS.map(function (y) {
  var train = S.filter(ee.Filter.neq('year', y));
  var test  = S.filter(ee.Filter.eq('year', y));
  var tBest = ee.Number(countsFor(train, 'train').map(addMetrics)
                 .sort('Youden', false).first().get('threshold'));
  var pos = test.filter(ee.Filter.eq('jhum', 1));
  var neg = test.filter(ee.Filter.eq('jhum', 0));
  return addMetrics(ee.Feature(null, {
    test_year: 'Y' + y, threshold_from_other_years: tBest,
    P: pos.size(), N: neg.size(),
    TP: pos.filter(jhumRule(tBest)).size(),
    FP: neg.filter(jhumRule(tBest)).size()
  }));
}));

print('OUTPUT 3 – Leave-one-year-out test (' + VAR + ')',
      ui.Chart.feature.byFeature(loyo, 'test_year',
        ['threshold_from_other_years', 'OA', 'PA_jhum', 'UA_jhum', 'F1', 'Youden'])
      .setChartType('Table'));

//--------------------------------------------------------------
// 5. EXPORT TABLES TO DRIVE (use these if the Console times out)
//--------------------------------------------------------------
var tag = VAR + '_1990_2025_' + VERSION;
Export.table.toDrive({collection: summary, description: 'MIZ_Jhum_C1_Summary_' + tag,
  fileNamePrefix: 'MIZ_Jhum_C1_Summary_' + tag, folder: DRIVE_DIR, fileFormat: 'CSV'});
Export.table.toDrive({collection: perYear, description: 'MIZ_Jhum_C2_ThresholdCurves_' + tag,
  fileNamePrefix: 'MIZ_Jhum_C2_ThresholdCurves_' + tag, folder: DRIVE_DIR, fileFormat: 'CSV'});
Export.table.toDrive({collection: loyo, description: 'MIZ_Jhum_C3_LeaveOneYearOut_' + tag,
  fileNamePrefix: 'MIZ_Jhum_C3_LeaveOneYearOut_' + tag, folder: DRIVE_DIR, fileFormat: 'CSV'});
