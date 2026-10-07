// LANDFIRE fuel / vegetation export for California
// Account: xl6340@nyu.edu (run in Code Editor while logged into this account)
//
// Data availability note:
// - Vegetation / fuel layers on GEE community catalog are VERSIONED snapshots,
//   not continuous annual stacks. Currently available: LF 2023 (version 2.4.0).
// - Official Earth Engine catalog also has older EVT/EVC/EVH: LANDFIRE/.../v1_4_0
//   (~LF 2016 Remap era).
// - Annual Disturbance IS yearly: 1999–2023 in
//   projects/sat-io/open-datasets/landfire/ANNUAL-DIST
//
// Best layers for shrub/grass vs forest (Line 229):
//   EVT  - Existing Vegetation Type (lifeform / EVT groups)
//   FVT  - Fuel Vegetation Type
//   FBFM40 - Scott & Burgan fuel models (GR*/GS*/SH* = grass/shrub)
// Note: CBD/CC/CH are forest CANOPY products, not statewide herb biomass.

// ====== 1. AOI ======
var states = ee.FeatureCollection('TIGER/2018/States');
var californiaAOI = states.filter(ee.Filter.eq('NAME', 'California')).geometry();
Map.centerObject(californiaAOI, 6);

// ====== 2. LF 2023 (v2.4.0) CONUS mosaics ======
// Naming: *_LC_240  -> LC = CONUS, 240 = version 2.4.0 (LF 2023)
var evt = ee.Image('projects/sat-io/open-datasets/landfire/VEGETATION/EVT/EVT_LC_240')
  .select('EVT').rename('EVT');
var fvt = ee.Image('projects/sat-io/open-datasets/landfire/FUEL/FVT/FVT_LC_240')
  .select('FVT').rename('FVT');
var fbfm40 = ee.Image('projects/sat-io/open-datasets/landfire/FUEL/FBFM40/F40_LC_240')
  .select('F40').rename('FBFM40');

// Optional forest canopy (NOT herbaceous biomass)
var cbd = ee.Image('projects/sat-io/open-datasets/landfire/FUEL/CBD/CBD_LC_240')
  .select('CBD').rename('CBD');
var cc = ee.Image('projects/sat-io/open-datasets/landfire/FUEL/CC/CC_LC_240')
  .select('CC').rename('CC');

// Optional older official EVT (~LF 2016 Remap / v1.4.0)
var evt_v14 = ee.ImageCollection('LANDFIRE/Vegetation/EVT/v1_4_0')
  .filter(ee.Filter.eq('system:index', 'CONUS'))
  .first()
  .select('EVT')
  .rename('EVT');

// Clip to CA
var evt_ca = evt.clip(californiaAOI);
var fvt_ca = fvt.clip(californiaAOI);
var f40_ca = fbfm40.clip(californiaAOI);
var cbd_ca = cbd.clip(californiaAOI);
var cc_ca = cc.clip(californiaAOI);
var evt14_ca = evt_v14.clip(californiaAOI);

// Quick map preview (EVT)
Map.addLayer(evt_ca, {min: 3001, max: 7999, palette: ['a6d96a', '1a9850', 'fee08b', 'd73027']}, 'EVT LF2023', false);
Map.addLayer(f40_ca, {min: 91, max: 204}, 'FBFM40 LF2023', true);

// ====== 3. Export helpers ======
function exportLF(image, description, folder) {
  Export.image.toDrive({
    image: image.toInt16(),
    description: description,
    folder: folder,
    scale: 30,
    region: californiaAOI,
    crs: 'EPSG:3310',
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF',
    formatOptions: {cloudOptimized: true}
  });
}

// Drive folder for xl6340@nyu.edu exports
var outFolder = 'GEE_LANDFIRE_CA';

exportLF(evt_ca, 'LANDFIRE_LF2023_EVT_CA', outFolder);
exportLF(fvt_ca, 'LANDFIRE_LF2023_FVT_CA', outFolder);
exportLF(f40_ca, 'LANDFIRE_LF2023_FBFM40_CA', outFolder);
// Uncomment if needed:
// exportLF(cbd_ca, 'LANDFIRE_LF2023_CBD_CA', outFolder);
// exportLF(cc_ca, 'LANDFIRE_LF2023_CC_CA', outFolder);
// exportLF(evt14_ca, 'LANDFIRE_v14_EVT_CA', outFolder);

print('Queued Drive exports for EVT / FVT / FBFM40 (LF 2023, CONUS→CA).');
print('Run Tasks tab while signed in as xl6340@nyu.edu.');
