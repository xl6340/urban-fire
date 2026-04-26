// ======== SCRIPT TO EXPORT ORIGINAL ANNUAL LAND COVER FOR CALIFORNIA ========

// --- 1. Define Region of Interest (ROI) ---
// Load US States boundaries and filter for California
var states = ee.FeatureCollection('TIGER/2018/States');
var california = states.filter(ee.Filter.eq('NAME', 'California')).geometry();

// Center the map on the ROI (optional)
Map.centerObject(california, 6);
Map.addLayer(california, {color: 'FF0000'}, 'California Border', false);

// --- 2. Define Time Range and Data ---
// (Section 2 for Reclassification was removed)
var startYear = 2001;
var endYear = 2024; // Will export all available years up to this one

// Load the MODIS Land Cover collection, select the original IGBP band
var collection = ee.ImageCollection('MODIS/061/MCD12Q1')
                   .filterDate(ee.Date.fromYMD(startYear, 1, 1), ee.Date.fromYMD(endYear, 12, 31))
                   .select('LC_Type1');

// --- 3. Prepare for Export Loop ---
// (Section 4 for Remapping was removed)

// Get a list of all images in the original collection
var imageList = collection.toList(collection.size());

// Get a list of the years available.
var systemIndices = collection.aggregate_array('system:index').getInfo();
var years = systemIndices.map(function(index) {
  return index.substring(0, 4); // Get year from index (e.g., '2001_01_01')
});

// --- 4. Create Export Tasks ---

// Loop through the list of years and create an export task for each one.
// You must run these tasks from the 'Tasks' tab in the Code Editor.
for (var i = 0; i < years.length; i++) {
  var year = years[i];
  var image = ee.Image(imageList.get(i));
  var fileName = 'lc' + year; // Changed filename slightly

  Export.image.toDrive({
    // image is already just the 'LC_Type1' band
    image: image.clip(california).toUint8(), 
    description: fileName,
    fileNamePrefix: fileName,
    folder: 'GEE_lc_original', // Changed folder
    region: california,    // The geometry to export
    scale: 500,            // 500 meters in the new projection
    crs: 'EPSG:3310',      // Use California Albers projection
    maxPixels: 1e13,
    formatOptions: {
      cloudOptimized: true
    }
  });
}

// --- 5. Optional: Add Most Recent Layer to Map for Verification ---
// --- MODIFIED: Visualization for the 17 IGBP classes ---
var igbpVis = {
  min: 1,
  max: 17,
  palette: [
    '05450a', // 1  Evergreen Needleleaf Forests
    '086a10', // 2  Evergreen Broadleaf Forests
    '54a708', // 3  Deciduous Needleleaf Forests
    '78d203', // 4  Deciduous Broadleaf Forests
    '009900', // 5  Mixed Forests
    'c6b044', // 6  Closed Shrublands
    'dcd159', // 7  Open Shrublands
    'dade48', // 8  Woody Savannas
    'fbff13', // 9  Savannas
    'b6ff05', // 10 Grasslands
    '27ff87', // 11 Permanent Wetlands
    'c24f44', // 12 Croplands
    'a5a5a5', // 13 Urban and Built-up Lands
    'ff6d4c', // 14 Cropland/Natural Vegetation Mosaics
    '69fff8', // 15 Snow and Ice
    'f9ffa4', // 16 Barren
    '1c0dff'  // 17 Water Bodies
  ],
};

// Get the most recent image from the ORIGINAL collection
var mostRecent = collection.sort('system:time_start', false).first();
Map.addLayer(
  mostRecent.clip(california).toUint8(), 
  igbpVis, // Use the new 17-class visualization
  'Original Land Cover (' + years[years.length - 1] + ')'
);

print('Script finished. Check the "Tasks" tab to run your exports.', years);