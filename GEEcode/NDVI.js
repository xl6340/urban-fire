// ====== 1. Define Area of Interest (AOI) ======

// Load the U.S. states feature collection.
var states = ee.FeatureCollection('TIGER/2018/States');

// Filter the collection to get the geometry for California.
var californiaAOI = states.filter(ee.Filter.eq('NAME', 'California')).geometry();
Map.centerObject(californiaAOI, 6); // Center the map on California.


// ====== 2. Define QA Masking Function ======

// This function checks the 'SummaryQA' band to keep only "Good Data" pixels.
// This works for both MOD13A1 and MYD13A1.
function maskGoodQA(image) {
  var qa = image.select('SummaryQA');
  // The 'eq(0)' part corresponds to the "Good Data" rank key.
  var goodDataMask = qa.eq(0);
  return image.updateMask(goodDataMask);
}


// ====== 3. Batch Export for Years 2000-2024 (MODIFIED) ======

// Create a list of years from 2001 to 2024.
var years = ee.List.sequence(2001, 2024);

years.evaluate(function(yearList) {
  yearList.forEach(function(year) {
    var startDate = ee.Date.fromYMD(year, 1, 1);
    var endDate = ee.Date.fromYMD(year, 12, 31);

    // --- NEW: Load Terra (MOD) collection ---
    var modisTerra = ee.ImageCollection('MODIS/061/MOD13A1')
                        .filterBounds(californiaAOI)
                        .filterDate(startDate, endDate)
                        .map(maskGoodQA); // Apply the QA mask

    // --- NEW: Load Aqua (MYD) collection ---
    var modisAqua = ee.ImageCollection('MODIS/061/MYD13A1')
                       .filterBounds(californiaAOI)
                       .filterDate(startDate, endDate)
                       .map(maskGoodQA); // Apply the QA mask
                       
    // --- NEW: Merge the two collections ---
    var modisCombined = modisTerra.merge(modisAqua);

    // Create a max composite from the COMBINED collection.
    var annualComposite = modisCombined.select('NDVI').max();

    // **MODIFICATION 1: Remove values lower than 0.**
    // Create a mask to keep only pixels with values >= 0.
    var positiveMask = annualComposite.gte(0);

    // **MODIFICATION 2: Apply scale factor.**
    // Apply the mask and the scale factor to get NDVI from 0-1.
    var annualNdviScaled = annualComposite
                            .updateMask(positiveMask) // Apply mask
                            .multiply(0.0001)       // Apply scale factor
                            .rename('NDVI');        // Rename the band for clarity

    // Clip the final, scaled image to the California boundary.
    var annualNdviCalifornia = annualNdviScaled.clip(californiaAOI);

    // **MODIFICATION 3: Export data from 2000 to 2024.**
    // Export the processed image to your Google Drive.
    Export.image.toDrive({
      image: annualNdviCalifornia,
      description: 'ndviMax' + year, // Filename for the export.
      folder: 'GEE_ndviMax',         // Folder in your Google Drive.
      scale: 500,                      // Resolution in meters.
      region: californiaAOI,
      crs: 'EPSG:3310',      // <-- THE FIX: Use California Albers projection
      maxPixels: 1e10,                 // Handle large export size.
      fileFormat: 'GeoTIFF',
      formatOptions: {
        cloudOptimized: true // More efficient format
      }
    });
  });
});


// ====== 4. Optional Visualization (MODIFIED) ======
// Use this section to preview the final year's data on the map.
var lastYear = 2024;
var visStartDate = String(lastYear) + '-01-01';
var visEndDate = String(lastYear) + '-12-31';

// --- NEW: Load and mask Terra (MOD) for visualization ---
var visTerra = ee.ImageCollection('MODIS/061/MOD13A1')
                    .filterBounds(californiaAOI)
                    .filterDate(visStartDate, visEndDate)
                    .map(maskGoodQA);

// --- NEW: Load and mask Aqua (MYD) for visualization ---
var visAqua = ee.ImageCollection('MODIS/061/MYD13A1')
                   .filterBounds(californiaAOI)
                   .filterDate(visStartDate, visEndDate)
                   .map(maskGoodQA);

// --- NEW: Merge and get max composite ---
var visCombined = visTerra.merge(visAqua);
var visComposite = visCombined.select('NDVI').max();

// --- Continue with scaling and clipping as before ---
var visPositiveMask = visComposite.gte(0);
var visScaled = visComposite.updateMask(visPositiveMask).multiply(0.0001);
var visClipped = visScaled.clip(californiaAOI);

// Define visualization parameters for scaled 0-1 NDVI.
var ndviVisScaled = {
  min: 0.0,
  max: 1.0,
  palette: [
    'FFFFFF', 'CE7E45', 'DF923D', 'F1B555', 'FCD163', '99B718', '74A901',
    '66A000', '529400', '3E8601', '207401', '056201', '004C00', '023B01'
  ],
};

Map.addLayer(visClipped, ndviVisScaled, 'Scaled NDVI (Combined) ' + lastYear);


// ====== 5. Check NDVI Distribution (NEW CODE) ======
// This section will automatically use the updated 'visClipped' image
var ndviHistogram = ui.Chart.image.histogram({
  image: visClipped,       // Use the final processed image for visualization.
  region: californiaAOI,
  scale: 500,              // Use a consistent scale for a representative sample.
  maxBuckets: 64           // Adjust the number of bins in the histogram.
}).setOptions({
  title: 'NDVI Distribution in California (Combined) (' + lastYear + ')',
  hAxis: { title: 'NDVI Value' },
  vAxis: { title: 'Pixel Count' }
});

// Print the histogram to the console.
print(ndviHistogram);