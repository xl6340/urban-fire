// ====== 1. Define Area of Interest (AOI) ======

// Load the U.S. states feature collection to get the California boundary.
var states = ee.FeatureCollection('TIGER/2018/States');
var californiaAOI = states.filter(ee.Filter.eq('NAME', 'California')).geometry();

// Center the map on California.
Map.centerObject(californiaAOI, 6);


// ====== 2. Load and Prepare Elevation Data ======

// Import the NASADEM dataset and select the elevation band.
var dataset = ee.Image('NASA/NASADEM_HGT/001');
var elevation = dataset.select('elevation');

// Clip the global elevation data to the California boundary.
var elevationCalifornia = elevation.clip(californiaAOI);


// ====== 3. Visualize the Data on the Map (Optional) ======

// Set visualization properties for the elevation.
var elevationVis = {
  min: 0,
  max: 4000, // Max elevation in California is ~4400m.
  palette: ['#006633', '#E5FFCC', '#662A00', '#D8D8D8', '#FFFFFF']
};

// Add the clipped California elevation layer to the map.
Map.addLayer(elevationCalifornia, elevationVis, 'California Elevation');


// ====== 4. Export the Image to Google Drive ======

// Set up the export task.
Export.image.toDrive({
  image: elevationCalifornia,
  description: 'NASADEM_Elevation_California', // This will be the filename.
  folder: 'GEE_Exports', // The folder in your Google Drive to save the file.
  scale: 30, // Native resolution of NASADEM is 1 arc-second (~30 meters).
  region: californiaAOI,
  maxPixels: 1e10 // Use a high maxPixels value for large exports.
});