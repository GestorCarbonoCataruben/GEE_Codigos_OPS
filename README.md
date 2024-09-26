# GEE_Codigos_OPS SCRIPT #1 IDENTIFICACIÓN DE TURBERAS ANUAL 

Map.addLayer({
  eeObject: Datos, 
  visParams: {min: 0, max: 1},
  name: 'Entrenamiento Analisis', 
  shown: 0, 
});

// Localizacion
var area = ee.FeatureCollection(roi);


// Establezca el área de estudio como centro del mapa.
Map.centerObject({object: area, zoom: 12});


//---------------------------------------------- OPTRAM DATOS SAR ----------------------------------------

var time_start = ee.Date('2017');
var time_end = ee.Date('2024');
var time_dif = time_end.difference(time_start, 'day');

var sen1 = ee.ImageCollection("COPERNICUS/S1_GRD")
.filterDate(time_start, time_end)
.filterBounds(area)
.filter(ee.Filter.listContains('transmitterReceiverPolarisation','VV'))
.filter(ee.Filter.eq('instrumentMode','IW')).select('VV');

//print(sen1.aggregate_array('orbitProperties_pass').distinct());

var asc = sen1
.filter(ee.Filter.eq('orbitProperties_pass','ASCENDING'));

var des = sen1
.filter(ee.Filter.eq('orbitProperties_pass','DESCENDING'));


var list_dates = ee.List.sequence(0, time_dif, 10).map(function(interval){
  return ee.Date.fromYMD(2017, 1, 1).advance(interval, 'days');
  });

var asc_10days = ee.ImageCollection(list_dates.map(function(dates){
  var start_date  = ee.Date(dates)
  var end_date = start_date.advance(10, 'days')
  var composite = asc.filterDate(start_date, end_date).mean()
  var bands = composite.bandNames().size()
  return composite
  .set('system:time_start', start_date.millis())
  .set('system:time_end', end_date.millis())
  .set('band_number', bands)
  })).filter(ee.Filter.eq('band_number',1));
  
var des_10days = ee.ImageCollection(list_dates.map(function(dates){
  var start_date  = ee.Date(dates)
  var end_date = start_date.advance(10, 'days')
  var composite = des.filterDate(start_date, end_date).mean()
  var bands = composite.bandNames().size()
  return composite
  .set('system:time_start', start_date.millis())
  .set('system:time_end', end_date.millis())
  .set('band_number', bands)
  })).filter(ee.Filter.eq('band_number',1));
  

var asc_sigma = asc_10days.map(function(img){
  var sigma = ee.Image(10).pow(img.divide(10)).rename('sigma')
  var speckel = sigma.focalMean(30, 'square', 'meters');
  return speckel
  .copyProperties(img, img.propertyNames())
  });
  
var des_sigma = des_10days.map(function(img){
  var sigma = ee.Image(10).pow(img.divide(10)).rename('sigma')
  var speckel = sigma.focalMean(30, 'square', 'meters');
  return speckel
  .copyProperties(img, img.propertyNames())
  });
  
  
var water_mask = ee.ImageCollection("GOOGLE/DYNAMICWORLD/V1").select('label')
.filterDate(time_start, time_end)
.filterBounds(area).mode().eq(0).not()

var urban_mask = ee.ImageCollection("GOOGLE/DYNAMICWORLD/V1").select('label')
.filterDate(time_start, time_end)
.filterBounds(area).mode().eq(6).not()


var asc_min = asc_sigma.min()
var asc_max = asc_sigma.max()

var asc_sm = asc_sigma.map(function(img){
  var index = (img.subtract(asc_min)).divide(asc_max.subtract(asc_min))
  var date = img.date().format('YYYY-MM-dd')
  return index.multiply(water_mask).multiply(urban_mask).rename('sm_asc')
  .copyProperties(img,['system:time_start'])
  .set('date', ee.String(date))
  })
  
Map.addLayer(asc_sm.filterDate('2020','2021').toBands().clip(area),[],'sm_asc',false)

var des_min = des_sigma.min()
var des_max = des_sigma.max()

var des_sm = des_sigma.map(function(img){
  var index = (img.subtract(des_min)).divide(des_max.subtract(des_min))
  var date = img.date().format('YYYY-MM-dd')
  return index.multiply(water_mask).multiply(urban_mask).rename('sm_des')
  .copyProperties(img,['system:time_start'])
  .set('date', ee.String(date))
  })
  
Map.addLayer(des_sm.filterDate('2020','2021').toBands().clip(roi),[],'sm_des',false)

//----------------------- DATOS PROMEDIO SOIL MESUREMENT SAR ANUAL ------------------------------
var asc_sm_mean = asc_sm.mean().clip(roi);
var des_sm_mean = des_sm.mean().clip(roi);



//------------------------ OPTRAM DATOS OPTICOS ------------------------------------------------

var sentinel = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
.filterBounds(roi)
.filterDate('2020','2021')
.filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',65));

function temporal_collection(collection, start, count, interval, unit){
  var seq = ee.List.sequence(0, ee.Number(count).subtract(1));
  var origin_date = ee.Date(start);
  return ee.ImageCollection(seq.map(function(i){
    var start_date = origin_date.advance(ee.Number(interval).multiply(i), unit);
    var end_date = origin_date.advance(ee.Number(interval).multiply(ee.Number(i).add(1)), unit);
    return collection.filterDate(start_date, end_date).median()
    .set('system:time_start', start_date.millis())
    .set('system:time_end', end_date.millis())
    }))
  }

var monthly = temporal_collection(sentinel, '2020', 12, 1, 'month');

var parameters = monthly.map(function(img){
  var bands = img.select('B.*').multiply(0.0001);
  var ndvi = bands.normalizedDifference(['B8','B4']).rename('ndvi');
  var str = bands.expression('((1 - swir) ** 2)/(2 * swir)',{'swir': bands.select('B12')}).rename('str');
  var stack = ee.Image.cat([ndvi, str]);
  return stack
  .copyProperties(img, img.propertyNames())
  });

var str_full_cover = parameters.map(function(img){
  var thr_full = img.select('ndvi').gt(0.3);
  var str_full = img.select('str').updateMask(thr_full);
  return str_full
  .copyProperties(img, img.propertyNames())
  });
  
var str_bareland = parameters.map(function(img){
  var thr_bare = img.select('ndvi').gte(0).and(img.select('ndvi').lt(0.2));
  var str_bare = img.select('str').updateMask(thr_bare)
  return str_bare
  .copyProperties(img, img.propertyNames())
  })
  
var scale = 30 
var vw = ee.Number(str_full_cover.max().reduceRegion({
  reducer: ee.Reducer.max(), geometry: roi, scale: scale
  }).values().get(0));
  
var vd = ee.Number(str_full_cover.min().reduceRegion({
  reducer: ee.Reducer.min(), geometry: roi, scale: scale
  }).values().get(0));
  

var iw = ee.Number(str_bareland.max().reduceRegion({
  reducer: ee.Reducer.max(), geometry: roi, scale: scale
  }).values().get(0)); 
  

var id = ee.Number(str_bareland.min().reduceRegion({
  reducer: ee.Reducer.min(), geometry: roi, scale: scale
  }).values().get(0)); 

var sw = vw.subtract(iw);
var sd = vd.subtract(id);

var sm = parameters.map(function(img){
  var mask = img.select('ndvi').lt(-0.1).not()
  var index = img.expression('(id + sd * ndvi - str) / (id - iw + (sd - sw) * ndvi)',
  {'id': id, 'sd': sd, 'ndvi': img.select('ndvi'), 'str': img.select('str'), 'iw': iw, 'sw': sw})
  .rename('soil_moisture');
  return index.multiply(1000).updateMask(mask)
  .copyProperties(img, img.propertyNames())
  })
  
print(
  ui.Chart.image.series(sm, muestra, ee.Reducer.first(), scale, 'system:time_start')
  )
  
Map.addLayer(sm.toBands().clip(roi),[],'sm',false);

var sm_mean = sm.first().clip(roi);
Map.addLayer(sm_mean,[],'sm',false);

//------------------------------------CLASIFICACIÓN SUPERVISADA----------------------------------
var peatland = ee.Image(Datos).select('b1').clip(area);
var clases = peatland.sample({
  region: area,
  scale: 30,
  numPixels: 500, 
  seed: 0, 
  geometries: true
});

// Divide el conjunto de entrenamiento
var muestraConRandom = clases.randomColumn();                  // Agregar una columna 'random' para la división
var datosEntrenamiento = muestraConRandom.filter(ee.Filter.lte('random', 0.7));   // ~70% datos para entrenamiento
var datosValidacion = muestraConRandom.filter(ee.Filter.gt('random', 0.7));       // ~30% datos para validación

//Define both optical and SAR to train your data
var opt_sar = ee.Image.cat(asc_sm_mean, des_sm_mean, sm_mean);
var bands_opt_sar = ['sm_asc', 'sm_des', 'soil_moisture']
var training_opt_sar = opt_sar.select(bands_opt_sar).sampleRegions({
  collection: datosEntrenamiento,
  properties: ['b1'],
  scale: 30 });
  
//Train the classifier
var classifier = ee.Classifier.smileRandomForest(50).train({
  features: training_opt_sar, 
  classProperty: 'b1',
  inputProperties: bands_opt_sar 
});


// Calculate variable importance
var importance = ee.Dictionary(classifier.explain().get('importance'));

// Calculate relative importance
var sum = importance.values().reduce(ee.Reducer.sum());

var relativeImportance = importance.map(function(key, val) {
   return (ee.Number(val).multiply(100)).divide(sum);
  });
print(relativeImportance);

// Create a FeatureCollection so we can chart it
var importanceFc = ee.FeatureCollection([
  ee.Feature(null, relativeImportance)
]);

var chart = ui.Chart.feature.byProperty({
  features: importanceFc
}).setOptions({
      title: 'Feature Importance',
      vAxis: {title: 'Importance'},
      hAxis: {title: 'Feature'}
  });
print(chart);


var test = opt_sar.select(bands_opt_sar).sampleRegions({
  collection: datosValidacion,
  properties: ['b1'],
  scale: 30
});

// Tune the numberOfTrees parameter.
var numTreesList = ee.List.sequence(10, 150, 10);

var accuracies = numTreesList.map(function(numTrees) {
  var classifier = ee.Classifier.smileRandomForest(numTrees)
      .train({
        features: training_opt_sar,
        classProperty: 'b1',
        inputProperties: bands_opt_sar
      });

  // Here we are classifying a table instead of an image
  // Classifiers work on both images and tables
  return test
    .classify(classifier)
    .errorMatrix('b1', 'classification')
    .accuracy();
});

var chart = ui.Chart.array.values({
  array: ee.Array(accuracies),
  axis: 0,
  xLabels: numTreesList
  }).setOptions({
      title: 'Hyperparameter Tuning for the numberOfTrees Parameters',
      vAxis: {title: 'Validation Accuracy'},
      hAxis: {title: 'Number of Trees', gridlines: {count: 15}}
  });
print(chart);

// Tuning Multiple Parameters
// We can tune many parameters together using
// nested map() functions
// Let's tune 2 parameters
// numTrees and bagFraction 
var numTreesList = ee.List.sequence(10, 150, 10);
var bagFractionList = ee.List.sequence(0.1, 0.9, 0.1);

var accuracies = numTreesList.map(function(numTrees) {
  return bagFractionList.map(function(bagFraction) {
     var classifier = ee.Classifier.smileRandomForest({
       numberOfTrees: numTrees,
       bagFraction: bagFraction
     })
      .train({
        features: training_opt_sar,
        classProperty: 'b1',
        inputProperties: bands_opt_sar
      });

    // Here we are classifying a table instead of an image
    // Classifiers work on both images and tables
    var accuracy = test
      .classify(classifier)
      .errorMatrix('b1', 'classification')
      .accuracy();
    return ee.Feature(null, {'accuracy': accuracy,
      'numberOfTrees': numTrees,
      'bagFraction': bagFraction});
  });
}).flatten();
var resultFc = ee.FeatureCollection(accuracies);

// Alternatively we can automatically pick the parameters
// that result in the highest accuracy
var resultFcSorted = resultFc.sort('accuracy', false);
var highestAccuracyFeature = resultFcSorted.first();
print(highestAccuracyFeature);
var highestAccuracy = highestAccuracyFeature.getNumber('accuracy');
var optimalNumTrees = highestAccuracyFeature.getNumber('numberOfTrees');
var optimalBagFraction = highestAccuracyFeature.getNumber('bagFraction');
//var properties = resultFc.getInfo();


// Use the optimal parameters in a model and perform final classification
var optimalModel = ee.Classifier.smileRandomForest({
  numberOfTrees: optimalNumTrees,
  bagFraction: optimalBagFraction
}).train({
  features: training_opt_sar, 
  classProperty: 'b1',
  inputProperties: bands_opt_sar
});

var visClass =({
  min: 0,
  max: 1,
  palette: ['yellow', 'green']
});
var finalClassification = opt_sar.select(bands_opt_sar).classify(optimalModel);
//var bosque = finalClassification.updateMask(finalClassification.neq(0));
Map.addLayer(finalClassification, visClass, 'Clasificacion Optimizada', 1);


// Clasificar el conjunto de datos de validación usando el clasificador entrenado
var classifiedValidation = test.classify(optimalModel);

// Calcular la matriz de confusión
var errorMatrix = classifiedValidation.errorMatrix('b1', 'classification');

// Calcular el índice Kappa
var kappa = errorMatrix.kappa();

// Imprimir la matriz de confusión y el índice Kappa
print('Matriz de Confusión:', errorMatrix);
print('Índice Kappa:', kappa);

# SCRIPT ANALISIS DE TURBERAS CLASIFICACIÓN SUPERVISADA

Map.addLayer({
  eeObject: Datos, 
  visParams: {min: 0, max: 1},
  name: 'Entrenamiento Analisis', 
  shown: 0, 
});

// Localizacion
var area = ee.FeatureCollection(roi);

// Define a Geometry object.
var geometry = ee.Geometry(roi);

// Apply the area method to the Geometry object.
var geometryArea = geometry.area({'maxError': 1}).divide(10000);

// Print the result to the console.
print('geometry.area(...) =', geometryArea);

// Establezca el área de estudio como centro del mapa.
Map.centerObject({object: area, zoom: 12});

//---------------------------------------------- OPTRAM DATOS SAR ----------------------------------------

var time_start = ee.Date('2017');
var time_end = ee.Date('2024');
var time_dif = time_end.difference(time_start, 'day');

var sen1 = ee.ImageCollection("COPERNICUS/S1_GRD")
.filterDate(time_start, time_end)
.filterBounds(area)
.filter(ee.Filter.listContains('transmitterReceiverPolarisation','VV'))
.filter(ee.Filter.eq('instrumentMode','IW')).select('VV');

//print(sen1.aggregate_array('orbitProperties_pass').distinct());

var asc = sen1
.filter(ee.Filter.eq('orbitProperties_pass','ASCENDING'));

var des = sen1
.filter(ee.Filter.eq('orbitProperties_pass','DESCENDING'));


var list_dates = ee.List.sequence(0, time_dif, 10).map(function(interval){
  return ee.Date.fromYMD(2017, 1, 1).advance(interval, 'days');
  });

var asc_10days = ee.ImageCollection(list_dates.map(function(dates){
  var start_date  = ee.Date(dates)
  var end_date = start_date.advance(10, 'days')
  var composite = asc.filterDate(start_date, end_date).mean()
  var bands = composite.bandNames().size()
  return composite
  .set('system:time_start', start_date.millis())
  .set('system:time_end', end_date.millis())
  .set('band_number', bands)
  })).filter(ee.Filter.eq('band_number',1));
  
var des_10days = ee.ImageCollection(list_dates.map(function(dates){
  var start_date  = ee.Date(dates)
  var end_date = start_date.advance(10, 'days')
  var composite = des.filterDate(start_date, end_date).mean()
  var bands = composite.bandNames().size()
  return composite
  .set('system:time_start', start_date.millis())
  .set('system:time_end', end_date.millis())
  .set('band_number', bands)
  })).filter(ee.Filter.eq('band_number',1));
  

var asc_sigma = asc_10days.map(function(img){
  var sigma = ee.Image(10).pow(img.divide(10)).rename('sigma')
  var speckel = sigma.focalMean(30, 'square', 'meters');
  return speckel
  .copyProperties(img, img.propertyNames())
  });
  
var des_sigma = des_10days.map(function(img){
  var sigma = ee.Image(10).pow(img.divide(10)).rename('sigma')
  var speckel = sigma.focalMean(30, 'square', 'meters');
  return speckel
  .copyProperties(img, img.propertyNames())
  });
  
  
var water_mask = ee.ImageCollection("GOOGLE/DYNAMICWORLD/V1").select('label')
.filterDate(time_start, time_end)
.filterBounds(area).mode().eq(0).not()

var urban_mask = ee.ImageCollection("GOOGLE/DYNAMICWORLD/V1").select('label')
.filterDate(time_start, time_end)
.filterBounds(area).mode().eq(6).not()


var asc_min = asc_sigma.min()
var asc_max = asc_sigma.max()

var asc_sm = asc_sigma.map(function(img){
  var index = (img.subtract(asc_min)).divide(asc_max.subtract(asc_min))
  var date = img.date().format('YYYY-MM-dd')
  return index.multiply(water_mask).multiply(urban_mask).rename('sm_asc')
  .copyProperties(img,['system:time_start'])
  .set('date', ee.String(date))
  })
  
Map.addLayer(asc_sm.filterDate('2020','2021').toBands().clip(area),[],'sm_asc',false)

var des_min = des_sigma.min()
var des_max = des_sigma.max()

var des_sm = des_sigma.map(function(img){
  var index = (img.subtract(des_min)).divide(des_max.subtract(des_min))
  var date = img.date().format('YYYY-MM-dd')
  return index.multiply(water_mask).multiply(urban_mask).rename('sm_des')
  .copyProperties(img,['system:time_start'])
  .set('date', ee.String(date))
  })
  
Map.addLayer(des_sm.filterDate('2020','2021').toBands().clip(roi),[],'sm_des',false)

//----------------------- DATOS PROMEDIO SOIL MESUREMENT SAR ANUAL TEMPORADA SECA ------------------------------
var asc_sm_mean = asc_sm.filterDate('2019-12-01', '2020-04-30').mean().clip(roi);
var des_sm_mean = des_sm.filterDate('2019-12-01', '2020-04-30').mean().clip(roi);

Map.addLayer(asc_sm_mean, {min : 0, max: 0.7272940857048874}, 'asc_sm_mean', false);
Map.addLayer(des_sm_mean, {min : 0, max: 0.6719401096318013}, 'des_sm_mean', false);

//-------------------------------------- OPTRAM DATOS OPTICOS ------------------------------------------------
//Establezca la fecha o parámetros de inicio y fin de la búsqueda:

var sen_time_start = '2019-12-01';   // <--- fecha antes del inicio del proyecto.
var sen_time_end = '2020-04-30';     // <--- fecha inicio del proyecto.
// Establecer el umbral máximo de cobertura de nubes (0-100)

// The threshold for masking; values between 0.50 and 0.65 generally work well.
// Higher values will remove thin clouds, haze & cirrus shadows.
var CLEAR_THRESHOLD = 0.65

// Use 'cs' or 'cs_cdf', depending on your use case; see docs for guidance.
var QA_BAND = 'cs_cdf';

//**********************************************************************************************************************************//
//                                 PART 5: CARGAR COLECCIONES DE DATOS DE LAS IMAGENES

// Harmonized Sentinel-2 Level 2A collection.
var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED');

// Cloud Score+ image collection. Note Cloud Score+ is produced from Sentinel-2
// Level 1C data and can be applied to either L1C or L2A collections.
var csPlus = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED');

// Make a clear median composite.
var sentinel = s2
    .filterBounds(area)
    .filterDate(sen_time_start, sen_time_end)
    .linkCollection(csPlus, [QA_BAND])
    .map(function(img) {
      return img.updateMask(img.select(QA_BAND).gte(CLEAR_THRESHOLD));
    }).map(function(image){return image.clip(area)});


var parameters = sentinel.map(function(img){
  var bands = img.select('B.*').multiply(0.0001);
  var ndvi = bands.normalizedDifference(['B8','B4']).rename('ndvi');
  var str = bands.expression('((1 - swir) ** 2)/(2 * swir)',{'swir': bands.select('B12')}).rename('str');
  var stack = ee.Image.cat([ndvi, str]);
  return stack
  .copyProperties(img, img.propertyNames())
  });

var str_full_cover = parameters.map(function(img){
  var thr_full = img.select('ndvi').gt(0.3);
  var str_full = img.select('str').updateMask(thr_full);
  return str_full
  .copyProperties(img, img.propertyNames())
  });
  
var str_bareland = parameters.map(function(img){
  var thr_bare = img.select('ndvi').gte(0).and(img.select('ndvi').lt(0.2));
  var str_bare = img.select('str').updateMask(thr_bare)
  return str_bare
  .copyProperties(img, img.propertyNames())
  })
  
var scale = 30 
var vw = ee.Number(str_full_cover.max().reduceRegion({
  reducer: ee.Reducer.max(), geometry: roi, scale: scale
  }).values().get(0));
  
var vd = ee.Number(str_full_cover.min().reduceRegion({
  reducer: ee.Reducer.min(), geometry: roi, scale: scale
  }).values().get(0));
  

var iw = ee.Number(str_bareland.max().reduceRegion({
  reducer: ee.Reducer.max(), geometry: roi, scale: scale
  }).values().get(0)); 
  

var id = ee.Number(str_bareland.min().reduceRegion({
  reducer: ee.Reducer.min(), geometry: roi, scale: scale
  }).values().get(0)); 

var sw = vw.subtract(iw);
var sd = vd.subtract(id);

var sm = parameters.map(function(img){
  var mask = img.select('ndvi').lt(-0.1).not()
  var index = img.expression('(id + sd * ndvi - str) / (id - iw + (sd - sw) * ndvi)',
  {'id': id, 'sd': sd, 'ndvi': img.select('ndvi'), 'str': img.select('str'), 'iw': iw, 'sw': sw})
  .rename('soil_moisture');
  return index.multiply(1000).updateMask(mask)
  .copyProperties(img, img.propertyNames())
  })
  
print(
  ui.Chart.image.series(sm, muestra, ee.Reducer.first(), scale, 'system:time_start')
  )
  
Map.addLayer(sm.toBands().clip(roi),[],'sm_optical',false);

var sm_mean = sm.mean().clip(roi);
Map.addLayer(sm_mean,{min: 8.762528630450698, max : 274.7074238998578},'sm_opticalFirst',false);

//------------------------------------CLASIFICACIÓN SUPERVISADA----------------------------------

var clases = Turberas.merge(NoTurberas);

// Divide el conjunto de entrenamiento
var muestraConRandom = clases.randomColumn();                  // Agregar una columna 'random' para la división
var datosEntrenamiento = muestraConRandom.filter(ee.Filter.lte('random', 0.7));   // ~70% datos para entrenamiento
var datosValidacion = muestraConRandom.filter(ee.Filter.gt('random', 0.7));       // ~30% datos para validación

//Define both optical and SAR to train your data
var opt_sar = ee.Image.cat(asc_sm_mean, des_sm_mean, sm_mean);
var bands_opt_sar = ['sm_asc', 'sm_des', 'soil_moisture']
var training_opt_sar = opt_sar.select(bands_opt_sar).sampleRegions({
  collection: datosEntrenamiento,
  properties: ['b1'],
  scale: 30 });
  
//Train the classifier
var classifier = ee.Classifier.smileRandomForest(50).train({
  features: training_opt_sar, 
  classProperty: 'b1',
  inputProperties: bands_opt_sar 
});


// Calculate variable importance
var importance = ee.Dictionary(classifier.explain().get('importance'));

// Calculate relative importance
var sum = importance.values().reduce(ee.Reducer.sum());

var relativeImportance = importance.map(function(key, val) {
   return (ee.Number(val).multiply(100)).divide(sum);
  });
print(relativeImportance);

// Create a FeatureCollection so we can chart it
var importanceFc = ee.FeatureCollection([
  ee.Feature(null, relativeImportance)
]);

var chart = ui.Chart.feature.byProperty({
  features: importanceFc
}).setOptions({
      title: 'Feature Importance',
      vAxis: {title: 'Importance'},
      hAxis: {title: 'Feature'}
  });
print(chart);


var test = opt_sar.select(bands_opt_sar).sampleRegions({
  collection: datosValidacion,
  properties: ['b1'],
  scale: 30
});

// Tune the numberOfTrees parameter.
var numTreesList = ee.List.sequence(10, 150, 10);

var accuracies = numTreesList.map(function(numTrees) {
  var classifier = ee.Classifier.smileRandomForest(numTrees)
      .train({
        features: training_opt_sar,
        classProperty: 'b1',
        inputProperties: bands_opt_sar
      });

  // Here we are classifying a table instead of an image
  // Classifiers work on both images and tables
  return test
    .classify(classifier)
    .errorMatrix('b1', 'classification')
    .accuracy();
});

var chart = ui.Chart.array.values({
  array: ee.Array(accuracies),
  axis: 0,
  xLabels: numTreesList
  }).setOptions({
      title: 'Hyperparameter Tuning for the numberOfTrees Parameters',
      vAxis: {title: 'Validation Accuracy'},
      hAxis: {title: 'Number of Trees', gridlines: {count: 15}}
  });
print(chart);

// Tuning Multiple Parameters
// We can tune many parameters together using
// nested map() functions
// Let's tune 2 parameters
// numTrees and bagFraction 
var numTreesList = ee.List.sequence(10, 150, 10);
var bagFractionList = ee.List.sequence(0.1, 0.9, 0.1);

var accuracies = numTreesList.map(function(numTrees) {
  return bagFractionList.map(function(bagFraction) {
     var classifier = ee.Classifier.smileRandomForest({
       numberOfTrees: numTrees,
       bagFraction: bagFraction
     })
      .train({
        features: training_opt_sar,
        classProperty: 'b1',
        inputProperties: bands_opt_sar
      });

    // Here we are classifying a table instead of an image
    // Classifiers work on both images and tables
    var accuracy = test
      .classify(classifier)
      .errorMatrix('b1', 'classification')
      .accuracy();
    return ee.Feature(null, {'accuracy': accuracy,
      'numberOfTrees': numTrees,
      'bagFraction': bagFraction});
  });
}).flatten();
var resultFc = ee.FeatureCollection(accuracies);

// Alternatively we can automatically pick the parameters
// that result in the highest accuracy
var resultFcSorted = resultFc.sort('accuracy', false);
var highestAccuracyFeature = resultFcSorted.first();
print(highestAccuracyFeature);
var highestAccuracy = highestAccuracyFeature.getNumber('accuracy');
var optimalNumTrees = highestAccuracyFeature.getNumber('numberOfTrees');
var optimalBagFraction = highestAccuracyFeature.getNumber('bagFraction');
//var properties = resultFc.getInfo();


// Use the optimal parameters in a model and perform final classification
var optimalModel = ee.Classifier.smileRandomForest({
  numberOfTrees: optimalNumTrees,
  bagFraction: optimalBagFraction
}).train({
  features: training_opt_sar, 
  classProperty: 'b1',
  inputProperties: bands_opt_sar
});

var visClass =({
  min: 0,
  max: 1,
  palette: ['yellow', 'green']
});
var finalClassification = opt_sar.select(bands_opt_sar).classify(optimalModel);
//var bosque = finalClassification.updateMask(finalClassification.neq(0));
Map.addLayer(finalClassification, visClass, 'Clasificacion Optimizada', 1);


// Clasificar el conjunto de datos de validación usando el clasificador entrenado
var classifiedValidation = test.classify(optimalModel);

// Calcular la matriz de confusión
var errorMatrix = classifiedValidation.errorMatrix('b1', 'classification');

// Calcular el índice Kappa
var kappa = errorMatrix.kappa();

// Imprimir la matriz de confusión y el índice Kappa
print('Matriz de Confusión:', errorMatrix);
print('Índice Kappa:', kappa);




