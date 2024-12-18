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

# ANALISIS HUMEDALES

Map.addLayer({
  eeObject: ubicacion,
  name: 'Municipio de Santa Luzia do Itanhy',
  shown: 0
});

// Localizacion
var area = ee.FeatureCollection(ubicacion);

// Establezca el área de estudio como centro del mapa.
Map.centerObject(area);

var dataset = ee.Image('JRC/GSW1_4/GlobalSurfaceWater').clip(area);

var seasonality = dataset.select('seasonality').lte(6); // Seleccionamos la banda 'occurrence' y aplicamos la máscara

// Aplicamos la máscara a la imagen del dataset
var seasonalityUpdate = dataset.updateMask(seasonality).toUint16();

var visualizationSeasonality = {
  bands: ['seasonality'],
  min: 0,
  max: 6,
    palette: ['ffffff', 'ffbbbb']
};

//Map.addLayer(seasonalityUpdate, visualizationSeasonality, 'seasonality');

var transition = {
  bands: ['transition'],
  min: 0,
  max: 10,
  palette: ['ffffff', '#0000ff', '#22b14c', '#d1102d', '#99d9ea', '#b5e61d', '#e6a1aa', '#ff7f27', '#ffc90e', '#7f7f7f', '#c3c3c3']
};

Map.addLayer(seasonalityUpdate, transition, 'Transition');

Export.image.toDrive({
  image: seasonalityUpdate.select('transition'), 
  description: 'Humedales_Brasil',
  folder: 'DESCARGAS_GEE', // Cambiar folder
  fileNamePrefix: 'Humedales_Brasil', // Cambiar nombre predio
  region: area, 
  scale: 30, 
  maxPixels: 1e10
});

# CODIGOS MES DE OCTUBRE

## CODIGO NUMERO 1

//===============================================================================================================================
//                                     MODELO CLASIFICACIÓN BnB 
//===============================================================================================================================


//===============================================================================================================================//

//*******************************************************************************************************************************
//                            PART 1: SELECCIONA,DIBUJA IMPORTA EL AREA A ANALIZAR
// incialmente debemos importar o dibujar una geometria del area donde se encuentran los predios viables para interpretación 
// esta se puede hacer con la herramienta de dibujo de geometrias o importando un shapefile o archivo CSV esto con el fin de deli-
// mitar la zona de estudio, una vez hecho esto procedemos a colocarle nombre en la zona de archivos importados.
Map.addLayer({
  eeObject:aoi3, 
  name: 'Areas Analisis', 
  shown: 0, 
});

/******************************************************************************************************************************/
//                            PART 2: IDENTIFICAR Y EJECUTAR ENTRADAS DE USUARIO


// Localizacion
var area = ee.FeatureCollection(aoi3);


// Establezca el área de estudio como centro del mapa.
Map.centerObject(area);

//*******************************************************************************************************************************
//                                  PART 4: ESTABLECE LAS FECHAS DE BÚSQUEDA

/*
Establezca las fechas de inicio y finalización de un período ANTES de la fecha de INICIO del proyecto, esto se recomienda que para landsat 8 sea de 3 meses
por análisis mientras que para sentinel 2 sea de cada mes esto la disponibilidad de imágenes que captan los satélites y en el detalle 
del analisis.

*/


//Establezca la fecha o parámetros de inicio y fin de la búsqueda:

var time_start = '2020-01-01';   // <--- fecha antes del inicio del proyecto.
var time_end = '2020-12-30';     // <--- fecha inicio del proyecto.
// Establecer el umbral máximo de cobertura de nubes (0-100)

// The threshold for masking; values between 0.50 and 0.65 generally work well.
// Higher values will remove thin clouds, haze & cirrus shadows.
var CLEAR_THRESHOLD = 0.65;

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
var collection = s2
    .filterBounds(area)
    .filterDate(time_start, time_end)
    .linkCollection(csPlus, [QA_BAND])
    .map(function(img) {
      return img.updateMask(img.select(QA_BAND).gte(CLEAR_THRESHOLD));
    }).map(function(image){return image.clip(area)});

print(collection);
var vis = {bands: ['B4', 'B3', 'B3'], min: 224 , max: 1709, gamma: 1};
Map.addLayer(collection, vis,'Coleccion de imagenes', 0);

var image = collection.median();
Map.addLayer(image, {bands: ['B4', 'B3', 'B2'], min: 223.5 , max: 1445, gamma: 1},'Mosaico', 1);
/******************************************************************************************************************************/
//                                     PART 6: AGREGAR INDICES ESPECTRALES

var addIndices = function(image) {
  var NDVI = image.normalizedDifference(['B8', 'B4']).rename(['NDVI']);
  return image.addBands(NDVI)};

var composite = addIndices(image);

Map.addLayer(composite, {min: 0.15055880479471886, max: 0.8577648766328012, bands: ['NDVI']},'NDVI', 0);

// Reclasificar el NDVI usando una función de mapeo condicional
var NDVIReclass = composite.select('NDVI').expression(
  "b(0) > 0.68 ? 4 : " +                    // Valores mayores a 0.68 (Áreas de bosque)
  "b(0) > 0.64 && b(0) <= 0.68 ? 3 : " +    // Valores entre 0.64 y 0.68
  "b(0) > 0.38 && b(0) <= 0.65 ? 2 : " +    // Valores entre 0.38 y 0.65
  "b(0) > 0.06 && b(0) <= 0.38 ? 1 : " +    // Valores entre 0.06 y 0.38
  "b(0) <= 0.06 ? 0 : 0"                    // Valores menores o iguales a 0.06
).rename('NDVI_Reclass').clip(area);

Map.addLayer(NDVIReclass, {min: 0, max: 4, palette: ['FFFF00', 'F2FF00', 'F2FF00', 'F2FF00', '006400']},'NDVIReclass', 0);

/*
var addIndices = function(image) {
    var ARVI = image.expression ('((NIR - (2 * RED) + BLUE) / (NIR + (2 * RED) + BLUE))',{
        'NIR': image.select ('B8'),
        'BLUE': image.select ('B2'), 
        'RED': image.select ('B4')
    }).rename("ARVI");
  return image.addBands(ARVI)};

var composite = addIndices(image);

Map.addLayer(composite, {min: 0.15055880479471886, max: 0.8577648766328012, bands: ['ARVI']},'ARVI', 0);

// Reclasificar el NDVI usando una función de mapeo condicional
var ARVIReclass = composite.select('ARVI').expression(
  "b(0) > 0.8 && b(0) <= 1 ? 2 : " +      //Vegetación 
  "b(0) > 0.6 && b(0) <= 0.8 ? 1 : " +    // Vegetación sana
  "b(0) > 0.6 && b(0) <= -1 ? 0 : 0"   // Sequía, superficies sin agua                   // Valores menores o iguales a 0.06
).rename('NDVI_Reclass').clip(area);

Map.addLayer(ARVIReclass, {min: 0, max: 2, palette: ['FFFF00', '#3fca29', '006400']},'ARVIReclass', 1);
*/


Export.image.toDrive({
  image: NDVIReclass, 
  description: 'Insumo_Paramuno', 
  folder: 'DESCARGAS_GEE',
  fileNamePrefix: 'Insumo_Paramuno',
  region: area, 
  scale: 10, 
  maxPixels: 1e10
  
});

## CODIGO NUMERO 2


//===============================================================================================================================
//                                                 MODELO MONITOREO AREAS DE INCENDIO
//===============================================================================================================================
// Este codigo muestra la relación de quemado normalizada se aplicará a las imágenes de antes y después de un incendio forestal.
// Esto calculando la diferencia después (dNBR) se deriva Burn Severity, que muestra el espacio impacto de la perturbación. 
// Las imágenes utilizadas en este proceso provienen de Sentinel-2 o Landsat 8.
//===============================================================================================================================//

//***********************************************UBICA EL ÁREA DE ESTUDIO********************************************************//



//*******************************************************************************************************************************
//                            PART 1: SELECCIONA,DIBUJA IMPORTA EL AREA A ANALIZAR
// incialmente debemos importar o dibujar una geometria del area donde queremos realizar el analisis por afectaciones de fuego,
// esta se puede hacer con la herramienta de dibujo de geometrias o importando un shapefile o archivo CSV esto con el fin de deli-
// mitar la zona de estudio, una vez hecho esto procedemos a colocarle nombre en la zona de archivos importados y la llamamos table.

//*******************************************************************************************************************************
//                                  PART 2: ESTABLECE LAS FECHAS DE ANALISIS

//Establezca las fechas de inicio y finalización de un período ANTES del incendio, esto se recomienda que para landsat 8 sea de 3 meses
//por análisis mientras que para sentinel 2 sea de cada mes esto la disponibilidad de imágenes que captan los satélites y en el detalle 
//del analisis.

//Establezca la fecha o parámetros de inicio y fin de los incendios
var prefire_start = '2024-01-01';   // <--- fecha de Inicio del incendio
var prefire_end = '2024-05-07';     // <--- fecha Fin del incendio

// Establecer la fecha o parametros despues del incendio
var postfire_start = '2024-01-01'; // <--- fecha después del incendio
var postfire_end = '2024-05-07';   // <--- fecha después del incendio

//******************************************************************************************************************************
//                                     PART 3: SELECCIONA EL SATELITE A USAR

// Puede seleccionar imágenes de detección remota de dos sensores satelitales disponibles.
// Considere los detalles de cada misión a continuación para elegir los datos adecuados para sus necesidades.

// Landsat 8  (L8)                       |  Sentinel-2   (S2)
//-------------------------------------------------------------------------------------------
// lanzamiento:     Febrero  11, 2015    |  Junio 23, 2015 & Marzo 7, 2017
// tasa de repetición: 16 dias           |  5 dias (desde 2017)
// resolucion:      30 metros            |  10 metros
// ventajas:    serie temporal más larga | 9 veces mayor detalle espacial
//    archivo de exportación más pequeño | mayor probabilidad de imágenes sin nubes

//Selecciona uno de los siguientes:   'L8'  o 'S2' 

var platform = 'S2';        // <--- asigna tu satelite a la variable "" de plataforma recuerda que puedes elegir si 'L8'  o 'S2'
                                 // si deseas usar el satelite landsat 8 escribe L8 o si quieres trabajar con sentinel 2 
                                 // escribe S2.

//*******************************************************************************************************************************

//                               PARTE 4: IDENTIFICAR Y EJECUTAR ENTRADAS DE USUARIO

//---------------------------------- Traducción de las entradas del usuario-------------------------------------------------------

// Imprimir plataforma Satélite y fechas a consola
if (platform == 'S2' | platform == 's2') {
  var ImCol = 'COPERNICUS/S2_SR_HARMONIZED';  // <---Satelite usado para sentinel 2
  var pl = 'Sentinel-2';
} else {
  var ImCol = 'LANDSAT/LC08/C01/T1_SR'; //<---Satelite usado para landsat 8
  var pl = 'Landsat 8';
}
print(ee.String('Satelite seleccionado para el análisis: ').cat(pl));
print(ee.String('Incendio ocurrido entre ').cat(prefire_end).cat(' y ').cat(postfire_start));

// Localizacion
var area = ee.FeatureCollection(ALTAMIRA);

// Establezca el área de estudio como centro del mapa.
Map.centerObject(area);
Map.addLayer({
  eeObject: area, 
  name: 'Areas Analisis', 
  shown: 0, 
});


//----------------------- Seleccionar imágenes por hora, ubicación y nubosidad ---------------------------------------------------
var imagery = ee.ImageCollection(ImCol);

// Establecer el umbral máximo de cobertura de nubes (0-100)
var maxCloudCover = 30; // <--- Ingrese el % de nubeS, entre mas % aumentan la disponiblidad de imagenes pero tambien la nubes en ellas

// En las siguientes líneas, las imágenes se recopilarán en un ImageCollection, según el
// ubicación de nuestra área de estudio, un marco de tiempo determinado y la proporción de cobertura de nubes.
var prefireImCol = ee.ImageCollection(imagery
    // Filtrar por fechas.
    .filterDate(prefire_start, prefire_end)
    // Filtrar por ubicación.
    .filterBounds(area)
    // Filtrar por nubosidad.
    .filter(ee.Filter.lt('CLOUD_COVER', maxCloudCover)));

// Selecciona todas las imágenes que se superponen con el área de estudio de un marco de tiempo dado
var postfireImCol = ee.ImageCollection(imagery
    // Filtrar por fechas.
    .filterDate(postfire_start, postfire_end)
    // Filtrar por ubicación.
    .filterBounds(area)
    // Filtrar por nubosidad.
    .filter(ee.Filter.lt('CLOUD_COVER', maxCloudCover)));

// Agrega las imágenes recortadas a la consola de la derecha
print("Pre-fire Image Collection: ", prefireImCol); 
print("Post-fire Image Collection: ", postfireImCol);

//**********************************************************************************************************************************
//                           PART 5: APLICAR UNA MASCARA DE NUBES,NIEVE Y SOMBRAS

// Función para enmascarar nubes de la banda de calidad de píxeles de los datos de Sentinel-2 SR.
function maskS2sr(image) {
  // Los bits 10 y 11 son nubes y cirros, respectivamente.
  var cloudBitMask = ee.Number(2).pow(10).int();
  var cirrusBitMask = ee.Number(2).pow(11).int();
  // Obtener la banda de control de calidad de píxeles.
  var qa = image.select('QA60');
  // Todas las banderas deben establecerse en cero, lo que indica condiciones claras.
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
      .and(qa.bitwiseAnd(cirrusBitMask).eq(0));
  // Devuelve la imagen enmascarada, escalada a reflectancia TOA, sin las bandas QA.
  return image.updateMask(mask)
      .copyProperties(image, ["system:time_start"]);
}

// Función para enmascarar nubes de la banda de calidad de píxeles de los datos Landsat 8 SR.
function maskL8sr(image) {
  // Los bits 3 y 5 son sombra de nube y nube, respectivamente.
  var cloudShadowBitMask = 1 << 3;
  var cloudsBitMask = 1 << 5;
  var snowBitMask = 1 << 4;
  // Obtener la banda de control de calidad de píxeles.
  var qa = image.select('pixel_qa');
  // Todas las banderas deben establecerse en cero, lo que indica condiciones claras.
  var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
      .and(qa.bitwiseAnd(cloudsBitMask).eq(0))
      .and(qa.bitwiseAnd(snowBitMask).eq(0));
  // Devuelve la imagen enmascarada, escalada a reflectancia TOA, sin las bandas QA.
  return image.updateMask(mask)
      .select("B[0-9]*")
      .copyProperties(image, ["system:time_start"]);
}

// Aplicar máscara de nube específica de la plataforma
if (platform == 'S2' | platform == 's2') {
  var prefire_CM_ImCol = prefireImCol.map(maskS2sr);
  var postfire_CM_ImCol = postfireImCol.map(maskS2sr);
} else {
  var prefire_CM_ImCol = prefireImCol.map(maskL8sr);
  var postfire_CM_ImCol = postfireImCol.map(maskL8sr);
}


//**************************************************************************************************************************************
//                                              PARTE 6: IMÁGEN INDIVIDUAL Y CLIP DE IMÁGENES AL ÁREA DE ESTUDIO

var pre_mos = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20220331T151659_20220331T152018_T18NZL').clip(area); ////---> SE DEBE AJUSTAR A UNA IMÁGEN UNICA
var post_mos = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20220823T151721_20220823T151717_T18NZL').clip(area);

var pre_cm_mos = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20220331T151659_20220331T152018_T18NZL').clip(area);
var post_cm_mos = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20220823T151721_20220823T151717_T18NZL').clip(area);

// Agrega las imágenes recortadas a la consola
print("Imagen en color verdadero antes del incendio: ", pre_mos); 
print("Imagen de color verdadero posterior al incendio: ", post_mos);

//****************************************************************************************************************************************
//                                                PARTE 7: CALCULAR NBR PARA IMÁGENES ANTES Y POST-FUEGO

// Aplicar NBR específico de la plataforma = (NIR-SWIR2) / (NIR+SWIR2)
if (platform == 'S2' | platform == 's2') {
  var preNBR = pre_cm_mos.normalizedDifference(['B8','B12']);
  var postNBR = post_cm_mos.normalizedDifference(['B8','B12']);
} else {
  var preNBR = pre_cm_mos.normalizedDifference(['SR_B5', 'SR_B7']); 
  var postNBR = post_cm_mos.normalizedDifference(['SR_B5', 'SR_B7']);
}
// Agregue las imágenes NBR a la consola 
print("Ratio de quema normalizado previo al incendio: ", preNBR);
print("Proporción de area quemada normalizada posterior al incendio: ", postNBR);

//***************************************************************************************************************************************
//                                               PART 8: CALCULE LA DIFERENCIA ENTRE IMÁGENES ANTES Y POST-FUEGO

// Se llama el resultado delta NBR o dNBR
var dNBR_unscaled = preNBR.subtract(postNBR);

// Escale el producto según los estándares de USGS
var dNBR = dNBR_unscaled.multiply(1000);

// Agrega la imagen de diferencia a la consola
print("Diferencia en la relación de quema normalizada: ", dNBR);

//***************************************************************************************************************************************
//                                                      PART 9: AGREGAR CAPAS AL MAPA

// Agregar límite.
Map.addLayer(area.draw({color: 'ffffff', strokeWidth: 5}), {},'Area de estudio');

//---------------------------------- Imágenes en color verdadero ------------------------------------

// Aplicar parámetros de visualización específicos de la plataforma para imágenes en color verdadero
if (platform == 'S2' | platform == 's2') {
  var vis = {bands: ['B4', 'B3', 'B2'], max: 2000, gamma: 1.5};
} else {
  var vis = {bands: ['B4', 'B3', 'B2'], min: 0, max: 4000, gamma: 1.5};
}

// Agrega las imágenes en color verdadero al mapa.
Map.addLayer(pre_mos, vis,'Imagen previa al incendio');
Map.addLayer(post_mos, vis,'Imagen posterior al incendio');

// Agregue las imágenes en color verdadero al mapa con la mascara de nubes.
Map.addLayer(pre_cm_mos, vis,'Imagen en color verdadero antes del disparo - Nubes enmascaradas');
Map.addLayer(post_cm_mos, vis,'Imagen de color verdadero posterior al incendio - Nubes enmascaradas');

//--------------------------- Producto de relación de quemado: escala de grises -------------------------------

var grey = ['white', 'black'];
Map.addLayer(preNBR, {min: -1, max: 1, palette: grey}, 'Relación de quema normalizada antes del incendio');
Map.addLayer(postNBR, {min: -1, max: 1, palette: grey}, 'Relación de quema normalizada despues del incendio');
Map.addLayer(dNBR, {min: -1000, max: 1000, palette: grey}, 'escala de grises dNBR');

//------------------------- Producto de relación de quema - Clasificación USGS ----------------------------

// Definir un estilo SLD de intervalos discretos para aplicar a la imagen.
var sld_intervals =
  '<RasterSymbolizer>' +
    '<ColorMap type="intervals" extended="false" >' +
      '<ColorMapEntry color="#ffffff" quantity="-500" label="-500"/>' +
      '<ColorMapEntry color="#7a8737" quantity="-250" label="-250" />' +
      '<ColorMapEntry color="#acbe4d" quantity="-100" label="-100" />' +
      '<ColorMapEntry color="#0ae042" quantity="100" label="100" />' +
      '<ColorMapEntry color="#fff70b" quantity="270" label="270" />' +
      '<ColorMapEntry color="#ffaf38" quantity="440" label="440" />' +
      '<ColorMapEntry color="#ff641b" quantity="660" label="660" />' +
      '<ColorMapEntry color="#a41fd6" quantity="2000" label="2000" />' +
    '</ColorMap>' +
  '</RasterSymbolizer>';

// Agregue la imagen al mapa utilizando tanto la rampa de color como los esquemas de intervalo.
Map.addLayer(dNBR.sldStyle(sld_intervals), {}, 'dNBR CLASIFICADO');

// Seperate result into 8 burn severity classes
var thresholds = ee.Image([-1000, -251, -101, 99, 269, 439, 659, 2000]);
var classified = dNBR.lt(thresholds).reduce('sum').toInt();

//***************************************************************************************************************************************
//                                           PART 10: AGREGAR ESTADÍSTICAS DE ÁREA QUEMADA

// Contar el número de píxeles en toda la capa
var allpix =  classified.updateMask(classified);  // enmascara toda la capa
var pixstats = allpix.reduceRegion({
  reducer: ee.Reducer.count(),               // cuenta píxeles en una sola clase
  geometry: area,
  scale: 30
  });
var allpixels = ee.Number(pixstats.get('sum')); // extrae el recuento de píxeles como un número


// crea una lista vacía para almacenar valores de área 
var arealist = [];

// crea una función para derivar la extensión de una clase de gravedad de quemado
// los argumentos son el número de clase y el nombre de clase
var areacount = function(cnr, name) {
 var singleMask =  classified.updateMask(classified.eq(cnr));  // enmascara una sola clase
 var stats = singleMask.reduceRegion({
  reducer: ee.Reducer.count(),               // cuenta píxeles en una sola clase
  geometry: area,
  scale: 30
  });
var pix =  ee.Number(stats.get('sum'));
var hect = pix.multiply(900).divide(10000);               // Píxel Landsat = 30m x 30m --> 900 m2
var perc = pix.divide(allpixels).multiply(10000).round().divide(100);   // obtiene el porcentaje de área por clase y redondea a 2 decimales
arealist.push({Class: name, Pixels: pix, Hectares: hect, Percentage: perc});
};

// clases de gravedad en diferente orden
var names2 = ['NA', 'Gravedad alta', 'Gravedad moderada-alta',
'Gravedad moderada-baja', 'Gravedad baja', 'No quemado', 'Rebrote mejorado, bajo', 'Rebrote mejorado, alto'];

//ejecutar funcion para cada clase
for (var i = 0; i < 8; i++) {
  areacount(i, names2[i]);
  }

print('Área quemada por clase de gravedad', arealist, '--> haga clic en la lista de objetos para clases individuales');

//***************************************************************************************************************************************
//                                               PART 11: AÑADIR UNA LEYENDA
// Establecer la posición del panel
var legend = ui.Panel({
  style: {
    position: 'bottom-left',
    padding: '8px 15px'
  }});
 
// Crear título de leyenda
var legendTitle = ui.Label({
  value: 'Clases de areas afectadas por severidad',
  style: {fontWeight: 'bold',
    fontSize: '18px',
    margin: '0 0 4px 0',
    padding: '0'
    }});
 
// Agregar el título al panel
legend.add(legendTitle);
 
// Crea y aplica estilo a 1 fila de la leyenda.
var makeRow = function(color, name) {
 
      // Cree la etiqueta que en realidad es el cuadro de color.
      var colorBox = ui.Label({
        style: {
          backgroundColor: '#' + color,
          // Usa padding para dar la altura y el ancho de la caja.
          padding: '8px',
          margin: '0 0 4px 0'
        }});
 
      // Cree la etiqueta rellena con el texto de descripción.
      var description = ui.Label({
        value: name,
        style: {margin: '0 0 4px 6px'}
      });
 
      // Devuelve el panel
      return ui.Panel({
        widgets: [colorBox, description],
        layout: ui.Panel.Layout.Flow('horizontal')
      })};
 
// Paleta con los colores
//            verde oliva, verde amarillento, verde lima brillante, amarillo, naranja, naranja oscuro, violeta, Blanco
var palette =['7a8737', 'acbe4d', '0ae042', 'fff70b', 'ffaf38', 'ff641b', 'a41fd6', 'ffffff'];
 
// Nombre de la leyenda
var names = ['7-Rebrote mejorado, alto', '6-Rebrote mejorado, bajo', '5-Sin quemar', '4-Gravedad baja',
'3-Gravedad moderada-baja', '2-Gravedad moderada-alta', '1-Gravedad alta', '0-NA'];
 
// Agregar color y nombres
for (var i = 0; i < 8; i++) {
  legend.add(makeRow(palette[i], names[i]));
  }  
 
// Agregar leyenda al mapa (alternativamente, también puede imprimir la leyenda en la consola)
Map.add(legend);

//***************************************************************************************************************************************
//                                               PART 12: PREPARAR EXPORTACIÓN DE ARCHIVO
//Asegurate que al momento de descargar el archivo le ingreses el codigo del sistema de coordenadas que deseas, como se llamara
// el archivo, en que carpeta se va alojar y como aparecera el nombre del archivo en el panel de descargas
var id = dNBR.id().getInfo();
      
Export.image.toDrive({image: classified, scale: 30, description:'dNBR-Clasificacion' , fileNamePrefix: 'dNBR-Clasificacion',region: area, maxPixels: 1e10}); // <--- Descarga la clasificacion
Export.image.toDrive({image: post_cm_mos.visualize(vis), scale: 30, description: 'dNBR-ColorverdaderoDespuesdelincendio', fileNamePrefix: 'dNBR-Color verdadero Despues del incendio ',region: area, maxPixels: 1e10}); // <--- Descarga la imagen despues del incendio

Map.addLayer(roi);

# CODIGOS MES DE NOVIEMBRE

## CODIGO NUMERO 1

Map.addLayer({
  eeObject:CO2BIOP2, 
  name: 'CO2BIO_P2', 
  shown: 0, 
});

//var dt = ind_dt.filter(ee.Filter.eq('DISTRICT','GOALPARA'));

var dt = ee.FeatureCollection(CO2BIOP2);
// Establezca el área de estudio como centro del mapa.
Map.centerObject(dt);

var collection = ee.ImageCollection('COPERNICUS/S1_GRD')
        .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
        .filter(ee.Filter.eq('instrumentMode', 'IW'))
        .filter(ee.Filter.or(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'),ee.Filter.eq('orbitProperties_pass', 'ASCENDING')));
        
        
var before = collection.filter(ee.Filter.date('2021-01-01', '2022-12-30')).filterBounds(dt);    
var after = collection.filter(ee.Filter.date('2022-01-01', '2023-12-30')).filterBounds(dt);

print(before);
print(after);

var before_image = before.select('VH').median().clip(dt);
var after_image = after.select('VH').median().clip(dt);

var before_filtered = ee.Image(toDB(RefinedLee(toNatural(before_image))));
var after_filtered = ee.Image(toDB(RefinedLee(toNatural(after_image))));

var flood = before_filtered.gt(-20).and(after_filtered.lt(-20));
var flood_mask = flood.updateMask(flood.eq(1));

var water = before_filtered.lt(-20).and(after_filtered.lt(-20));
var water_mask = water.updateMask(water.eq(1));

Map.centerObject(dt);
//Map.addLayer(before_image,{min:-25,max:0},'before')
//Map.addLayer(after_image,{min:-25,max:0},'after')
Map.addLayer(before_filtered,{min:-25,max:0},'before_filtered');
Map.addLayer(after_filtered,{min:-25,max:0},'after_filtered');
Map.addLayer(flood_mask,{palette:['Yellow']},'Flood_Inundation');
Map.addLayer(water_mask,{palette:['Blue']},'Water');

print('Total District Area (Ha)', dt.geometry().area().divide(10000));

var stats = flood_mask.multiply(ee.Image.pixelArea()).reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: dt,
  scale: 10,
  maxPixels: 1e13,
  tileScale: 16
});

print(stats);
/*var flood_area = ee.Number(stats.get('sum')).divide(10000).round();
print('Flooded Area (Ha)', flood_area);*/


//=============================   EXPORTAR LA INUNDACIÓN EN FORMATO TIFF ================================

Export.image.toDrive({
  image: flood_mask, 
  description: 'Inundacion',
  folder: 'DESCARGAS_GEE',                                                              // Cambiar folder
  fileNamePrefix: 'Inundacion',                                                   // Cambiar nombre predio
  region: dt, 
  scale: 10, 
  maxPixels: 1e11

});

Export.image.toDrive({
  image: water_mask, 
  description: 'Cuerpos_de_agua',
  folder: 'DESCARGAS_GEE',                                                              // Cambiar folder
  fileNamePrefix: 'Cuerpos_de_agua',                                                   // Cambiar nombre predio
  region: dt, 
  scale: 10, 
  maxPixels: 1e11

});
// set position of panel
var legend = ui.Panel({
  style: {
    position: 'bottom-left',
    padding: '8px 15px'
  }
});
 
// Create legend title
var legendTitle = ui.Label({
  value: 'Legend',
  style: {
    fontWeight: 'bold',
    fontSize: '18px',
    margin: '0 0 4px 0',
    padding: '0'
    }
});
 
// Add the title to the panel
legend.add(legendTitle);
 
// Creates and styles 1 row of the legend.
var makeRow = function(color, name) {
 
      // Create the label that is actually the colored box.
      var colorBox = ui.Label({
        style: {
          backgroundColor: color,
          // Use padding to give the box height and width.
          padding: '8px',
          margin: '0 0 4px 0'
        }
      });
 
      // Create the label filled with the description text.
      var description = ui.Label({
        value: name,
        style: {margin: '0 0 4px 6px'}
      });
 
      // return the panel
      return ui.Panel({
        widgets: [colorBox, description],
        layout: ui.Panel.Layout.Flow('horizontal')
      });
};
 
//  Palette with the colors
var palette =['Yellow', 'Blue'];
 
// name of the legend
var names = ['Flood Inundation','Water body'];
 
// Add color and and names
for (var i = 0; i <2; i++) {
  legend.add(makeRow(palette[i], names[i]));
  }  
 
// add legend to map (alternatively you can also print the legend to the console)
Map.add(legend);

// set position of panel
var title = ui.Panel({
  style: {
    position: 'top-center',
    padding: '8px 15px'
  }
});
 
// Create legend title
var mapTitle = ui.Label({
  value: 'Inundación en el área de interés 2018 - 2021',
  style: {
    fontWeight: 'bold',
    fontSize: '18px',
    margin: '0 0 4px 0',
    padding: '0'
    }
});
 
// Add the title to the panel
title.add(mapTitle);
Map.add(title);


//############################
// Speckle Filtering Functions
//############################

// Function to convert from d
function toNatural(img) {
  return ee.Image(10.0).pow(img.select(0).divide(10.0));
}

//Function to convert to dB
function toDB(img) {
  return ee.Image(img).log10().multiply(10.0);
}

//Apllying a Refined Lee Speckle filter as coded in the SNAP 3.0 S1TBX:

//https://github.com/senbox-org/s1tbx/blob/master/s1tbx-op-sar-processing/src/main/java/org/esa/s1tbx/sar/gpf/filtering/SpeckleFilters/RefinedLee.java
//Adapted by Guido Lemoine

// by Guido Lemoine
function RefinedLee(img) {
  // img must be in natural units, i.e. not in dB!
  // Set up 3x3 kernels 
  var weights3 = ee.List.repeat(ee.List.repeat(1,3),3);
  var kernel3 = ee.Kernel.fixed(3,3, weights3, 1, 1, false);

  var mean3 = img.reduceNeighborhood(ee.Reducer.mean(), kernel3);
  var variance3 = img.reduceNeighborhood(ee.Reducer.variance(), kernel3);

  // Use a sample of the 3x3 windows inside a 7x7 windows to determine gradients and directions
  var sample_weights = ee.List([[0,0,0,0,0,0,0], [0,1,0,1,0,1,0],[0,0,0,0,0,0,0], [0,1,0,1,0,1,0], [0,0,0,0,0,0,0], [0,1,0,1,0,1,0],[0,0,0,0,0,0,0]]);

  var sample_kernel = ee.Kernel.fixed(7,7, sample_weights, 3,3, false);

  // Calculate mean and variance for the sampled windows and store as 9 bands
  var sample_mean = mean3.neighborhoodToBands(sample_kernel); 
  var sample_var = variance3.neighborhoodToBands(sample_kernel);

  // Determine the 4 gradients for the sampled windows
  var gradients = sample_mean.select(1).subtract(sample_mean.select(7)).abs();
  gradients = gradients.addBands(sample_mean.select(6).subtract(sample_mean.select(2)).abs());
  gradients = gradients.addBands(sample_mean.select(3).subtract(sample_mean.select(5)).abs());
  gradients = gradients.addBands(sample_mean.select(0).subtract(sample_mean.select(8)).abs());

  // And find the maximum gradient amongst gradient bands
  var max_gradient = gradients.reduce(ee.Reducer.max());

  // Create a mask for band pixels that are the maximum gradient
  var gradmask = gradients.eq(max_gradient);

  // duplicate gradmask bands: each gradient represents 2 directions
  gradmask = gradmask.addBands(gradmask);

  // Determine the 8 directions
  var directions = sample_mean.select(1).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(7))).multiply(1);
  directions = directions.addBands(sample_mean.select(6).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(2))).multiply(2));
  directions = directions.addBands(sample_mean.select(3).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(5))).multiply(3));
  directions = directions.addBands(sample_mean.select(0).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(8))).multiply(4));
  // The next 4 are the not() of the previous 4
  directions = directions.addBands(directions.select(0).not().multiply(5));
  directions = directions.addBands(directions.select(1).not().multiply(6));
  directions = directions.addBands(directions.select(2).not().multiply(7));
  directions = directions.addBands(directions.select(3).not().multiply(8));

  // Mask all values that are not 1-8
  directions = directions.updateMask(gradmask);

  // "collapse" the stack into a singe band image (due to masking, each pixel has just one value (1-8) in it's directional band, and is otherwise masked)
  directions = directions.reduce(ee.Reducer.sum());  

  //var pal = ['ffffff','ff0000','ffff00', '00ff00', '00ffff', '0000ff', 'ff00ff', '000000'];
  //Map.addLayer(directions.reduce(ee.Reducer.sum()), {min:1, max:8, palette: pal}, 'Directions', false);

  var sample_stats = sample_var.divide(sample_mean.multiply(sample_mean));

  // Calculate localNoiseVariance
  var sigmaV = sample_stats.toArray().arraySort().arraySlice(0,0,5).arrayReduce(ee.Reducer.mean(), [0]);

  // Set up the 7*7 kernels for directional statistics
  var rect_weights = ee.List.repeat(ee.List.repeat(0,7),3).cat(ee.List.repeat(ee.List.repeat(1,7),4));

  var diag_weights = ee.List([[1,0,0,0,0,0,0], [1,1,0,0,0,0,0], [1,1,1,0,0,0,0], 
    [1,1,1,1,0,0,0], [1,1,1,1,1,0,0], [1,1,1,1,1,1,0], [1,1,1,1,1,1,1]]);

  var rect_kernel = ee.Kernel.fixed(7,7, rect_weights, 3, 3, false);
  var diag_kernel = ee.Kernel.fixed(7,7, diag_weights, 3, 3, false);

  // Create stacks for mean and variance using the original kernels. Mask with relevant direction.
  var dir_mean = img.reduceNeighborhood(ee.Reducer.mean(), rect_kernel).updateMask(directions.eq(1));
  var dir_var = img.reduceNeighborhood(ee.Reducer.variance(), rect_kernel).updateMask(directions.eq(1));

  dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), diag_kernel).updateMask(directions.eq(2)));
  dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), diag_kernel).updateMask(directions.eq(2)));

  // and add the bands for rotated kernels
  for (var i=1; i<4; i++) {
    dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)));
    dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)));
    dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)));
    dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)));
  }

  // "collapse" the stack into a single band image (due to masking, each pixel has just one value in it's directional band, and is otherwise masked)
  dir_mean = dir_mean.reduce(ee.Reducer.sum());
  dir_var = dir_var.reduce(ee.Reducer.sum());

  // A finally generate the filtered value
  var varX = dir_var.subtract(dir_mean.multiply(dir_mean).multiply(sigmaV)).divide(sigmaV.add(1.0));

  var b = varX.divide(dir_var);

  var result = dir_mean.add(b.multiply(img.subtract(dir_mean)));
  return(result.arrayFlatten([['sum']]));
}

## CÓDIGO NÚMERO 2

//===============================================================================================================================
//                                MODELO DESCARGA DATOS NDWI PARA UNA REGIÓN DE INTERÉS
//===============================================================================================================================


//===============================================================================================================================//

//*******************************************************************************************************************************
//                            PART 1: SELECCIONA,DIBUJA IMPORTA EL AREA A ANALIZAR
// incialmente debemos importar o dibujar una geometria del area donde se encuentran los predios viables para interpretación 
// esta se puede hacer con la herramienta de dibujo de geometrias o importando un shapefile o archivo CSV esto con el fin de deli-
// mitar la zona de estudio, una vez hecho esto procedemos a colocarle nombre en la zona de archivos importados.
Map.addLayer({
  eeObject:EAM, 
  name: 'EAM', 
  shown: 0, 
});

Map.addLayer({
  eeObject:ZH, 
  name: 'Zonas Humedas', 
  shown: 0, 
});

//throw('stop'); // ---> PONER EN COMENTARIOS CON: //, EN EL MOMENTO DE TENER LA UBICACIÓN DEL PREDIO PARA LA DESCARGA

//===============================================================================================================================
//                                MODELO DESCARGA DATOS NDWI PARA UNA REGIÓN DE INTERÉS
//===============================================================================================================================


//===============================================================================================================================//

//*******************************************************************************************************************************/
//                         PART 1: DEFINIR LOS PARÁMETROS DEL FILTRO DE COLECCIÓN Y DE LA MÁSCARA DE NUBE
// Parameter	Type	Description
/*
AOI             ==  ee.Geometry	Area of interest
START_DATE      ==  string	Fecha de inicio de la recogida de imágenes (inclusive)
END_DATE        ==  string	Fecha de finalización de la recogida de imágenes (excluyente)
CLOUD_FILTER    ==  integer Porcentaje máximo de cobertura de nubes permitido en la colección de imágenes
CLD_PRB_THRESH	==  integer Probabilidad de nubes (%); los valores superiores se consideran nubes
NIR_DRK_THRESH	==  float Reflectancia en el infrarrojo cercano; los valores inferiores se consideran sombra potencial de nubes
CLD_PRJ_DIST    ==  float Maximum Distancia máxima (km) para buscar sombras de nubes desde los bordes de las nubes
BUFFER          ==  integer Distancia (m) para dilatar el borde de los objetos identificados como nubes
*/


var AOI = ee.FeatureCollection(ROI);
var START_DATE = '2020-01-01';
var END_DATE = '2020-12-30';
var CLOUD_FILTER = 50;
var CLD_PRB_THRESH = 50;
var NIR_DRK_THRESH = 0.15;
var CLD_PRJ_DIST = 1;
var BUFFER = 50;


//Center the map view on a given AOI with the zool level.
Map.centerObject(AOI); 

//**********************************************************************************************************************************
//                                       PART 2: ESTABLECIMIENTO DE LAS FUNCIÓNES

/* 
Defina una función para filtrar las colecciones SR y s2cloudless según los parámetros de área de interés y fecha 
y, a continuación, únalas en la propiedad system:index. El resultado es una copia de la colección SR en la que cada 
imagen tiene una nueva propiedad 's2cloudless' cuyo valor es la imagen s2cloudless correspondiente. 
*/

function get_s2_sr_cld_col(aoi, start_date, end_date) {
    // # Import and filter S2 SR.
    var s2_sr_col = (ee.ImageCollection('COPERNICUS/S2_SR')
        .filterBounds(aoi)
        .filterDate(start_date, end_date)
        .filter(ee.Filter.lte('CLOUDY_PIXEL_PERCENTAGE', CLOUD_FILTER)))
        .map(function(image1){return image1.clip(AOI)});

    // # Import and filter s2cloudless.
    var s2_cloudless_col = (ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY')
        .filterBounds(aoi)
        .filterDate(start_date, end_date))
        .map(function(image1){return image1.clip(AOI)});
        
    // # Join the filtered s2cloudless collection to the SR collection by the 'system:index' property.
    return ee.ImageCollection(ee.Join.saveFirst('s2cloudless').apply({
        'primary': s2_sr_col,
        'secondary': s2_cloudless_col,
        'condition': ee.Filter.equals({
            'leftField': 'system:index',
            'rightField': 'system:index'
        })
    }));
}


function add_cloud_bands(img) {
    // # Get s2cloudless image, subset the probability band.
    var cld_prb = ee.Image(img.get('s2cloudless')).select('probability');

    // # Condition s2cloudless by the probability threshold value.
    var is_cloud = cld_prb.gt(CLD_PRB_THRESH).rename('clouds');

    // # Add the cloud probability layer and cloud mask as image bands.
    return img.addBands(ee.Image([cld_prb, is_cloud]));
    
}


function add_shadow_bands(img) {
    // # Identify water pixels from the SCL band.
    var not_water = img.select('SCL').neq(6);

    // # Identify dark NIR pixels that are not water (potential cloud shadow pixels).
    var SR_BAND_SCALE = 1e4;
    var dark_pixels = img.select('B8').lt(NIR_DRK_THRESH*SR_BAND_SCALE).multiply(not_water).rename('dark_pixels');

    // # Determine the direction to project cloud shadow from clouds (assumes UTM projection).
    var shadow_azimuth = ee.Number(90).subtract(ee.Number(img.get('MEAN_SOLAR_AZIMUTH_ANGLE')));

    // # Project shadows from clouds for the distance specified by the CLD_PRJ_DIST input.
    var cld_proj = (img.select('clouds').directionalDistanceTransform(shadow_azimuth, CLD_PRJ_DIST*10)
        .reproject({'crs': img.select(0).projection(), 'scale': 100})
        .select('distance')
        .mask()
        .rename('cloud_transform'));

    // # Identify the intersection of dark pixels with cloud shadow projection.
    var shadows = cld_proj.multiply(dark_pixels).rename('shadows');

    // # Add dark pixels, cloud projection, and identified shadows as image bands.
    return img.addBands(ee.Image([dark_pixels, cld_proj, shadows]));
}


function add_cld_shdw_mask(img) {
    // # Add cloud component bands.
    var img_cloud = add_cloud_bands(img);

    // # Add cloud shadow component bands.
    var img_cloud_shadow = add_shadow_bands(img_cloud);

    // # Combine cloud and shadow mask, set cloud and shadow as value 1, else 0.
    var is_cld_shdw = img_cloud_shadow.select('clouds').add(img_cloud_shadow.select('shadows')).gt(0);

    // # Remove small cloud-shadow patches and dilate remaining pixels by BUFFER input.
    // # 20 m scale is for speed, and assumes clouds don't require 10 m precision.
    is_cld_shdw = (is_cld_shdw.focal_min(2).focal_max(BUFFER*2/20)
        .reproject({'crs': img.select([0]).projection(), 'scale': 20})
        .rename('cloudmask'));

    // # Add the final cloud-shadow mask to the image.
    return img_cloud_shadow.addBands(is_cld_shdw);
}



function apply_cld_shdw_mask(img) {
    // # Subset the cloudmask band and invert it so clouds/shadow are 0, else 1.
    var not_cld_shdw = img.select('cloudmask').not();

    // # Subset reflectance bands and update their masks, return the result.
    return img.select('B.*').updateMask(not_cld_shdw);
}

var s2_sr_cld_col = get_s2_sr_cld_col(AOI, START_DATE, END_DATE);

var s2_sr_median = (s2_sr_cld_col.map(add_cld_shdw_mask)
                             .map(apply_cld_shdw_mask)
                             .mosaic());

Map.addLayer(s2_sr_median, {bands: ["B4","B3","B2"],
gamma: 1,
max: 2744.2910152719487,
min: 459.9036445041072,
opacity: 1}, 'S2 cloud-free mosaic');

/******************************************************************************************************************************/
//                                     PART 6: AGREGAR INDICES ESPECTRALES

//Computes the normalized difference between
var addIndices = function(image) {
  var NDWI = image.normalizedDifference(['B3', 'B8']).rename(['NDWI']);
  return image.addBands(NDWI)};

var composite = addIndices(s2_sr_median);
var image = composite.clip(AOI);

Map.addLayer(image, {min: -0.2648029338212383, max: 0.42513437420015376, bands: ['NDWI']},'NDWI', 0);


// Reclasificar el NDWI usando una función de mapeo condicional
var NDWIReclass = image.select('NDWI').expression(
  "b(0) > 0.6 && b(0) <= 1 ? 3 : " +      // Superficie del agua,
  "b(0) > 0 && b(0) <= 0.6 ? 2 : " +    // Inundación, humedad,
  "b(0) > -0.3 && b(0) <= 0 ? 1 : " +    // Sequía moderada, superficies sin agua,
  "b(0) > -1 && b(0) <= -0.3 ? 0 : 0"   // Sequía, superficies sin agua
).rename('NDWI_Reclass').clip(AOI);

Map.addLayer(NDWIReclass, {min: 0, max: 3, palette: ['FFFF00', 'F2FF00', 'FEB607', '0086FF']},'NDWIReclass', 0);

var NDWIReclassMasked = NDWIReclass.updateMask(NDWIReclass.gt(1));

Map.addLayer(NDWIReclassMasked, {min: 1, max: 3, palette: ['FFFF00', 'F2FF00', 'FEB607', '0086FF']}, 'NDWIReclassMasked', 1);

Export.image.toDrive({
  image: NDWIReclassMasked, 
  description: 'NDWI_EAM_2018',
  folder: 'DESCARGAS_GEE',                                                              // Cambiar folder
  fileNamePrefix: 'NDWI_EAM_2018',                                                   // Cambiar nombre predio
  region: AOI, 
  scale: 30, 
  maxPixels: 1e11

});




