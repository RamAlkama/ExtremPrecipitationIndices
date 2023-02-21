/* CHIRPS

Precipitation extrem index  
----------------------------
CDD: Consecutive Dry Days 
CWD: Consecutive Wet Days 
PRCPTOT: annual total precip  (mm/year)
R10: number of heavy precipitation days (nbr days with P>10mm/day)
R20: number of heavy precipitation days (nbr days with P>20mm/day)
R95p: Annual total PRCP when RR > 95th percentile (mm/year)
R99p: Annual total PRCP when RR > 99th percentile (mm/year)
SDII:	Simple daily intensity index	Annual total precipitation divided by the number of wet days (mm/day)
Pmax: daily max precip (mm/day)
Pmax5 : maximum 5 days precipitation (mm/5days)
*/
//-----------------
// Functions

function drySpells(img, list){
  // get previous image
  var prev = ee.Image(ee.List(list).get(-1));
  // find areas gt precipitation threshold (gt==0, lt==1)
  var dry = img.select('precipitation').lt(precipThresh);
  // find areas  precipitation threshold (>1)
  var wet = img.select('precipitation').gt(precipThresh);
  // add previous day counter to today's counter
  var accum = prev.select('counter').add(dry).rename('counter');
  // add previous day counter to today's counter
  var accumWet = prev.select('counterw').add(wet).rename('counterw');
  // create a result image for iteration
  // precip < thresh will equal the accumulation of counters
  // otherwise it will equal zero
  var out = img.select('precipitation').addBands(img.select('counter').where(dry.eq(1),accum)).uint8()
                                       .addBands(img.select('counterw').where(wet.eq(1),accumWet)).uint8();
  return ee.List(list).add(out);
}


// Define function to mask image by quality metric.
function qualityMask10(img) {
  // Make boolean image where high quality pixels are value 1 and all else 0.
  var mask = img.select('precipitation').gt(10);
  // Update the image mask and return it.
  return img.updateMask(mask);
}

function qualityMask20(img) { var mask = img.select('precipitation').gt(20);
                              return img.updateMask(mask);
                             }
                             
function qualityMask1(img) { var mask = img.select('precipitation').gt(1);
                              return img.updateMask(mask);
                             }

// Define a function to calculate the number of valid observations for a given
// pixel's time series. Takes a collection of masked images.
function countValidPixels(collection) {
  // For each image in the collection return the mask; returns an image
  // collection.
  return collection.map(function(img) {
    return img.select(0).mask();
  })
  // Sum the masks; this gives valid pixel count.
  .sum()
  // Optionally mask pixels that have 0 observation over the give time series.
  //.selfMask();
}



function countPercentile95(img) {var mask = img.select('precipitation').gt(precipitation_p95);
                               var prec = img.select('precipitation')
                               return prec.updateMask(mask);
                               }

function countPercentile99(img) {var mask = img.select('precipitation').gt(precipitation_p99);
                               var prec = img.select('precipitation')
                               return prec.updateMask(mask);
                               }  
// ----------------
// load Algeria boundaries
var countries = ee.FeatureCollection("USDOS/LSIB_SIMPLE/2017");
var algerie =countries.filter(ee.Filter.eq('country_co','AG'));
// DATA
var collection = ee.ImageCollection("UCSB-CHG/CHIRPS/DAILY");
print (collection.select("precipitation").filterDate('2020-01-01', '2020-01-06'))
//get projection of the CHIRPS dataset
var proj = ee.Image(collection.first()).projection();
// first mask days <1mm of prec then compute 95 and 99 percentiles
var p95p99=collection.select("precipitation").filterDate('1981-01-01', '2021-01-01').map(qualityMask1).reduce(ee.Reducer.percentile([95,99]));
print (p95p99)
var precipitation_p95=p95p99.select('precipitation_p95');
var precipitation_p99=p95p99.select('precipitation_p99');

for(var yr=1981; yr<2021; yr++) {
// Define time range
var startyear = yr;
var endyear = yr+1;

var startmonth = 1;
var endmonth = 12;

// Set date in ee date format
var startdate = ee.Date.fromYMD(startyear,1,1);
var enddate = ee.Date.fromYMD(endyear,1,1);

// Filter data
var datain_t = collection.filterDate(startdate, enddate)
  .filter(ee.Filter.calendarRange(startmonth,endmonth, 'month'))
  .select("precipitation").map(function(img){
     return img.addBands(ee.Image(0).reproject(proj).uint8().rename('counter'))
               .addBands(ee.Image(0).reproject(proj).uint8().rename('counterw'));
  })
  .sort('system:time_start');


// DIFFERENT SCALE, I NEED COUNTER HAS THE SAME AS PRECIPITATION (5565.974539663679)
var scale_p = ee.Image(datain_t.first().select('precipitation')).projection().nominalScale()
//print(scale_p,'scale_p')

var scale_c = ee.Image(datain_t.first().select('counter')).projection().nominalScale()
//print(scale_c,'scale_c')

// // START 
var dataset = datain_t
.filterDate(startdate,enddate)
.sort('system:time_start:');
//print(dataset,"dataset");

var precipThresh = 1; // mm

// create first image for iteration
var first = ee.List([ee.Image(dataset.first())]);

// apply dry speall iteration function
var CD = ee.ImageCollection.fromImages(
    dataset.iterate(drySpells,first)
).max(); // get the max value


var col=collection.select("precipitation").filterDate(startdate, enddate);
var Pmax=col.max()
var CDD= CD.select('counter').rename('CDD');
var CWD= CD.select('counterw').rename('CWD');
var PRCPTOT=col.sum().rename('PRCPTOT');

// Apply quality mask to all images in the collection.
var VegCol = col.map(qualityMask10);
// Apply the valid pixel counting function.
var R10 = countValidPixels(VegCol).rename('R10');

// Do same for R20
var VegCol = col.map(qualityMask20);
var R20 = countValidPixels(VegCol).rename('R20');
// Do same for R1
var VegCol = col.map(qualityMask1);
var R1 = countValidPixels(VegCol).rename('R1');

var R95p=col.map(countPercentile95).sum().rename('R95p');
var R99p=col.map(countPercentile99).sum().rename('R99p');
var SDII=PRCPTOT.divide(R1).rename('SDII')

// ---- maximum sum of consicutive 5 days  ----
// this is the window size in days
var dwindow = 5;
// just calculating number of windows so that i can map over it
// i could go for iterate with a break condition but i prefer map
// as i can compute parallelly 
var numberOfWindows = enddate.difference(startdate,'day').subtract(dwindow-1).toInt();
// generating a sequence that can be used as an indicator for my window step
var sequence = ee.List.sequence(0, numberOfWindows); // inclusive series so the number of windows will be correct

// mapping over the sequence
sequence = sequence.map(function(num){
  // just casting element of sequence to a number object
  num = ee.Number(num);
  // finding the start and end point of my windows
  var windowStart = startdate.advance(num, 'day');
  var windowEnd = startdate.advance(num.add(dwindow), 'day');
  // selecting images that fall within those windows
  var subset = col.filterDate(windowStart,windowEnd);
  // calculating the sum prec of that window
  return subset.sum().set('system:time_start',windowStart);
});

// converting list of sum images to imagecollection 
var composites = ee.ImageCollection.fromImages(sequence);

// calculating the max image among those 5 days composites
var Pmax5 = composites.max();

//

var Indices=CDD.addBands(CWD).addBands(PRCPTOT).addBands(R10).addBands(R20)
               .addBands(R95p).addBands(R99p).addBands(SDII).addBands(Pmax).addBands(Pmax5);

Export.image.toDrive({
  image: Indices.clip(algerie).float(),
description: 'CHIRPS_PRCPindex_'+yr, 
folder: 'CHIRPS',
  crs: 'EPSG:4326',
  region: ee.Geometry.Polygon([-9, 38, 0, 38, 12, 38, 12, 19, 0, 19, -9, 19], null, false),
  maxPixels: 1e13,
  scale:'5565.974539663679'
  }); 
}
// display results
Map.setCenter(3.53,29, 5);
Map.addLayer(CDD.clip(algerie),{min:0,max:90,palette:'#9ecae1,#ffffff,#ffeda0,#feb24c,#f03b20'},'CDD');
Map.addLayer(CWD.clip(algerie),{min:0,max:10,palette:'#9ecae1,#ffffff,#ffeda0,#feb24c,#f03b20'},'CWD');
Map.addLayer(PRCPTOT.clip(algerie),{min:0,max:700,palette:'#9ecae1,#ffffff,#ffeda0,#feb24c,#f03b20'},'PRCPTOT');
Map.addLayer(R10.clip(algerie),{min:0,max:20,palette:'#9ecae1,#ffffff,#ffeda0,#feb24c,#f03b20'},'R10');
Map.addLayer(R20.clip(algerie),{min:0,max:10,palette:'#9ecae1,#ffffff,#ffeda0,#feb24c,#f03b20'},'R20');
Map.addLayer(R95p.clip(algerie),{min:0,max:200,palette:'#9ecae1,#ffffff,#ffeda0,#feb24c,#f03b20'},'R95p');
Map.addLayer(R99p.clip(algerie),{min:0,max:100,palette:'#9ecae1,#ffffff,#ffeda0,#feb24c,#f03b20'},'R99p');
Map.addLayer(SDII.clip(algerie),{min:0,max:20,palette:'#9ecae1,#ffffff,#ffeda0,#feb24c,#f03b20'},'SDII');

