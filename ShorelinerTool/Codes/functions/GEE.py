import ee
import os
from urllib.request import urlretrieve
import shutil
import zipfile
from osgeo import gdal
global inputs
#%% General

def index_f2(name):
    if name == 'SCoWI':
        return([SCoWIL5,SCoWIL7,SCoWIL8,SCoWIL9,SCoWIS2])
    elif name == 'MNDWI':
        return([MNDWIL5,MNDWIL7,MNDWIL8,MNDWIL9,MNDWIS2])
    elif name == 'NDVI':
        return([NDVIL5,NDVIL7,NDVIL8,NDVIL9,NDVIS2])
    else:
        return('Error')

def index_f(name):
    if name == 'SCoWI':
        return(SCoWI)
    elif name == 'MNDWI':
        return(MNDWI)
    elif name == 'NDWI':
        return(NDWI)
    else:
        return('Error')

def resample(satlist,inputs):
    process = inputs['Preprocessing']
    if 'Bicubic' in process:
        satlist = satlist.map(bicubicResample)
    if 'Bilinear' in process:
        satlist = satlist.map(bilinearResample)
    return satlist

def preprocess(satlist,inputs):
    process = inputs['Preprocessing']
    if 'Pansharpening' in process:
        satlist = satlist.map(pansharpening)
    return satlist

#%% Preprocess
def bicubicResample(img):
    return img.resample('bicubic')

def bilinearResample(img):
    return img.resample('bilinear')

def pansharpening(image):
    hsv = image.select(['red', 'green', 'blue']).rgbToHsv();
    hsv2= image.select(['nir', 'green', 'blue']).rgbToHsv();
    hsv3 = image.select(['swir1', 'green', 'blue']).rgbToHsv();
    hsv4 = image.select(['swir2', 'green', 'blue']).rgbToHsv();
    
    sharpened = ee.Image.cat([
            hsv.select('hue'), hsv.select('saturation'), image.select('panchromatic')
            ]).hsvToRgb();
    sharpened2 = ee.Image.cat([
            hsv2.select('hue'), hsv2.select('saturation'), image.select('panchromatic')
            ]).hsvToRgb();
    sharpened3 = ee.Image.cat([
            hsv3.select('hue'), hsv3.select('saturation'), image.select('panchromatic')
            ]).hsvToRgb();
    sharpened4 = ee.Image.cat([
            hsv4.select('hue'), hsv4.select('saturation'), image.select('panchromatic')
            ]).hsvToRgb();

    blue = sharpened.select(['blue'],['blue'])
    green = sharpened.select(['green'],['green'])
    red = sharpened.select(['red'],['red'])
    nir = sharpened2.select(['red'],['nir'])
    swir1 = sharpened3.select(['red'],['swir1'])
    swir2 = sharpened4.select(['red'],['swir2'])
    
    newImage = ee.Image.cat([red,green,blue,nir,swir1,swir2])
    return newImage.copyProperties(image,['system:time_start','satellite','SCENE_CENTER_TIME'])

def normalizeBandNamesS2(image):
    newImage = image.select(['B2','B3','B4','B8','B11','B12'],['blue','green','red','nir','swir1','swir2'])
    return newImage.set({'satellite':'S2', 'date':image.date()}).setDefaultProjection('EPSG:3857').reproject(
        **{
      'crs':'EPSG:3857',
      'scale': 10}).float().copyProperties(image,['SCENE_CENTER_TIME'])
    
def normalizeBandNamesL9(image):
    newImage = image.select(['B2','B3','B4','B5','B6','B7','B8'],['blue','green','red','nir','swir1','swir2','panchromatic'])
    return newImage.set({'satellite':'L9', 'date':image.date()}).setDefaultProjection('EPSG:3857').reproject(
        **{
      'crs':'EPSG:3857',
      'scale': 15}).float().copyProperties(image,['SCENE_CENTER_TIME'])

def normalizeBandNamesL8(image):
    newImage = image.select(['B2','B3','B4','B5','B6','B7','B8'],['blue','green','red','nir','swir1','swir2','panchromatic'])
    return newImage.set({'satellite':'L8', 'date':image.date()}).setDefaultProjection('EPSG:3857').reproject(
        **{
      'crs':'EPSG:3857',
      'scale': 15}).float().copyProperties(image,['system:time_start','SCENE_CENTER_TIME'])

def normalizeBandNamesL7(image):
    newImage = image.select(['B1','B2','B3','B4','B5','B7','B8'],['blue','green','red','nir','swir1','swir2','panchromatic'])
    return newImage.set({'satellite':'L7', 'date':image.date()}).setDefaultProjection('EPSG:3857').reproject(
        **{
      'crs':'EPSG:3857',
      'scale': 15}).float().copyProperties(image,['SCENE_CENTER_TIME'])

def normalizeBandNamesL5(image):
    newImage = image.select(['B1','B2','B3','B4','B5','B7'],['blue','green','red','nir','swir1','swir2'])
    return newImage.set({'satellite':'L5', 'date':image.date()}).setDefaultProjection('EPSG:3857').reproject(
        **{
      'crs':'EPSG:3857',
      'scale': 30}).float().copyProperties(image,['SCENE_CENTER_TIME'])

def maskLandsat(image):
  # Bits 3 and 5 are cloud shadow and cloud, respectively.
  cloudShadowBitMask = 1 << 3;
  cloudsBitMask = 1 << 5;
  # Get the pixel QA band.
  qa = image.select('QA_PIXEL');
  #  Both flags should be set to zero, indicating clear conditions.
  maskShadow = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
  maskCloud = qa.bitwiseAnd(cloudsBitMask).eq(0);
  #Return the masked image, scaled to reflectance, without the QA bands.
  return image.updateMask(maskShadow).updateMask(maskCloud).copyProperties(image, ['SCENE_CENTER_TIME'])

def get_s2_sr_cld_col(aoi, start_date, end_date,cloud_cover_treshold = 50):
    # Import and filter S2 SR.
    s2_sr_col = (ee.ImageCollection('COPERNICUS/S2_HARMONIZED')
        .filterBounds(aoi)
        .filterDate(start_date, end_date)
        .filter(ee.Filter.lte('CLOUDY_PIXEL_PERCENTAGE', cloud_cover_treshold)))

    # Import and filter s2cloudless.
    s2_cloudless_col = (ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY')
        .filterBounds(aoi)
        .filterDate(start_date, end_date))

    # Join the filtered s2cloudless collection to the SR collection by the 'system:index' property.
    return s2_sr_col.combine(s2_cloudless_col)


def add_cloud_maskS2(img):
    # Get s2cloudless image, subset the probability band.
    cld_prb = img.select('probability')

    # Condition s2cloudless by the probability threshold value.
    is_cloud = cld_prb.gt(75).reproject(**{'crs': 'EPSG:3857', 'scale': 20})
    
    no_cloud_img = is_cloud.select('probability').Not()
    # Add the cloud probability layer and cloud mask as image bands.
    return img.select('B.*').updateMask(no_cloud_img).addBands(ee.Image([cld_prb, is_cloud]))
#%% Indexes
def SCoWI(image):
    return image.expression('(BLUE + 2 * GREEN - 2 * NIR - 0.75 * SWIR1 - 0.5 * SWIR2)',#
    {
      'BLUE': image.select('blue'),
      'GREEN': image.select('green'),
      'NIR': image.select('nir'),
      'SWIR1': image.select('swir1'),
      'SWIR2': image.select('swir2')
    }).rename('SCoWI').copyProperties(image,['system:time_start','satellite','SCENE_CENTER_TIME'])

def MNDWI(image):
    return image.expression('(SWIR1-GREEN)/(SWIR1+GREEN)',#
    {
      'GREEN': image.select('green'),
      'SWIR1': image.select('swir1')
    }).rename('MNDWI').copyProperties(image,['system:time_start','satellite','SCENE_CENTER_TIME'])

def NDVI(image):
    return image.expression('(NIR-RED)/(NIR+GREEN)',#
    {
      'RED': image.select('red'),
      'NIR': image.select('nir')
    }).rename('NDVI').copyProperties(image,['system:time_start','satellite','SCENE_CENTER_TIME'])

def SBI(image):
    return image.expression('BLUE+2*GREEN-1.75*NIR',#
    {
      'BLUE': image.select('blue'),
      'GREEN': image.select('green'),
      'NIR': image.select('nir')
    }).rename('SBI').copyProperties(image,['system:time_start','satellite','SCENE_CENTER_TIME'])


def NDWI(image):
    return image.expression('(NIR-GREEN)/(NIR+GREEN)',#
    {
      'GREEN': image.select('green'),
      'NIR': image.select('nir')
    }).rename('NDWI').copyProperties(image,['system:time_start','satellite','SCENE_CENTER_TIME'])


def RGB(image):
        #image = image.resample('bicubic')
        
        red = image.select(['red'])
        blue = image.select(['blue'])
        green = image.select(['green'])
        
        new_image = ee.Image.cat([red,green,blue])
        return new_image.copyProperties(image,['system:time_start','satellite','SCENE_CENTER_TIME'])
#%% Data/Folder 

def folderCreation(inputs):
    
    if inputs['Waterline']:
        try:
            os.mkdir(inputs['WaterlineIndex'])
        except:
            pass
    if inputs['SandBar']:
        try:
            os.mkdir(inputs['SandBarIndex'])
        except:
            pass
    if inputs['Vegetation']:
        try:
            os.mkdir(inputs['VegetationIndex'])
        except:
            pass

def saveImage(imggee, i, filepath, index, poly):
    
    img = imggee.getInfo()
    if img['properties']['satellite']=='S2':
        name=img['properties']['system:index'][-22:-7]+'_'+img['properties']['satellite']
    else:
        name=(img['properties']['system:index'][-8:]+'T'+img['properties']['SCENE_CENTER_TIME'][:2]+
            img['properties']['SCENE_CENTER_TIME'][3:5] +img['properties']['SCENE_CENTER_TIME'][6:8]+
            '_'+img['properties']['satellite'])
    print(name)
    cond1=not(True in [name in j for j in os.listdir(filepath)])
    if (cond1):
    
        print('Downloading...')
        url = ee.data.makeDownloadUrl(ee.data.getDownloadId({
            'image': imggee,
            'bands': index,
            'region': poly,
            'name': name
            }))
    
        try:
            local_zip, headers = urlretrieve(url)
            # move zipfile from temp folder to data folder
            dest_file = os.path.join(filepath, 'imagezip')
            shutil.move(local_zip,dest_file)
            # unzip file
            with zipfile.ZipFile(dest_file) as local_zipfile:
                for fn in local_zipfile.namelist():
                    local_zipfile.extract(fn, filepath)
            # filepath + filename to single bands
                fn_tifs = [os.path.join(filepath,_) for _ in local_zipfile.namelist()]
        # stack bands into single .tif
            outds = gdal.BuildVRT(os.path.join(filepath,'stacked'+str(i)+'.vrt'), fn_tifs, separate=True)
            outds = gdal.Translate(os.path.join(filepath,name+'.tif'), outds)
        # delete single-band files
            for fn in fn_tifs: os.remove(fn)
        # delete .vrt file
            os.remove(os.path.join(filepath,'stacked'+str(i)+'.vrt'))
        # delete zipfile
            os.remove(dest_file)
        # delete data.tif.aux (not sure why this is created)
            if os.path.exists(os.path.join(filepath,'data.tif.aux')):
                os.remove(os.path.join(filepath,'data.tif.aux'))
        except:
            pass

#%%SCoWI

def SCoWIS2(image):
        #image = image.resample('bicubic')
        return image.expression('(BLUE + 2 * GREEN - 2 * NIR - 0.75 * SWIR1 - 0.5 * SWIR2)',#'/(B2 + 2 * B3 + 2 * B8 + 0.75 * B11 + 0.5 * B12)',
        {
            'BLUE': image.select('B2'),
            'GREEN': image.select('B3'),
            'NIR': image.select('B8'),
            'SWIR1': image.select('B11').resample('bicubic'),
            'SWIR2': image.select('B12').resample('bicubic')
  }).rename('SCoWI').set({'satellite':'S2', 'date':image.date()}).setDefaultProjection('EPSG:3857').reproject(
      **{
    'crs':'EPSG:3857',
    'scale': 10}
      ).float().copyProperties(image,['SCENE_CENTER_TIME'])

def SCoWIL9(image):
        hsv = image.resample('bicubic').select(['B5', 'B3', 'B2']).rgbToHsv();
        hsv2= image.resample('bicubic').select(['B6', 'B3', 'B2']).rgbToHsv();
        hsv3 = image.resample('bicubic').select(['B7', 'B3', 'B2']).rgbToHsv();
        
        sharpened = ee.Image.cat([
                hsv.select('hue'), hsv.select('saturation'), image.select('B8')
                ]).hsvToRgb();
        sharpened2 = ee.Image.cat([
                hsv2.select('hue'), hsv2.select('saturation'), image.select('B8')
                ]).hsvToRgb();
        sharpened3 = ee.Image.cat([
                hsv3.select('hue'), hsv3.select('saturation'), image.select('B8')
                ]).hsvToRgb();
        
        #return sharpened
        newImage = image.select(['B2'],['B2'])
        


        blue = sharpened.select(['blue'],['blue'])
        newImage = newImage.addBands(blue)

        green = sharpened.select(['green'],['green'])
        newImage = newImage.addBands(green)

        nir = sharpened.select(['red'],['nir'])
        newImage = newImage.addBands(nir)

        swir1 = sharpened2.select(['red'],['swir1'])
        newImage = newImage.addBands(swir1)

        swir2 = sharpened3.select(['red'],['swir2'])
        newImage = newImage.addBands(swir2)
            
        return newImage.expression('(BLUE + 2 * GREEN - 2 * NIR - 0.75 * SWIR1 - 0.5 * SWIR2)',#
        {
          'BLUE': newImage.select('blue'),
          'GREEN': newImage.select('green'),
          'NIR': newImage.select('nir'),
          'SWIR1': newImage.select('swir1'),
          'SWIR2': newImage.select('swir2')
        }).rename('SCoWI').set({'satellite':'L9', 'date':image.date()}).setDefaultProjection('EPSG:3857').reproject(
              **{
            'crs':'EPSG:3857',
            'scale': 15}
              ).float().copyProperties(image,['SCENE_CENTER_TIME'])

def SCoWIL7(image):
        hsv = image.resample('bicubic').select(['B4', 'B2', 'B1']).rgbToHsv();
        hsv2= image.resample('bicubic').select(['B5', 'B2', 'B1']).rgbToHsv();
        hsv3 = image.resample('bicubic').select(['B7', 'B2', 'B1']).rgbToHsv();
        
        sharpened = ee.Image.cat([
                hsv.select('hue'), hsv.select('saturation'), image.select('B8')
                ]).hsvToRgb();
        sharpened2 = ee.Image.cat([
                hsv2.select('hue'), hsv2.select('saturation'), image.select('B8')
                ]).hsvToRgb();
        sharpened3 = ee.Image.cat([
                hsv3.select('hue'), hsv3.select('saturation'), image.select('B8')
                ]).hsvToRgb();
        
        #return sharpened
        newImage = image.select(['B2'],['B2'])
        


        blue = sharpened.select(['blue'],['blue'])
        newImage = newImage.addBands(blue)

        green = sharpened.select(['green'],['green'])
        newImage = newImage.addBands(green)

        nir = sharpened.select(['red'],['nir'])
        newImage = newImage.addBands(nir)

        swir1 = sharpened2.select(['red'],['swir1'])
        newImage = newImage.addBands(swir1)

        swir2 = sharpened3.select(['red'],['swir2'])
        newImage = newImage.addBands(swir2)
        
        
        return newImage.expression('(BLUE + 2 * GREEN - 2 * NIR - 0.75 * SWIR1 - 0.5 * SWIR2)',#
        {
          'BLUE': newImage.select('blue'),
          'GREEN': newImage.select('green'),
          'NIR': newImage.select('nir'),
          'SWIR1': newImage.select('swir1'),
          'SWIR2': newImage.select('swir2')
        }).rename('SCoWI').set({'satellite':'L7', 'date':image.date()}).setDefaultProjection('EPSG:3857').reproject(
              **{
            'crs':'EPSG:3857',
            'scale': 15}
              ).float().copyProperties(image,['SCENE_CENTER_TIME'])
            
def SCoWIL5(image):
        image = image.resample('bicubic')
        return image.expression('(BLUE + 2 * GREEN - 2 * NIR - 0.75 * SWIR1 - 0.5 * SWIR2)',#
        {
          'BLUE': image.select('B1'),
          'GREEN': image.select('B2'),
          'NIR': image.select('B4'),
          'SWIR1': image.select('B5'),
          'SWIR2': image.select('B7')
        }).rename('SCoWI').set({'satellite':'L5', 'date':image.date()}).setDefaultProjection('EPSG:3857').reproject(
             **{
           'crs':'EPSG:3857',
           'scale': 30
             }).float().copyProperties(image,['SCENE_CENTER_TIME'])

def SCoWIL8(image):
        hsv = image.resample('bicubic').select(['B5', 'B3', 'B2']).rgbToHsv();
        hsv2= image.resample('bicubic').select(['B6', 'B3', 'B2']).rgbToHsv();
        hsv3 = image.resample('bicubic').select(['B7', 'B3', 'B2']).rgbToHsv();
        
        sharpened = ee.Image.cat([
                hsv.select('hue'), hsv.select('saturation'), image.select('B8')
                ]).hsvToRgb();
        sharpened2 = ee.Image.cat([
                hsv2.select('hue'), hsv2.select('saturation'), image.select('B8')
                ]).hsvToRgb();
        sharpened3 = ee.Image.cat([
                hsv3.select('hue'), hsv3.select('saturation'), image.select('B8')
                ]).hsvToRgb();
        
        #return hsv
        newImage = image.select(['B2'],['B2'])
        


        blue = sharpened.select(['blue'],['blue'])
        newImage = newImage.addBands(blue)

        green = sharpened.select(['green'],['green'])
        newImage = newImage.addBands(green)

        nir = sharpened.select(['red'],['nir'])
        newImage = newImage.addBands(nir)

        swir1 = sharpened2.select(['red'],['swir1'])
        newImage = newImage.addBands(swir1)

        swir2 = sharpened3.select(['red'],['swir2'])
        newImage = newImage.addBands(swir2)
        
        
        return ee.Image(newImage.expression('(BLUE + 2 * GREEN - 2 * NIR - 0.75 * SWIR1 - 0.5 * SWIR2)',#'/(B2 + 2 * B3 + 2 * B8 + 0.75 * B11 + 0.5 * B12)',
        {
          'BLUE': newImage.select('blue'),
          'GREEN': newImage.select('green'),
          'NIR': newImage.select('nir'),
          'SWIR1': newImage.select('swir1'),
          'SWIR2': newImage.select('swir2')
        }).rename('SCoWI').set({'satellite':'L8', 'date':image.date()}).setDefaultProjection('EPSG:3857').reproject(
              **{
            'crs':'EPSG:3857',
            'scale': 15}
              ).float().copyProperties(image,['SCENE_CENTER_TIME']))

#%%MNDWI

def MNDWIS2(image):
        return image.expression('(SWIR1 - GREEN)/(SWIR1 + GREEN)',#'/(B2 + 2 * B3 + 2 * B8 + 0.75 * B11 + 0.5 * B12)',
        {
            'GREEN': image.select('B3'),
            'SWIR1': image.select('B11'),
  }).rename('MNDWI').set({'satellite':'S2', 'date':image.date()}).setDefaultProjection('EPSG:32610').reproject(
      **{
    'crs': 'EPSG:3857',
    'scale': 20}
      ).float().copyProperties(image,['SCENE_CENTER_TIME'])

def MNDWIL9(image):
        return image.expression('(SWIR1 - GREEN)/(SWIR1 + GREEN)',#
        {
          'GREEN': image.select('B3'),
          'SWIR1': image.select('B6'),
        }).rename('MNDWI').set({'satellite':'L9', 'date':image.date()}).setDefaultProjection('EPSG:32610').reproject(
             **{
           'crs': 'EPSG:3857',
           'scale': 30}
             ).float().copyProperties(image,['SCENE_CENTER_TIME'])
def MNDWIL7(image):
        return image.expression('(SWIR1 - GREEN)/(SWIR1 + GREEN)',#
        {
          'GREEN': image.select('B2'),
          'SWIR1': image.select('B5'),
        }).rename('MNDWI').set({'satellite':'L7', 'date':image.date()}).setDefaultProjection('EPSG:32610').reproject(
             **{
           'crs': 'EPSG:3857',
           'scale': 30}
             ).float().copyProperties(image,['SCENE_CENTER_TIME'])
            
def MNDWIL5(image):
        return image.expression('(SWIR1 - GREEN)/(SWIR1 + GREEN)',#
        {
          'GREEN': image.select('B2'),
          'SWIR1': image.select('B5'),
        }).rename('MNDWI').set({'satellite':'L5', 'date':image.date()}).setDefaultProjection('EPSG:32610').reproject(
             **{
           'crs': 'EPSG:3857',
           'scale': 30}
             ).float().copyProperties(image,['SCENE_CENTER_TIME'])
def MNDWIL8(image):
       return image.expression('(SWIR1 - GREEN)/(SWIR1 + GREEN)',#'/(B2 + 2 * B3 + 2 * B8 + 0.75 * B11 + 0.5 * B12)',
        {
          'GREEN': image.select('B3'),
          'SWIR1': image.select('B6'),
        }).rename('MNDWI').set({'satellite':'L8', 'date':image.date()}).setDefaultProjection('EPSG:32610').reproject(
             **{
           'crs': 'EPSG:3857',
           'scale': 30}
             ).float().copyProperties(image,['SCENE_CENTER_TIME'])


#%%NDVI
def NDVIS2(image):
        return image.expression('(NIR - RED)/(NIR + RED)',#'/(B2 + 2 * B3 + 2 * B8 + 0.75 * B11 + 0.5 * B12)',
        {
            'NIR': image.select('B8'),
            'RED': image.select('B4'),
  }).rename('NDVI').set({'satellite':'S2', 'date':image.date()}).setDefaultProjection('EPSG:32610').reproject(
      **{
    'crs': 'EPSG:3857',
    'scale': 20}
      ).float().copyProperties(image,['SCENE_CENTER_TIME'])

def NDVIL9(image):
        return image.expression('(NIR - RED)/(NIR + RED)',#
        {
          'NIR': image.select('B5'),
          'RED': image.select('B4'),
        }).rename('NDVI').set({'satellite':'L9', 'date':image.date()}).setDefaultProjection('EPSG:32610').reproject(
             **{
           'crs': 'EPSG:3857',
           'scale': 30}
             ).float().copyProperties(image,['SCENE_CENTER_TIME'])
def NDVIL7(image):
        return image.expression('(NIR - RED)/(NIR + RED)',#
        {
          'NIR': image.select('B4'),
          'RED': image.select('B3'),
        }).rename('NDVI').set({'satellite':'L7', 'date':image.date()}).setDefaultProjection('EPSG:32610').reproject(
             **{
           'crs': 'EPSG:3857',
           'scale': 30}
             ).float().copyProperties(image,['SCENE_CENTER_TIME'])
            
def NDVIL5(image):
        return image.expression('(NIR - RED)/(NIR + RED)',#
        {
          'NIR': image.select('B4'),
          'RED': image.select('B3'),
        }).rename('NDVI').set({'satellite':'L5', 'date':image.date()}).setDefaultProjection('EPSG:32610').reproject(
             **{
           'crs': 'EPSG:3857',
           'scale': 30}
             ).float().copyProperties(image,['SCENE_CENTER_TIME'])
def NDVIL8(image):
       return image.expression('(NIR - RED)/(NIR + RED)',#'/(B2 + 2 * B3 + 2 * B8 + 0.75 * B11 + 0.5 * B12)',
        {
          'NIR': image.select('B5'),
          'RED': image.select('B4'),
        }).rename('NDVI').set({'satellite':'L8', 'date':image.date()}).setDefaultProjection('EPSG:32610').reproject(
             **{
           'crs': 'EPSG:3857',
           'scale': 30}
             ).float().copyProperties(image,['SCENE_CENTER_TIME'])

#%%RGB

def RGBS2(image):
        #image = image.resample('bicubic')
        
        red = image.select(['B4'],['red'])
        blue = image.select(['B2'],['blue'])
        green = image.select(['B3'],['green'])
        
        new_image = ee.Image.cat([red,green,blue])
        return new_image.set({'satellite':'S2', 'date':image.date()}).setDefaultProjection('EPSG:3857').reproject(
      **{
    'crs':'EPSG:3857',
    'scale': 10}
      ).float().resample('bicubic').copyProperties(image,['SCENE_CENTER_TIME'])

def RGBL9(image):
    
        image = image.resample('bicubic')
        
        hsv = image.select(['B4', 'B3', 'B2']).rgbToHsv();
        sharpened = ee.Image.cat([
                hsv.select('hue'), hsv.select('saturation'), image.select('B8')
                ]).hsvToRgb();
        
        
        new_image = ee.Image.cat(sharpened.select(['red','green','blue'],['red','green','blue']))
        return new_image.set({'satellite':'L9', 'date':image.date()}).setDefaultProjection('EPSG:3857').reproject(
              **{
            'crs':'EPSG:3857',
            'scale': 15}
              ).float().copyProperties(image,['SCENE_CENTER_TIME'])

def RGBL7(image):
    
        image = image.resample('bicubic')
        
        hsv = image.select(['B3', 'B2', 'B1']).rgbToHsv();
        sharpened = ee.Image.cat([
                hsv.select('hue'), hsv.select('saturation'), image.select('B8')
                ]).hsvToRgb();
        
        
        new_image = ee.Image.cat(sharpened.select(['red','green','blue'],['red','green','blue']))
        return new_image.set({'satellite':'L7', 'date':image.date()}).setDefaultProjection('EPSG:3857').reproject(
              **{
            'crs':'EPSG:3857',
            'scale': 15}
              ).float().copyProperties(image,['SCENE_CENTER_TIME'])
            
def RGBL5(image):
        red = image.select(['B3'],['red'])
        blue = image.select(['B1'],['blue'])
        green = image.select(['B2'],['green'])
        
        new_image = ee.Image.cat([red,green,blue])
        return new_image.set({'satellite':'L5', 'date':image.date()}).setDefaultProjection('EPSG:3857').reproject(
             **{
           'crs':'EPSG:3857',
           'scale': 30
             }).float().resample('bicubic').copyProperties(image,['SCENE_CENTER_TIME'])

def RGBL8(image):
        image = image.resample('bicubic')
        
        hsv = image.select(['B4', 'B3', 'B2']).rgbToHsv();
        sharpened = ee.Image.cat([
                hsv.select('hue'), hsv.select('saturation'), image.select('B8')
                ]).hsvToRgb();
        
        
        new_image = ee.Image.cat(sharpened.select(['red','green','blue'],['red','green','blue']))
        return new_image.set({'satellite':'L8', 'date':image.date()}).setDefaultProjection('EPSG:3857').reproject(
              **{
            'crs':'EPSG:3857',
            'scale': 15}
              ).float().copyProperties(image,['SCENE_CENTER_TIME'])
