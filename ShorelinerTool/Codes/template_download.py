#---------------------------------------------------------------------------------------------
#--------------- SCoWI maps generation and extraction from GEE -------------------------------
#--------------- Marcan Graffin (largely inspired by Kilian Vos' CoastSat tool) -------------
#---------------------------------------------------------------------------------------------
import ee
import os
import pickle
import sys
sys.path.insert(0, r"../../../Codes")
import yaml
from yaml.loader import SafeLoader
from functions import Tools, GEE


global inputs
file_inputs = 'config.yaml'
inputs = yaml.load(open(os.path.join('../',file_inputs),'rb'),Loader = SafeLoader)

if not(inputs['Waterline'] or inputs['SandBar'] or inputs['Vegetataion']):
    sys.exit('No Feature to extract, please modify the config.yaml file')


useS2,useL9,useL8,useL7,useL5=Tools.satMissions(inputs['Missions'])

ee.Initialize()

start,end = Tools.dates(inputs['Dates'])



filepath = './'
polygon=pickle.load(open('poly.json','rb'))
transects=pickle.load(open('transects.p','rb'))
if not(inputs['WaterlineIndex'] in os.listdir()) and inputs['Waterline']:
    os.mkdir(inputs['WaterlineIndex'])
    wl_fL5,wl_fL7,wl_fL8,wl_fL9,wl_fS2 = GEE.index_f(inputs['WaterlineIndex'])
if not(inputs['SandBarIndex'] in os.listdir()) and inputs['SandBar']:
    os.mkdir(inputs['SandBarIndex'])
    sb_fL5,sb_fL7,sb_fL8,sb_fL9,sb_fS2 = GEE.index_f(inputs['SandBarIndex'])
if not(inputs['VegetationIndex'] in os.listdir()) and inputs['Vegetation']:
    os.mkdir(inputs['VegetationIndex'])
    vg_fL5,vg_fL7,vg_fL8,vg_fL9,vg_fS2 = GEE.index_f(inputs['VegetationIndex'])

if not('RGB' in os.listdir()) and inputs['RGB']:
    os.mkdir('RGB')
# polygon=pickle.load(open(r'C:\Users\Marca\Desktop\THESE\Shoreliner\data_in\poly_gray.pkl','rb'))
# polygon=[[[[-75.81,36.06],
#             [-75.81,36.15],
#             [-75.69,36.15],
#             [-75.69,36.06]]]]

global polygon_geom
polygon_geom=ee.Geometry.Polygon(polygon)
polygon_tmp = ee.FeatureCollection(ee.Feature(polygon_geom))

if inputs['Waterline']:
    wl_data=ee.ImageCollection([])
if inputs['SandBar']:
    sb_data=ee.ImageCollection([])
if inputs['Vegetation']:
    vg_data=ee.ImageCollection([])
if inputs['RGB']:
    rgb_data=ee.ImageCollection([])


#Image Collection generation and SCoWI computation 
#%%S2 ini

if useS2:
    S2 = ee.ImageCollection("COPERNICUS/S2_HARMONIZED")
    S2 = S2.filterDate(start,end).filter(ee.Filter.bounds(polygon_tmp)).filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE",inputs['MaxCloudCover']))
    if inputs['CloudMasks']:
        S2 = S2.map(GEE.add_cloud_maskS2)
    
    if inputs['Waterline']:
        wl_data=wl_data.merge(S2.map(wl_fS2))
    if inputs['SandBar']:
        sb_data=sb_data.merge(S2.map(sb_fS2))
    if inputs['Vegetation']:
        vg_data.merge(S2.map(vg_fS2))
    if inputs['RGB']:
        rgb_data = rgb_data.merge(S2.map(GEE.RGBS2))
#%%L8 ini

if useL8:
    L8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_TOA")
    L8 = L8.filterDate(start,end).filter(ee.Filter.bounds(polygon_tmp)).filter(ee.Filter.lt('CLOUD_COVER',inputs['MaxCloudCover']))#.map(mask_clouds_landsat)
    if inputs['CloudMasks']:
        L8 = L8.map(GEE.maskLandsat)
    
    if inputs['Waterline']:
        wl_data=wl_data.merge(L8.map(wl_fL8))
    if inputs['SandBar']:
        sb_data=sb_data.merge(L8.map(sb_fL8))
    if inputs['Vegetation']:
        vg_data.merge(L8.map(vg_fL8))
    if inputs['RGB']:
        rgb_data = rgb_data.merge(L8.map(GEE.RGBL8))
#%%L5 ini

if useL5:
    L5 = ee.ImageCollection("LANDSAT/LT05/C02/T1_TOA")
    L5 = L5.filterDate(start,end).filter(ee.Filter.bounds(polygon_tmp)).filter(ee.Filter.lt('CLOUD_COVER',50))#.map(mask_clouds_landsat)
    if inputs['CloudMasks']:
        L5 = L5.map(GEE.maskLandsat)
    
    if inputs['Waterline']:
        wl_data=wl_data.merge(L5.map(wl_fL5))
    if inputs['SandBar']:
        sb_data=sb_data.merge(L5.map(sb_fL5))
    if inputs['Vegetation']:
        vg_data.merge(L5.map(vg_fL5))
    if inputs['RGB']:
        rgb_data = rgb_data.merge(L5.map(GEE.RGBL5))
#%%L7 ini

if useL7:
    L7 = ee.ImageCollection("LANDSAT/LE07/C02/T1_TOA")
    L7 = L7.filterDate(start,end).filter(ee.Filter.bounds(polygon_tmp)).filter(ee.Filter.lt('CLOUD_COVER',50))#.map(mask_clouds_landsat)
    if inputs['CloudMasks']:
        L7 = L7.map(GEE.maskLandsat)
    
    if inputs['Waterline']:
        wl_data=wl_data.merge(L7.map(wl_fL7))
    if inputs['SandBar']:
        sb_data=sb_data.merge(L7.map(sb_fL7))
    if inputs['Vegetation']:
        vg_data.merge(L7.map(vg_fL7))
    if inputs['RGB']:
        rgb_data = rgb_data.merge(L7.map(GEE.RGBL7))
#%%L9 ini

if useL9:
    L9 = ee.ImageCollection("LANDSAT/LC09/C02/T1_TOA")
    L9 = L9.filterDate(start,end).filter(ee.Filter.bounds(polygon_tmp)).filter(ee.Filter.lt('CLOUD_COVER',50))#.map(mask_clouds_landsat)
    if inputs['CloudMasks']:
        L9 = L9.map(GEE.maskLandsat)
    
    if inputs['Waterline']:
        wl_data=wl_data.merge(L9.map(wl_fL9))
    if inputs['SandBar']:
        sb_data=sb_data.merge(L9.map(sb_fL9))
    if inputs['Vegetation']:
        vg_data.merge(L9.map(vg_fL9))
    if inputs['RGB']:
        rgb_data = rgb_data.merge(L9.map(GEE.RGBL9))

#%% IMAGE SAVING ON THE LOCAL COMPUTER/SERVER
if inputs['Waterline']:
    sizeDataset = wl_data.size()
elif inputs['SandBar']:
    sizeDataset = sb_data.size()
elif inputs['Vegetation']:
    sizeDataset = sb_data.size()



if inputs['Waterline']:
    wl_list=wl_data.toList(sizeDataset)
if inputs['SandBar']:
    sb_list = sb_data.toList(sizeDataset)
if inputs['Vegetation']:
    vg_list=vg_data.toList(sizeDataset)
if inputs['RGB']:
    rgb_list=rgb_data.toList(sizeDataset)
length=sizeDataset.getInfo()


filepath_wl = os.path.join(os.getcwd(),inputs['Index'])
filepath_sb = os.path.join(os.getcwd(),inputs['SandBar'])
filepath_vg = os.path.join(os.getcwd(),inputs['Vegetation'])
filepath_rgb = os.path.join(os.getcwd(),'RGB')

for i in range(length):
    
    if inputs['Waterline']:
        img = ee.Image(wl_list.get(i)).getInfo()
        GEE.saveImage(img, i, filepath_wl, inputs['WaterlineIndex'], polygon_geom)
    
    if inputs['SandBar']:
        img = ee.Image(sb_list.get(i)).getInfo()
        GEE.saveImage(img, i, filepath_sb, inputs['SandBarIndex'], polygon_geom)
    
    if inputs['Vegetation']:
        img = ee.Image(vg_list.get(i)).getInfo()
        GEE.saveImage(img, i, filepath_sb, inputs['VegetationIndex'], polygon_geom)
        
    if inputs['RGB']:
        img = ee.Image(rgb_list.get(i)).getInfo()
        GEE.saveImage(img, i, filepath_sb, ['red','green','blue'], polygon_geom)
        
