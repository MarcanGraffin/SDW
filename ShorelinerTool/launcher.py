import os
import pickle
import shutil
import geopandas as gpd
from shapely.geometry import Polygon
import yaml
from yaml.loader import SafeLoader
import sys
sys.path.insert(0, r".\Codes")
from functions import Tools

force=True

path_codes = './Codes'
path_inputs  = './Inputs'
file_inputs = 'config.yaml'
inputs = yaml.load(open(os.path.join(path_inputs,file_inputs),'rb'),Loader = SafeLoader)

path_project = os.path.join('./Projects',inputs['Project'])

try :
    os.mkdir(path_project)
    print('Project '+inputs['Project']+' has been created in : '+inputs['PathProject'])
except:
    if force:
        pass
    else:
        print('The project : '+inputs['Project']+' already exists, if you wanna continue put force=True in this script.')
        raise
    
shutil.copy(os.path.join(path_inputs,'config.yaml'),os.path.join(path_project,'config.yaml'))
transects = pickle.load(open(os.path.join(inputs['PathTransect'][0],inputs['PathTransect'][1]),'rb'))
ROIs,TRs = Tools.polyFromTransects(transects,d_ref = inputs['SizeROI'])

if inputs['SaveShpROI']:
    polygons = [Polygon(coords[0]) for coords in ROIs]
    gdf = gpd.GeoDataFrame(geometry=polygons,crs = 4326)
    gdf.to_file(os.path.join(path_project,'ROIs.shp'))
    
c=0

for i in ROIs:
    TR_out = dict()
    for j in TRs[c]:
        TR_out[j] = transects[j]
        TR_out[j]['transect_proj'] = Tools.convert_epsg(TR_out[j]['transect'][:,::-1],4326,3857)[:,:2]
    tmp = 'poly_'+str(c)
    path_tmp = os.path.join(path_project,tmp)
    try:
        os.mkdir(path_tmp)
    except:
        continue
    pickle.dump(i,open(os.path.join(path_tmp,'poly.json'),'wb'))
    pickle.dump(TR_out,open(os.path.join(path_tmp,'transects.p'),'wb'))
    shutil.copy(os.path.join(path_codes,'template_download.py'),os.path.join(path_tmp,'download.py'))
    shutil.copy(os.path.join(path_codes,'template_main_process.py'),os.path.join(path_tmp,'main_process.py'))
    shutil.copy(os.path.join(path_codes,'template_ploting.py'),os.path.join(path_tmp,'ploting.py'))
    c+=1