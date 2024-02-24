import pickle
import numpy as np
import os
from osgeo import gdal
import yaml
from yaml.loader import SafeLoader
import datetime
import sys
sys.path.insert(0, r"../../../Codes")
from functions import Tools

file_inputs = 'config.yaml'
inputs = yaml.load(open(os.path.join('../',file_inputs),'rb'),Loader = SafeLoader)

path_index = os.path.join('./',inputs['WaterlineIndex'])

waterline = []
#%%Process
transects = pickle.load(open('transects.p','rb'))
for key in transects:
    transects[key]['satellite'] = {}

sat = []
dates = []
for i in os.listdir(path_index):
        print(i)
        sat.append(i[-6:-4])
        year = int(i[:4])
        month = int(i[4:6])
        day = int(i[6:8])
        
        date = datetime.datetime(year,month,day)
        dates.append(date)
        
        path_tmp_index = os.path.join(path_index,i)
        img = gdal.Open(path_tmp_index)
        SCoWI = img.GetRasterBand(1).ReadAsArray()
        gt = img.GetGeoTransform()
        t_otsu,hist = Tools.otsu(SCoWI)
        
        wl_tmp, wl_noproj_tmp = Tools.getWaterline(SCoWI,t_otsu,gt)
        waterline.append(wl_tmp)
        
for key in transects:
    transects[key]['satellite']['dates'] = np.array(dates)
    transects[key]['satellite']['sat'] = np.array(sat)

transects = Tools.computeIntersection(waterline,transects,sat,inputs)
pickle.dump(transects,open('transects.p','wb'))
