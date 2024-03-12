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
import matplotlib.pyplot as plt

file_inputs = 'config.yaml'
inputs = yaml.load(open(os.path.join('../',file_inputs),'rb'),Loader = SafeLoader)

path_index = os.path.join('./',inputs['WaterlineIndex'])
ploting=True
waterline = []
#%%Process
transects = pickle.load(open('transects.p','rb'))
for key in transects:
    transects[key]['satellite'] = {}

sat = []
dates = []
list_img = os.listdir(path_index)
list_img.sort()
c=-1
for i in list_img:
        c+=1
        if c%10==0 and ploting:
            plting=True
        else:
            plting=False
        print(i)
        
        year = int(i[:4])
        month = int(i[4:6])
        day = int(i[6:8])
        hour = int(i[9:11])
        minute = int(i[11:13])
        
        date = datetime.datetime(year,month,day,hour,minute)
        
        
        path_tmp_index = os.path.join(path_index,i)
        img = gdal.Open(path_tmp_index)
        try:
            SCoWI = img.GetRasterBand(1).ReadAsArray()
        except:
            continue
        idx_nan=0
        for l in range(len(SCoWI)):
            for h in range(len(SCoWI[l])):
                if SCoWI[l,h]==0.0 or np.isnan(SCoWI[l,h]):
                    SCoWI[l,h] = np.nan
                    idx_nan+=1
        rate = idx_nan/len(SCoWI.flatten())
        if rate<0.5:
            gt = img.GetGeoTransform()
            t_otsunaze,t_otsu,hist = Tools.refinedOtsu(SCoWI)
            
            wl_tmp, wl_noproj_tmp = Tools.getWaterline(SCoWI,t_otsu,gt,transects,date,inputs,i=i,ploting=plting)
            waterline.append(wl_tmp)
            dates.append(date)
            sat.append(i[-6:-4])
        
for key in transects:
    transects[key]['satellite']['dates'] = np.array(dates)
    transects[key]['satellite']['sat'] = np.array(sat)

transects = Tools.computeIntersection(waterline,transects,sat,inputs)

for i in transects:
    transects[i]['satellite'] = Tools.removeNaN(transects[i]['satellite'],inputs)
if ploting:
    for i in transects:
        Tools.quickCheck(transects[i],inputs)
        plt.title(i+'_noCorrections')
pickle.dump(transects,open('transects.p','wb'))
