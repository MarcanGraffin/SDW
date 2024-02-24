import matplotlib.pyplot as plt
import pickle
import numpy as np
import os
from osgeo import gdal
import xarray as xr
import yaml
from yaml.loader import SafeLoader
from skimage.transform import AffineTransform
import rasterio
import datetime
import pandas as pd
import sys
sys.path.insert(0, r"../../../Codes")
from functions import Tools

file_inputs = 'config.yaml'
inputs = yaml.load(open(os.path.join('../',file_inputs),'rb'),Loader = SafeLoader)

path_index = os.path.join('./',inputs['Index'])

if inputs['RGB']:
    path_rgb = os.path.join('./','RGB')

if inputs['Ploting']:
    plt.ion()
else:
    plt.ioff()

waterline = []
waterline_noproj = []


#%%Process



months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']



transects = pickle.load(open('transects.p','rb'))
for key in transects:
    transects[key]['satellite'] = {}

TR=dict()
if inputs['SpecifyTransects']:
    names = inputs['NameTransects']
    for j in names:
        TR[j] = transects[j]
else:
    key = list(transects.keys())[len(transects)//2]
    TR[key] = transects[key]

main = list(TR.keys())[0]
lonmax = max([TR[main]['transect_proj'][0,0],TR[main]['transect_proj'][1,0]])
lonmin = min([TR[main]['transect_proj'][0,0],TR[main]['transect_proj'][1,0]])
meanlon = (lonmax+lonmin)/2
latmax = max([TR[main]['transect_proj'][0,1],TR[main]['transect_proj'][1,1]])
latmin = min([TR[main]['transect_proj'][0,1],TR[main]['transect_proj'][1,1]])
meanlat = (latmax+latmin)/2
delta = 0.6*max([abs(lonmin-lonmax),abs(latmin-latmax)])




path_plots = os.path.join('./','plots')
if inputs['SaveFig']:
    try:
        os.mkdir(path_plots)
    except:
        pass


if inputs['RGB']:
    schema = """
    ABCC
    ABCC
    abEE
    """
else:
    schema = """
    .BCC
    .BCC
    .bEE
    """

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
        t_otsu,t_opti,hist = Tools.refinedOtsu(SCoWI)
        
        wl_tmp, wl_noproj_tmp = Tools.getWaterline(SCoWI,t_opti,gt)
        waterline.append(wl_tmp)
        waterline_noproj.append(wl_noproj_tmp)
        
for key in transects:
    transects[key]['satellite']['dates'] = np.array(dates)
    transects[key]['satellite']['sat'] = np.array(sat)

waterline = np.array(waterline)
transects = Tools.compute_intersection(waterline,transects,sat,inputs)

c=-1
for i in os.listdir(path_index):
    print(i)
    
    c+=1
    TS=dict()
    for key in TR:
        TS[key] = pd.series(TR[key]['satellite']['SDW_'+inputs['Index']][:c],TR[key]['satellite']['dates'][:c])
    
    year = int(i[:4])
    month = int(i[4:6])
    day = int(i[6:8])
    
    #date = datetime.datetime(year,month,day)
    
    
    path_tmp_index = os.path.join(path_index,i)
    img = gdal.Open(path_tmp_index)
    SCoWI = img.GetRasterBand(1).ReadAsArray()
    gt = img.GetGeoTransform()
    t_otsu,t_opti,hist = Tools.refinedOtsu(SCoWI)
    
    if inputs['SaveFig']:
        
        if sat=='S2':
            minval = -8000
            maxval = 8000
        else:
            minval=-1
            maxval=1
            
            
        fig,ax = plt.subplot_mosaic(schema,figsize=(20,10))
        a = xr.open_rasterio(path_tmp_index)
        ax['B'].set_aspect('equal')
        ax['B'].pcolormesh(a.x.values,a.y.values,a.values.squeeze(),cmap='Greys_r',vmin=minval,vmax=maxval)
        ax['B'].axis('off')
        ax['B'].set_title(inputs['Index'])
        ax['B'].plot([meanlon-delta,meanlon-delta,meanlon+delta,meanlon+delta,meanlon-delta],
                     [meanlat-delta,meanlat+delta,meanlat+delta,meanlat-delta,meanlat-delta],'-k')
        for key in TR:
            ax['B'].plot(TR[key]['transect_proj'][:,0],TR[key]['transect_proj'][:,1],'r')
        ax['B'].plot(waterline[c][:,0],waterline[c][:,1],'.',markersize=0.5,label=i[:-4])
        
        
        
        ax['b'].pcolormesh(a.x.values,a.y.values,a.values.squeeze(),cmap='Greys_r',vmin=minval,vmax=maxval)
        ax['b'].set_aspect('equal')
        ax['b'].set_ylim([meanlat-delta,meanlat+delta])
        ax['b'].set_xlim([meanlon-delta,meanlon+delta])
        ax['b'].axis('off')
        for key in TR:
            ax['b'].plot(TR[key]['transect_proj'][:,0],TR[key]['transect_proj'][:,1],'r')
        ax['b'].plot(waterline[c][:,0],waterline[c][:,1],'.',markersize=1,label=i[:-4])
        
        for key in TR:
            ax['C'].plot(TS[key],'*-',label = key)
        ax['C'].grid()
        ax['C'].set_xlabel('Time')
        ax['C'].set_ylabel('Waterline position (m)')
        ax['C'].legend()
        
        
        ax['E'].plot(hist[1][:-1],hist[0])
        maxy = ax['E'].get_ylim()[1]
        ax['E'].set_xlabel(inputs['Index']+' values')
        ax['E'].set_xlim([minval,maxval])
        ax['E'].plot([t_otsu,t_otsu],[0,maxy],'k',linewidth=2)
        ax['E'].plot([t_opti,t_opti],[0,maxy],'r',linewidth=2)
        ax['E'].set_ylim([0,maxy])
        ax['E'].grid()
        
        if inputs['RGB']:
            path_tmp_rgb = os.path.join(path_rgb,i)
            
            img_rgb = gdal.Open(path_tmp_rgb)
            red = img_rgb.GetRasterBand(1).ReadAsArray()
            red[np.isinf(red)] = 0
            
            
            green = img_rgb.GetRasterBand(2).ReadAsArray()
            green[np.isinf(green)] = 0
            
            blue = img_rgb.GetRasterBand(3).ReadAsArray()
            blue[np.isinf(blue)] = 0
            
            
            r = Tools.norm(red)
            g = Tools.norm(green)
            b = Tools.norm(blue)
            rgb = np.dstack((r,g,b))
            ax['A'].imshow(rgb)
            ax['A'].axis('off')
            ax['A'].set_title(months[month-1]+' - '+str(year))
            ax['A'].plot(waterline_noproj[c][:,1],waterline_noproj[c][:,0],'.',markersize=0.5,label=i[:-4])
            
            ax['a'].imshow(rgb)
            ax['a'].plot(waterline_noproj[c][:,1],waterline_noproj[c][:,0],'.',markersize=1,label=i[:-4])
            tmp_lon = abs(a.x-meanlon).values.tolist().index(min(abs(a.x-meanlon)))
            tmp_lat = abs(a.y-meanlat).values.tolist().index(min(abs(a.y-meanlat)))
            ax['a'].set_ylim([tmp_lat+delta//gt[1],tmp_lat-delta//gt[1]])
            ax['a'].set_xlim([tmp_lon-delta//gt[1],tmp_lon+delta//gt[1]])
            ax['a'].axis('off')
            
            plt.savefig(os.path.join(path_plots,i[:-4]+'.png'))