import numpy as np
import ee
from skimage.filters import threshold_otsu
from skimage.measure import find_contours
from skimage.transform import AffineTransform
from scipy.ndimage import median_filter
from shapely import geometry
from osgeo import osr
import pandas as pd

def dates(dates):
    start = dates[0]
    end = dates[1]
    eeStart = ee.Date.fromYMD(int(start[:4]),int(start[5:7]),int(start[8:10]))
    eeEnd = ee.Date.fromYMD(int(end[:4]),int(end[5:7]),int(end[8:10]))
    return eeStart,eeEnd

def satMissions(missions):
    S2,L9,L8,L7,L5 = False,False,False,False,False
    if 'S2' in missions:
        S2=True
    if 'L9' in missions:
        L9=True
    if 'L8' in missions:
        L8=True
    if 'L7' in missions:
        L7=True
    if 'L5' in missions:
        L5=True
    return S2,L9,L8,L7,L5

def polyFromTransects(transect, d_ref = 0.05):
    polygon=[]
    classi_transect=[[]]
    lon_transect=[[]]
    lat_transect=[[]]
    first_key=list(transect.keys())[0]
    #print(transect[first_key])
    lon_ref=(transect[first_key]['transect'][0][0]+transect[first_key]['transect'][1][0])/2
    lat_ref=(transect[first_key]['transect'][0][1]+transect[first_key]['transect'][1][1])/2 
    for i in transect:
      lon_tmp=(transect[i]['transect'][0][0]+transect[i]['transect'][1][0])/2
      lat_tmp=(transect[i]['transect'][0][1]+transect[i]['transect'][1][1])/2
      d=np.sqrt((lon_tmp-lon_ref)**2+(lat_tmp-lat_ref)**2)
      if d<d_ref:
        lon_transect[-1].append(transect[i]['transect'][0][0])
        lon_transect[-1].append(transect[i]['transect'][1][0])
        lat_transect[-1].append(transect[i]['transect'][0][1])
        lat_transect[-1].append(transect[i]['transect'][1][1])
        classi_transect[-1].append(i)
      else:
        lon_transect.append([transect[i]['transect'][0][0],transect[i]['transect'][1][0]])
        lat_transect.append([transect[i]['transect'][0][1],transect[i]['transect'][1][1]])
        classi_transect.append([i])
        lon_ref=lon_tmp
        lat_ref=lat_tmp
      
    indx_boxes=[]
    rate_overl=0.01
    for j in range(len(classi_transect)):
       if not(lon_transect[j] == [] or lat_transect[j] == []):
         indx_boxes.append(j)
         lon1=max(lon_transect[j])+rate_overl*(max(lon_transect[j])-min(lon_transect[j]))
         lon2=min(lon_transect[j])-rate_overl*(max(lon_transect[j])-min(lon_transect[j]))
         lat1=max(lat_transect[j])+rate_overl*(max(lat_transect[j])-min(lat_transect[j]))
         lat2=min(lat_transect[j])-rate_overl*(max(lat_transect[j])-min(lat_transect[j]))
         deltalat = d_ref-(lat1-lat2)
         deltalon = d_ref-(lon1-lon2)*np.cos(lat1*2*3.1415/360)
         print(deltalat,deltalon)
         delta = abs(deltalat-deltalon)
         if deltalon>0:
            lon1 += deltalon/2/np.cos(lat1*2*3.1415/360)
            lon2 -= deltalon/2/np.cos(lat1*2*3.1415/360)
         if deltalat>0:
             lat1 += deltalat/2
             lat2 -= deltalat/2
         print('idx : '+str(j)+' | lat : [ '+str(lat1)+' , '+str(lat2)+' ]')
         tmp=[[[lon1,lat1],
               [lon1,lat2],
               [lon2,lat2],
               [lon2,lat1],
               [lon1,lat1]]]
         polygon.append(tmp)
    return polygon, classi_transect

def norm(band):
    band_min, band_max = band.min(), band.max()
    return ((band - band_min)/(band_max - band_min))


def otsu(img,ax=[],val=256,ploting=False):
    img_val = img[np.logical_not(np.isinf(img))]
    hist=np.histogram(img_val,val,density=True)
    pdf=runmedian(hist[0],5)
    bins=hist[1]
    bins = np.array([bins[i]+(bins[i+1]-bins[i])/2 for i in range(len(bins)-1)])
    t_otsu=0
    try:
        t_otsu = threshold_otsu(img_val,val)
    except:
        print('fail_normal')
    if ploting :
        ax.plot(bins,pdf)
        maxy = ax.get_ylim()[1]
        ax.plot([t_otsu,t_otsu],[0,maxy],'k',linewidth=2)
        ax.grid()
        ax.set_xlabel('Index Value')
        ax.set_ylabel('Frequency')
        ax.set_ylim([0,maxy])
    return t_otsu,hist

def getWaterline(img,threshold,georef,i='    ',ax=[],MIN_LENGTH_SL=0,ploting=False):
    contours=find_contours(img,threshold)
    
    contours_out = [] #non_projected contours
    for j in contours:
        for k in j:
            contours_out.append(k)
    
    aff_mat=np.array([[georef[1], georef[2], georef[0]],
                            [georef[4], georef[5], georef[3]],
                            [0, 0, 1]])
    
    tform = AffineTransform(aff_mat)
    if type(contours) is list:
        points_converted = []
        # iterate over the list
        for l, arr in enumerate(contours): 
            tmp = arr[:,[1,0]]
            points_converted.append(tform(tmp))
    
    # if single np.array
    elif type(contours) is np.ndarray:
        tmp = contours[:,[1,0]]
        points_converted = tform(tmp)
    
    
    contours_coord=points_converted
    contours_long = [] #projected contours
    
    for l, wl in enumerate(contours_coord):
        coords = [(wl[k,0], wl[k,1]) for k in range(len(wl))]
        line = geometry.LineString(coords) # shapely LineString structure
        if len(coords)> 30:
            contours_long.append(wl)
            
    
    x_points = np.array([])
    y_points = np.array([])
    for k in range(len(contours_long)):
        x_points = np.append(x_points,contours_long[k][:,0])
        y_points = np.append(y_points,contours_long[k][:,1])
    shoreline = np.transpose(np.array([x_points,y_points]))
    contours_out = np.array(contours_out)
    return shoreline, contours_out

def computeIntersection(shorelines, transects, sat_id, inputs):
    """
    Computes the intersection between the 2D shorelines and the shore-normal.
    transects. Largely inspired by CoastSat compute_intersections function (Vos et al., 2019)
    """
    if 'Pansharpening' in inputs['Preprocessing']:
        along_dist = {'S2':5,'L5':15,'L7':7.5,'L8':7.5,'L9':7.5}
    else:
        along_dist = {'S2':5,'L5':15,'L7':15,'L8':15,'L9':15}
    # loop through shorelines and compute the median intersection    
    intersections = np.zeros((len(shorelines),len(transects)))
    for i in range(len(shorelines)):
        #print('shoreline n°'+str(i+1)+' over '+str(len(shorelines)))
        sl = shorelines[i]
        max_along = along_dist[sat_id[i]]
        for j,key in enumerate(list(transects.keys())): 
            
            # compute rotation matrix
            X0 = transects[key]['transect_proj'][0,0]
            Y0 = transects[key]['transect_proj'][0,1]
            X1 = transects[key]['transect_proj'][1,0]
            Y1 = transects[key]['transect_proj'][1,1]
            temp = np.array(transects[key]['transect_proj'][-1,:]) - np.array(transects[key]['transect_proj'][0,:])
            phi = np.arctan2(temp[1], temp[0])
            Mrot = np.array([[np.cos(phi), np.sin(phi)],[-np.sin(phi), np.cos(phi)]])
    
            # calculate point to line distance between shoreline points and the transect
            p1 = np.array([X0,Y0])
            p2 = np.array([X1,Y1])
            d_line = np.abs(np.cross(p2-p1,sl-p1)/np.linalg.norm(p2-p1))
            # calculate the distance between shoreline points and the origin of the transect
            d_origin = np.array([np.linalg.norm(sl[k,:] - p1) for k in range(len(sl))])
            # find the shoreline points that are close to the transects and to the origin
            # the distance to the origin is hard-coded here to 1 km 
            idx_along_dist = d_line <= max_along
            idx_cross_dist = d_origin <= ((X1-X0)**2+(Y1-Y0)**2)**.5
            # find the shoreline points that are in the direction of the transect (within 90 degrees)
            temp_sl = sl - p1
            scal = np.dot(temp_sl,p2-p1)
            idx_angle = scal>=0
            # combine the transects that are close in distance and close in orientation
            idx_close = np.where(np.logical_and(idx_along_dist,idx_cross_dist,idx_angle))[0]     
            # idx_close = np.where(idx_dist)[0]
            
            # in case there are no shoreline points close to the transect 
            if len(idx_close) == 0:
                intersections[i,j] = np.nan
            else:
                # change of base to shore-normal coordinate system
                xy_close = np.array([sl[idx_close,0],sl[idx_close,1]]) - np.tile(np.array([[X0],
                                   [Y0]]), (1,len(sl[idx_close])))
                xy_rot = np.matmul(Mrot, xy_close)
                # compute the median of the intersections along the transect
                if i==0 :
                    intersections[i,j] = np.nanmax(xy_rot[0,:])
                else:
                    if not(np.isnan(intersections[i-1,j])):
                        idx = (xy_rot[0,:] - intersections[i-1,j]).tolist().index(min(xy_rot[0,:] - intersections[i-1,j]))
                        intersections[i,j] = xy_rot[0,idx]
                    else:
                        intersections[i,j] = np.nanmax(xy_rot[0,:])
    for j,key in enumerate(list(transects.keys())):
        transects[key]['satellite']['SDW_'+inputs['Index']] = intersections[:,j]
    return transects

def convert_epsg(points, epsg_in, epsg_out):
    """
    Converts from one spatial reference to another using the epsg codes
    
    KV WRL 2018

    Arguments:
    -----------
    points: np.array or list of np.ndarray
        array with 2 columns (rows first and columns second)
    epsg_in: int
        epsg code of the spatial reference in which the input is
    epsg_out: int
        epsg code of the spatial reference in which the output will be            
                
    Returns:    
    -----------
    points_converted: np.array or list of np.array 
        converted coordinates from epsg_in to epsg_out
        
    """
    
    # define input and output spatial references
    inSpatialRef = osr.SpatialReference()
    inSpatialRef.ImportFromEPSG(epsg_in)
    outSpatialRef = osr.SpatialReference()
    outSpatialRef.ImportFromEPSG(epsg_out)
    # create a coordinates transform
    coordTransform = osr.CoordinateTransformation(inSpatialRef, outSpatialRef)
    # if list of arrays
    if type(points) is list:
        points_converted = []
        # iterate over the list
        for i, arr in enumerate(points): 
            points_converted.append(np.array(coordTransform.TransformPoints(arr)))
    # if single array
    elif type(points) is np.ndarray:
        points_converted = np.array(coordTransform.TransformPoints(points))  
    else:
        raise Exception('invalid input type')

    return points_converted

def runnanmean(X,n):
    if isinstance(X,pd.core.series.Series):
        return pd.Series(np.convolve(X, np.ones(n)/n, mode='same'),list(X.index))
    else:
        tmp = X
        for i in range(len(X)):
            a = max(0,i-n//2)
            b = min(len(X),i+n//2)
            idx = np.arange(a,b)
            try:
                X[i] = np.nanmean(tmp[idx])
            except:
                X[i] = np.nan
                continue
        return X

def runmedian(X,n):
    if isinstance(X,pd.core.series.Series):
        return pd.Series(median_filter(X[X != np.nan], n),list(X.index))
    else:
        tmp = X
        for i in range(len(X)):
            a = max(0,i-n//2)
            b = min(len(X),i+n//2)
            idx = np.arange(a,b)
            try:
                X[i] = np.nanmedian(tmp[idx])
            except:
                X[i] = np.nan
                continue
    return X