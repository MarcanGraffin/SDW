import numpy as np
import ee
from skimage.filters import threshold_otsu
from skimage.measure import find_contours
from skimage.transform import AffineTransform
from scipy.ndimage import median_filter
from shapely import geometry
from osgeo import osr
import pandas as pd
from scipy.stats import linregress
from astropy import timeseries
from scipy import signal, integrate, interpolate
from scipy.signal import find_peaks
import matplotlib.pyplot as plt

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

def refinedOtsu(img,ax=[],val=256,ploting=False):
    img_val = img[np.logical_not(np.isnan(img))]
    hist=np.histogram(img_val,val,density=True)
    pdf=runmedian(hist[0],10)
    bins=hist[1]
    bins = np.array([bins[i]+(bins[i+1]-bins[i])/2 for i in range(len(bins)-1)])
    t_otsu=0
    try:
        t_otsu = threshold_otsu(img_val,val)
    except:
        print('fail_normal')
    t_otsu_pdf = pdf[bins>=t_otsu][0]
    peaks, prop = find_peaks(pdf, height=(t_otsu_pdf, None))
    maxs_histo = bins[peaks]

    res_max_left = maxs_histo[maxs_histo < t_otsu]
    res_max_right = maxs_histo[maxs_histo >= t_otsu]

    if res_max_left.size==0:
        maxl = t_otsu 
    else:
        maxl = res_max_left[-1]

    if res_max_right.size==0:
        maxr = t_otsu
    else:
        maxr = res_max_right[0]

    # select the part of the histogram between 2 peaks around the threshold
    pdf_selected = pdf[(bins>maxl)&(bins<maxr)]
    bins_selected = bins[(bins>maxl)&(bins<maxr)]

    if bins_selected.size==0:
        t_opti = t_otsu
    else:
        min_pdf = np.argmin(pdf_selected)
        t_opti = bins_selected[min_pdf]
    if ploting :
        ax.plot(bins,pdf)
        maxy = ax.get_ylim()[1]
        ax.plot([t_otsu,t_otsu],[0,maxy],'k',linewidth=2)
        ax.plot([t_opti,t_opti],[0,maxy],'r',linewidth=2)
        ax.grid()
        ax.set_xlabel('Index Value')
        ax.set_ylabel('Frequency')
        ax.set_ylim([0,maxy])
    return t_otsu,t_opti,hist

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

def getWaterline(img,threshold,georef,transects,date,inputs,i='    ',ax=[],MIN_LENGTH_SL=0,ploting=False):
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
            if len(tmp)>10:
                points_converted.append(tform(tmp))
    # elif type(contours) is np.ndarray:
    #     tmp = contours[:,[1,0]]
    #     points_converted = tform(tmp)
        
        
        
        
    if ploting:
        plt.figure()
        L,l =  img.shape[0],img.shape[1]
        x = np.arange(georef[0],georef[0]+(L+1)*georef[1],georef[1])
        y = np.arange(georef[3],georef[3]+(l+1)*georef[5],georef[5])
        X,Y = np.meshgrid(x,y)
        quantiles=np.nanquantile(img.flatten(),[0.05,0.95])
        plt.pcolormesh(X,Y,img,cmap='Greys_r',vmin=quantiles[0],vmax=quantiles[1])
        for j in transects:
            lon = transects[j]['transect_proj'][:,0]
            lat = transects[j]['transect_proj'][:,1]
            plt.plot(lon,lat,'r')
            if 'situ' in transects[j]:
                d = ((lon[0]-lon[1])**2+(lat[0]-lat[1])**2)**.5
                tsitu = transects[j]['situ']['dates']
                timg = np.array([date])
                idx = findNearestTimes(range(len(tsitu)), tsitu, timg)[0]
                if not(np.isnan(idx)):
                    X = transects[j]['situ']['chainage'][idx]
                    Z = np.array(transects[j]['situ']['elevation'][idx])
                    Xsitu = getCrossPos(X,Z,inputs['MSLOffset'])
                    posx = lon[0] + Xsitu/d*(lon[1] - lon[0]) 
                    posy = lat[0] + Xsitu/d*(lat[1] - lat[0])
                    plt.plot([posx],[posy],'x',color='orange',markersize=10)
        for j in points_converted:
            plt.plot(j[:,0],j[:,1],'.b')
        plt.title(i)
        
    # if single np.array
    
    
    
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
    transects. Adapted from CoastSat compute_intersections function (Vos et al., 2019)
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
                intersections[i,j] = np.nanmax(xy_rot[0,:])
                        
    for j,key in enumerate(list(transects.keys())):
        transects[key]['satellite']['SDW_'+inputs['WaterlineIndex']] = intersections[:,j]
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

#%% QUICK CHECKS

def quickCheck(profile,inputs,wl=0.,var = 'SDW_', index='SCoWI', ploting = True):
    
    """
    use it as quickCheck(transects[PROFILEYOUWANNASEE] 
                         for visual comparison situ/sat
    """
    z = wl + inputs['MSLOffset']
    Xsat = profile['satellite'][var+index]
    tsat = profile['satellite']['dates']
    tsitu = []
    Z = profile['situ']['elevation']
    X = profile['situ']['chainage']
    Xsitu=[]
    idx_error=[]
    for i in range(len(X)):
        try:
            Xsitu.append(getCrossPos(X[i],Z[i],z))
            tsitu.append(profile['situ']['dates'][i])
        except:
            idx_error.append(i)
            continue
    for i in profile['situ']:
        profile['situ'][i] = np.delete(np.array(profile['situ'][i],dtype='object'),idx_error)
    profile['situ']['situ_'+var+index] = Xsitu
    profile['situ'] = IQR(profile['situ'],inputs,'situ_'+var+index,ratio=1)
    nXsat = findNearestTimes(Xsat, tsat, tsitu)
    idx_nan = np.arange(len(nXsat))[np.isnan(nXsat)]
    nXsat = np.delete(nXsat,idx_nan)
    Xsitu = np.delete(Xsitu,idx_nan)
    tsitu = np.delete(tsitu,idx_nan)
    out = getStats(nXsat,Xsitu)
    maxitmp = max([max(nXsat),max(Xsitu)])
    minitmp = min([min(nXsat),min(Xsitu)])
    mini = minitmp - 0.05*(maxitmp-minitmp)
    maxi = maxitmp + 0.05*(maxitmp-minitmp)
    if ploting:
        fig,ax = plt.subplot_mosaic("AB",figsize=((13,6)))
        ax['A'].plot(tsat,Xsat,'.k',label='satellite')
        ax['A'].plot(tsitu,Xsitu,'.r',label='situ')
        ax['A'].set_xlabel('Dates')
        ax['A'].set_ylabel('Cross-shore position (m)')
        ax['A'].legend()
        ax['B'].plot(nXsat,Xsitu,'.k')
        ax['B'].set_xlabel('Satellite (m)')
        ax['B'].set_ylabel('Situ (m)')
        ax['B'].set_xlim([mini,maxi])
        ax['B'].set_ylim([mini,maxi])
        
        tRMSE = 'RMSE = '+str(round(out['RMSE'],2))+' m'
        tR2 = 'R² = '+str(round(out['R2'],2))
        tbias = 'bias = '+str(round(out['bias'],2))+' m'
        tstd = r'$\sigma$ = '+str(round(out['STD'],2))+' m'
        allt = tR2+'\n'+tRMSE+'\n'+tbias+'\n'+tstd
        plt.text(minitmp,maxitmp - 0.15*(maxitmp-minitmp),allt)
    return(out)


#%% POST PROCESS
def removeNaN(TR,inputs,refmin = 0):
    X = TR['SDW_'+inputs['WaterlineIndex']].copy()
    nans = np.isnan(X)
    belowzero = X<refmin #means can't have a waterline more landward than the transect origin
    cond = np.logical_or(nans,belowzero)
    idx_nan = np.arange(len(X))[cond]
    for i in TR:
        TR[i] = np.delete(TR[i],idx_nan)
    return(TR)


def Mode(TR,inputs,n=0.65,valid=True):
    
    X = TR['SDW_'+inputs['WaterlineIndex']].copy()
    if valid:
        tmpX = X - np.nanmax([0,np.nanmedian(X)])
    else:
        tmpX = X - np.nanmedian(X)
    
    threshold = abs(threshold_otsu(tmpX))
    # print(threshold)
    # print(n*np.std(tmpX))
    idx_otsu = []
    if threshold > n*np.std(tmpX):
        for i in range(len(tmpX)):
            if tmpX[i]>threshold or tmpX[i]<-threshold: #we remove data both sides, but most certainly outliers will be in the tmpX[i]<-threshold  part.
                idx_otsu.append(i)
    for i in TR:
        TR[i] = np.delete(TR[i],idx_otsu)
    return(TR)



def IQR(TR,inputs, varname,val1 = 0.25, val2 = 0.75, ratio=1.5):
    
    """
    data[PROFILE]['satellite'] = IQR(data[PROFILE]['satellite'])
    Clean a SDW timeseries by removing data deviating too much
    from the distribution, clean also the corresponding dates and others...
    """
    X = TR[varname].copy()
    Q1 = np.quantile(X,val1)
    Q3 = np.quantile(X,val2)
    IQR = Q3-Q1
    valmax = Q3 + ratio*IQR
    valmin = Q1 - ratio*IQR
    
    idx_iqr=[]
    for i in range(len(X)):
        if X[i]>valmax or X[i]<valmin:
            idx_iqr.append(i)
            
    for i in TR:
        try:
            TR[i] = np.delete(TR[i],idx_iqr)
        except:
            TR[i] = np.delete(np.array(TR[i],dtype='object'),idx_iqr)
            continue
    
    return TR

def findNearestTimes(data_in,dates_in,dates_out):
    
    """
    new_data = findNearestTimes(old_data,old_dates,new_dates)
    Returns the data
    """
    data_out=[]
    tmp_in=np.array([dates_in[i].timestamp() for i in range(len(dates_in))])
    tmp_out=np.array([dates_out[i].timestamp() for i in range(len(dates_out))])
    
    for i in tmp_out:
        tmpdiff = abs(tmp_in-i)
        if min(tmpdiff)/3600/24<30:
            indx= np.argmin(tmpdiff)
            data_out.append(data_in[indx])
        else:
            data_out.append(np.nan)
    return data_out

def getStats(X,Y):
    
    """
    STATS = getStats(X,Y)
    returns the common validations stats when comparing 2 lists of the same length
    """
    X = np.array(X)
    Y = np.array(Y)
    
    out = dict()
    R = out['R'] = np.corrcoef(X,Y)[0][1]
    R2 = out['R2'] = R**2
    MSE = np.mean((X - Y) ** 2)
    RMSE = out['RMSE'] = np.sqrt(MSE)
    bias = out['bias'] = np.mean(X - Y)
    std = out['STD'] = (MSE - bias**2)**.5
    return out

def slopeFromProfile(X,Z,zref=0,window=0.5):
    
    """
    slope = slopeFromProfile(X,Z)
    returns the slope around the elevation Z=zref
    """
    
    idx = np.logical_and(Z>zref-window,Z<zref+window)
    x=X[idx]
    z=Z[idx]
    slope = linregress(x,z)[0]
    return -slope

def getCrossPos(X,Z,z):
    
    """
    x = getCrossPos(X,Z,z)
    returns the cross-shore (x) position of the intersection between
    the elevation z and the profile(X,Z)
    """
    
    X = np.array(X)
    Z = np.array(Z)
    idx = np.argmin(abs(Z-z))
    z1 = abs(Z[idx-1]-z)
    z2 = abs(Z[idx+1]-z)
    dX = X[idx+1]-X[idx-1]
    return(X[idx-1]+z1/(z1+z2)*dX)

#%% Water level correctionsù

def rangeSlopes(minSlope,maxSlope,deltaSlope=0.0025):
    return np.arange(max(minSlope,deltaSlope),maxSlope+deltaSlope,deltaSlope)

def wlCorrect(Xi,wl,slopes,zref=0.0):
    'apply waterlevel correction with a range of slopes'
    Xall = []
    for i in range(len(slopes)):
        # apply tidal correction
        tmpX=Xi.copy()
        tide_correction = (zref-wl)/slopes[i]
        tmpX += tide_correction
        Xall.append(tmpX)
    return Xall

def wlPeak(dates,wl):
    'find the high frequency peak in the tidal time-series'
    # create frequency grid
    t = np.array([_.timestamp() for _ in dates]).astype('float64')
    #days_in_year = 365.2425
    seconds_in_day = 24*3600
    time_step = 8*seconds_in_day
    freqs = getFrequencyGrid(t,time_step,50)
    # compute power spectrum
    ps_tide,_,_ = powerSpectrum(t,wl,freqs,[])
    # find peaks in spectrum
    idx_peaks,_ = signal.find_peaks(ps_tide, height=0)
    y_peaks = _['peak_heights']
    idx_peaks = idx_peaks[np.flipud(np.argsort(y_peaks))]
    # find the strongest peak at the high frequency (defined by freqs_cutoff[1])
    idx_max = idx_peaks[freqs[idx_peaks] > 1./(seconds_in_day*30)][0]
    # compute the frequencies around the max peak with some buffer (defined by buffer_coeff)
    freqs_max = [freqs[idx_max] - 1e-8, freqs[idx_max] + 1e-8]
    # make a plot of the spectrum
    return freqs_max

def getFrequencyGrid(time,time_step,n0):
    'define frequency grid for Lomb-Scargle transform'
    T = np.max(time) - np.min(time)
    fmin = 1/T
    fmax = 1/(2*time_step) # Niquist criterium
    df = 1/(n0*T)
    N = np.ceil((fmax - fmin)/df).astype(int)
    freqs = fmin + df * np.arange(N)
    return freqs

def powerSpectrum(t,y,freqs,idx_cut):
    'compute power spectrum and integrate'
    model = timeseries.LombScargle(t, y, dy=None, fit_mean=True, center_data=True, nterms=1, normalization='psd')
    ps = model.power(freqs)
    # integrate the entire power spectrum
    E = integrate.simpson(ps, x=freqs, even='avg')
    if len(idx_cut) == 0:
        idx_cut = np.ones(freqs.size).astype(bool)
    # integrate only frequencies above cut-off
    Ec = integrate.simpson(ps[idx_cut], x=freqs[idx_cut], even='avg')
    return ps, E, Ec

def integratePowerSpectrum(dates,Xall,freqMax):
    'integrate power spectrum at the frequency band of peak tidal signal'
    t = np.array([_.timestamp() for _ in dates]).astype('float64')
    seconds_in_day = 24*3600
    time_step = 8*seconds_in_day
    freqs = getFrequencyGrid(t,time_step,8)    
    beach_slopes = rangeSlopes(0.0025, 0.2, 0.0025)
    # integrate power spectrum
    idx_interval = np.logical_and(freqs >= freqMax[0], freqs <= freqMax[1]) 
    E = np.zeros(beach_slopes.size)
    for i in range(len(Xall)):
        ps, _, _ = powerSpectrum(t,Xall[i],freqs,[])
        E[i] = integrate.simpson(ps[idx_interval], x=freqs[idx_interval], even='avg')
    # calculate confidence interval
    delta = 0.0001
    prc = 0.05
    f = interpolate.interp1d(beach_slopes, E, kind='linear')
    beach_slopes_interp = rangeSlopes(0.0025,0.2-delta,delta)
    E_interp = f(beach_slopes_interp)
    # find values below minimum + 5%
    slopes_min = beach_slopes_interp[np.where(E_interp <= np.min(E)*(1+prc))[0]]
    if len(slopes_min) > 1:
        ci = [slopes_min[0],slopes_min[-1]]
    else:
        ci = [beach_slopes[np.argmin(E)],beach_slopes[np.argmin(E)]]
    
    return beach_slopes[np.argmin(E)], ci

