import pickle
import numpy as np
from osgeo import gdal, osr
from shapely import geometry

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

def smallest_rectangle(polygon):
    """
    Converts a polygon to the smallest rectangle polygon with sides parallel
    to coordinate axes.
     
    KV WRL 2020

    Arguments:
    -----------
    polygon: list of coordinates 
        pair of coordinates for 5 vertices, in clockwise order,
        first and last points must match     
                
    Returns:    
    -----------
    polygon: list of coordinates
        smallest rectangle polygon
        
    """
    
    multipoints = geometry.Polygon(polygon[0])
    polygon_geom = multipoints.envelope
    coords_polygon = np.array(polygon_geom.exterior.coords)
    polygon_rect = [[[_[0], _[1]] for _ in coords_polygon]]
    return polygon_rect

with open('C:/Users/Marca/Desktop/THESE/Shoreliner/data_in/transects_USWest.p','rb') as f:
    transect=pickle.load(f)

tmp=list(transect.values())
    #print(tmp)
    
maxX=max(tmp[0]['transect'][0][0],tmp[0]['transect'][1][0])
minX=min(tmp[0]['transect'][0][0],tmp[0]['transect'][1][0])
maxY=max(tmp[0]['transect'][0][1],tmp[0]['transect'][1][1])
minY=min(tmp[0]['transect'][0][1],tmp[0]['transect'][1][1])

print('x  :' + str(minX)+ ' and ' + str(maxX))
print('y  :' + str(minY)+ ' and ' + str(maxY))
transect_3857=dict()
for i in transect:
  transect_3857[i]=transect
for i in transect:
    tmp=transect[i]# -5900 because transect from 14,273 to 20,163 (last transect) induce bugs du to their geographical positions
    if max([tmp['transect'][0][0],tmp['transect'][1][0]])>maxX:
        maxX=max([tmp['transect'][0][0],tmp['transect'][1][0]])
    if max([tmp['transect'][0][1],tmp['transect'][1][1]])>maxY:
        maxY=max([tmp['transect'][0][1],tmp['transect'][1][1]])
    if min([tmp['transect'][0][0],tmp['transect'][1][0]])<minX:
        minX=min([tmp['transect'][0][0],tmp['transect'][1][0]])
    if min([tmp['transect'][0][1],tmp['transect'][1][1]])<minY:
        minY=min([tmp['transect'][0][1],tmp['transect'][1][1]])
        
print('x  :' + str(minX)+ ' and ' + str(maxX))
print('y  :' + str(minY)+ ' and ' + str(maxY))

# borders=convert_epsg(np.array([[maxX,maxY],[minX,minY]]), 3857, 3857)
# print(borders)
# deltaY1=np.arange(borders[1][1]-0.002,borders[0][1]-0.002,0.05)
# deltaY2=np.arange(borders[1][1]+0.052,borders[0][1]+0.052,0.05)
# deltaY=[[deltaY1[i],deltaY2[i]] for i in range(len(deltaY1))]
    
#    #print(deltaY)
#    classi_transect=[]
#    coordo_transect=[]
#    polygon=[]
#
#    list_all=[j for j in range(len(deltaY))]
#    print(list_all)
#    for j in list_all :#[range(len(deltaY))]:
#        classi_transect.append([])
#        coordo_transect.append([])
#        for i in transect:
#        # pt1=[transect[i][0][0],transect[i][1][0]]
#        # pt2=[transect[i][0][1],transect[i][1][1]]
#            pt1=SDS_tools.convert_epsg(np.array([[transect[i][0][0],transect[i][0][1]]]),3857,3857)[0]
#            pt2=SDS_tools.convert_epsg(np.array([[transect[i][1][0],transect[i][1][1]]]),3857,3857)[0]
#            if (pt1[1]>deltaY[j][0] and pt2[1]>deltaY[j][0] and pt1[1]<deltaY[j][1] and pt2[1]<deltaY[j][1]):
#                classi_transect[j].append(i)
#                coordo_transect[j].append(pt1[0])
#                coordo_transect[j].append(pt2[0])
#    indx_boxes=[]
#    print(pt1[0])
#    for j in list_all:
#        if not(coordo_transect[j]== []):
#            indx_boxes.append(j)
#            lon1=max(coordo_transect[j])+0.1*(max(coordo_transect[j])-min(coordo_transect[j]))
#            lon2=min(coordo_transect[j])-0.1*(max(coordo_transect[j])-min(coordo_transect[j]))
#            lat1=deltaY[j][0]
#            lat2=deltaY[j][1]
#            print('idx : '+str(j)+' | lat : [ '+str(lat1)+' , '+str(lat2)+' ]')
#            tmp = [[[lon1,lat1],
#                    [lon1,lat2],
#                    [lon2,lat2],
#                    [lon2,lat1],
#                    [lon1,lat1]]]
#            polygon.append(SDS_tools.smallest_rectangle(tmp))
polygon=[]
classi_transect=[[]]
lon_transect=[[]]
lat_transect=[[]]
lon_ref=(transect['N_1']['transect'][0][0]+transect['N_1']['transect'][1][0])/2
lat_ref=(transect['N_1']['transect'][0][1]+transect['N_1']['transect'][1][1])/2
d_ref=0.08
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
for j in range(len(classi_transect)):
   if not(lon_transect[j] == [] or lat_transect[j] == []):
     indx_boxes.append(j)
     lon1=max(lon_transect[j])+0.1*(max(lon_transect[j])-min(lon_transect[j]))
     lon2=min(lon_transect[j])-0.1*(max(lon_transect[j])-min(lon_transect[j]))
     lat1=max(lat_transect[j])+0.1*(max(lat_transect[j])-min(lat_transect[j]))
     lat2=min(lat_transect[j])-0.1*(max(lat_transect[j])-min(lat_transect[j]))
     print('idx : '+str(j)+' | lat : [ '+str(lat1)+' , '+str(lat2)+' ]')
     tmp=[[[lon1,lat1],
           [lon1,lat2],
           [lon2,lat2],
           [lon2,lat1],
           [lon1,lat1]]]
     polygon.append(smallest_rectangle(tmp))

pickle.dump(polygon,open(r'C:\Users\Marca\Desktop\THESE\Shoreliner\data_in\poly_WestCoast.pkl','wb'))
