import numpy as np
import pickle
import datetime
import sys
sys.path.insert(0, r"../../Codes")
from functions import Tools
import os
import yaml
from yaml.loader import SafeLoader
import pandas as pd

file_inputs = 'config.yaml'
inputs = yaml.load(open(os.path.join('./',file_inputs),'rb'),Loader = SafeLoader)
tag_idx = inputs['WaterlineIndex']

Zref = 0


transect = dict()

for i in os.listdir():
    if i[:5]=='poly_':
        pathtmp = os.path.join(os.getcwd(),i, 'transects.p')
        tmptransect = pickle.load(open(os.path.join(pathtmp),'rb'))
        for j in tmptransect:
            if 'satellite' in tmptransect[j]:
                transect[j] = tmptransect[j]


if inputs['IQR']:
    for i in transect:
        try:
            tmpX = transect[i]['satellite']['SDW_'+tag_idx].copy()
            IQRidx = Tools.IQR(tmpX)
            for j in transect[i]['satellite']:
                transect[i]['satellite'][j] = np.delete(transect[i]['satellite'][j],IQRidx)
        except:
            continue

if inputs['TideCorrection']:
    if inputs['TidePrediction']:
        print("need to implement FES2022")
    else:
        tidedata = pd.read_csv(inputs['PathTideData'])
        tidedates = np.array([datetime.datetime.strptime(tidedata.dates.values[i],
                        '%Y-%m-%d %H:%M:%S') for i in range(len(tidedata))])
        tidevalues = tidedata.tides.values
        c=-1
        for i in transect:
            c+=1
            print(i + '   ' + str(int(c/len(transect))))
            try:
                X = transect[i]['satellite']['SDW_'+tag_idx].copy()
                t = transect[i]['satellite']['dates'].copy()
            except:
                continue
            idxtide = np.logical_and(tidedates>t[0],tidedates<t[-1])
            wl = np.array(Tools.findNearestTimes(tidevalues[idxtide], tidedates[idxtide], t))
            if inputs['SlopeCalculation']:
                slopes = Tools.rangeSlopes(0.0, 0.2)
                Xall = Tools.wlCorrect(X,wl,slopes)
                freqMax = Tools.wlPeak(t,wl)
                finalSlope,uncSlope = Tools.integratePowerSpectrum(t,Xall,freqMax)
                transect[i]['post'] = dict()
                transect[i]['post']['slope'+tag_idx]=finalSlope
                transect[i]['post']['uncertaintySlope'+tag_idx]=uncSlope
            elif 'situ' in transect[i]:
                Xsitu = transect[i]['situ']['chainage']
                Zsitu = transect[i]['situ']['elevation']
                slopesitu=[]
                for j in range(len(Xsitu)):
                    try:
                        slopesitu.append(Tools.slopeFromProfile(Xsitu, Zsitu,zref=Zref+inputs['MSLOffset']))
                    except:
                        slopesitu.append(slopesitu[-1])
                        continue
                finalSlope = np.nanmean(slopesitu)
            elif inputs['userDefinedSlope']!=0:
                finalSlope = inputs['userDefinedSlope']
            else:
                print('Need slope data for water level correction')
                break
            correction = (zref-wl)/finalSlope
            X += correction
            transect[i]['satellite']['SDW_'+tag_idx+'_tCorr'] = X

pickle.dump(transect,open('transects_post.p','wb'))












