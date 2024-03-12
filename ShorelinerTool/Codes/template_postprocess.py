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
import matplotlib.pyplot as plt


file_inputs = 'config.yaml'
inputs = yaml.load(open(os.path.join('./',file_inputs),'rb'),Loader = SafeLoader)
tag_idx = inputs['WaterlineIndex']

zref = 0
ploting = True

transects = dict()

for i in os.listdir():
    if i[:5]=='poly_':
        pathtmp = os.path.join(os.getcwd(),i, 'transects.p')
        tmptransect = pickle.load(open(os.path.join(pathtmp),'rb'))
        for j in tmptransect:
            if 'satellite' in tmptransect[j]:
                transects[j] = tmptransect[j]


if inputs['IQR']:
    for i in transects:
        try:
            transects[i]['satellite'] = Tools.IQR(transects[i]['satellite'],inputs)
        except:
            continue

if inputs['modes']:
    for i in transects:
        try:
            transects[i]['satellite'] = Tools.Mode(transects[i]['satellite'],inputs)
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
        for i in transects:
            c+=1
            print(i + '   ' + str(int(c/len(transects))))
            try:
                X = transects[i]['satellite']['SDW_'+tag_idx].copy()
                t = transects[i]['satellite']['dates'].copy()
            except:
                continue
            idxtide = np.logical_and(tidedates>t[0],tidedates<t[-1])
            wl = np.array(Tools.findNearestTimes(tidevalues[idxtide], tidedates[idxtide], t))
            idx_nan = np.arange(len(wl))[np.isnan(wl)]
            wl = np.delete(wl,idx_nan)
            t = np.delete(t,idx_nan)
            X = np.delete(X,idx_nan)
            if inputs['SlopeCalculation']:
                slopes = Tools.rangeSlopes(0.0, 0.2)
                Xall = Tools.wlCorrect(X,wl,slopes)
                freqMax = Tools.wlPeak(t,wl)
                finalSlope,uncSlope = Tools.integratePowerSpectrum(t,Xall,freqMax)
                transects[i]['post'] = dict()
                transects[i]['post']['slope'+tag_idx]=finalSlope
                transects[i]['post']['uncertaintySlope'+tag_idx]=uncSlope
            elif 'situ' in transects[i]:
                Xsitu = transects[i]['situ']['chainage']
                Zsitu = transects[i]['situ']['elevation']
                slopesitu=[]
                for j in range(len(Xsitu)):
                    try:
                        slopesitu.append(Tools.slopeFromProfile(Xsitu[j], Zsitu[j],zref+inputs['MSLOffset']))
                    except:
                        slopesitu.append(slopesitu[-1])
                        continue
                finalSlope = np.nanmean(slopesitu)
            elif inputs['userDefinedSlope']!=0:
                finalSlope = inputs['userDefinedSlope']
            else:
                print('Need slope data for water level correction')
                break
            correction = (wl-zref)/finalSlope
            X += correction
            transects[i]['satellite']['tideCorrected_'+tag_idx] = X
if ploting:
    for i in transects:
        transects[i]['stats'] = Tools.quickCheck(transects[i],inputs,var='tideCorrected_')
        plt.title(i+'_corrected')
pickle.dump(transects,open('transects_post.p','wb'))












