# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 15:04:50 2022

@author: chen.li
"""
import os
import shutil
import numpy as np
import pandas as pd
from scipy import optimize

class Log():
    def __init__(self):
        pass
    
    def loginitial(self,path):
        if os.path.exists(path):
            os.remove(path)
            
    def Log(self,path,content=''):
        print(content+'\n')
        with open(path,'a+') as f:
            f.write(content+'\n')

class _const:
    density = 1000
    viscosity = 1e-6
    kappa = 0.41
    Cplus = 5

def powerfit(x,m):
    return m*x

def Info2DLoad(FileName,Foption):
    CellX = {}
    CellY = {}
    CellSize = {}
    CellUave = {}
    CellDepth = {}
    Info2DFIle = '2D-Info/'+FileName
    df = pd.read_csv(Info2DFIle)
    Data2D = df.values.tolist()
    
    if Foption == 2 or Foption == 3:
        for i in range(len(Data2D)-1):
            key = 'Cell_'+str(i+1)
            CellX[key] = (Data2D[i][0]+Data2D[i+1][0])*0.5*0.3048
            CellY[key] = (Data2D[i][1]+Data2D[i+1][1])*0.5*0.3048
            CellSize[key] = 0.3048*((Data2D[i][0]-Data2D[i+1][0])**2+(Data2D[i][1]-Data2D[i+1][1])**2)**0.5
            CellUave[key] = (Data2D[i][2]+Data2D[i+1][2])*0.5*0.3048
            CellDepth[key] = (Data2D[i][3]+Data2D[i+1][3])*0.5*0.3048
    elif Foption ==1:
        for i in range(len(Data2D)-1):
            key = 'Cell_'+str(i+1)
            CellX[key] = Data2D[i][0]*0.3048
            CellY[key] = Data2D[i][1]*0.3048
            CellSize[key] = -1
            CellUave[key] = Data2D[i][2]*0.3048
            CellDepth[key] = Data2D[i][3]*0.3048

    return CellX,CellY,CellSize,CellUave,CellDepth

def SampleU(CellX,WSE):
    route = '3D_CFD/sampleDictU_template'
    route2 = '3D_CFD/sampleDictU'
    Top = str(WSE * 0.3048) + ')'
    LineTemp = "$LineX$\n{type lineCellFace;axis z;start $start$;end $end$;}"
    i = 1
    LineSample = '\n'
    for key in CellX.keys():
        LineNo = 'Line'+str(i)
        start = '(0 ' + str(CellX[key]) + ' 0)'
        end = '(0 ' + str(CellX[key]) + ' ' + Top
        LineSample = LineSample + LineTemp.replace('$LineX$', LineNo)\
                                        .replace('$start$', start)\
                                        .replace('$end$', end) + '\n'
        i = i+1
    with open(route) as file:
        content = file.read()
        content = content.replace('$SampleLines$',LineSample)
    with open(route2,'w') as f:
        f.write(content)

def UerrorS(Utau,args):
    ys = args['ys']
    Uave = args ['Uave']
    h = args['h']
    m = args ['m']
    ks = args ['ks']
    Uerror = []
    Umax = Uave/((0.4)**m)
    for y in ys:
        yplus = yPluslog(y,Utau)
        Upower = PowerU(m,Uave,h,y)
        Ulog = SmoothLogU(Utau,yplus)
        Uerr = (Upower - Ulog)**2
        Uerror.append(Uerr)
    return np.sum(Uerror)

def UerrorR(Utau,args):
    ys = args['ys']
    Uave = args ['Uave']
    h = args['h']
    m = args ['m']
    ks = args ['ks']
    Uerror = []
    Umax = Uave/((0.4)**m)
    for y in ys:
        yplus = yPluslog(y,Utau)
        Upower = PowerU(m,Uave,h,y)
        Ulog = RoughLogU(y,Utau,ks)
        Uerr = (Upower - Ulog)**2
        Uerror.append(Uerr)
    return np.sum(Uerror)

def mIndex(Foption,CellDepth,WSE):
    UdataRoute = '3D_CFD/U_Profiles/'
    m = -1
    mfit = {}
    if Foption == 3:
        for key in CellDepth.keys():
            XU = []
            YH = []
            Udata = []
            LineNo = 'Line'+key.split('_')[1]+'_U.csv'
            LineFile = UdataRoute + LineNo
            df = pd.read_csv(LineFile)
            U3D_Data = df.values.tolist()
            if len(U3D_Data)>0:
                for i in range(len(U3D_Data)):
                    Umag = ((U3D_Data[i][1])**2+(U3D_Data[i][2])**2+(U3D_Data[i][3])**2)**0.5
                    if Umag > 1e-4:
                        Udata.append(Umag)
                        yh = 1-(WSE*0.3048-U3D_Data[i][0])/CellDepth[key]
                        if yh <1e-6:
                            yh = 1e-6
                        YH.append(np.log(yh))
                Umax = np.max(Udata)
                for U in Udata:
                    XU.append(np.log(U/Umax))
                mm = optimize.curve_fit(powerfit, XU, YH)[0]
                mfit[key] = 1/mm[0]
            else:
                mfit[key] = -1
    else:
        for key in CellDepth.keys():
            mfit[key] = 1/7
    return mfit

def PowerU(m,Uave,h,y):
    Umax = Uave/(0.4)**m
    U = Umax * (y/h)**m
    return U

def yPluslog(y,Utau):
    return y*Utau/_const.viscosity

def SmoothLogU(Utau,yplus):
    Ulog = Utau*(1/_const.kappa * np.log(yplus) + _const.Cplus)
    return Ulog

def RoughLogU(y,Utau,ks):
    ksplus = ks*Utau/_const.viscosity
    Cs = 0.253
    E = 9.8
    yPlus = y*Utau/_const.viscosity
    if ksplus < 2.25:
        B = 1
    elif ksplus >90:
        B = 1 + Cs*ksplus
    else:
        FactorA = np.sin(0.4258*(np.log(ksplus)-0.811))
        FactorB = (ksplus - 2.25)/87.75 + Cs*ksplus
        B = FactorB ** FactorA
    Ulog = Utau*(1.0/_const.kappa * np.log(E*yPlus/B))
    return Ulog