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

def Info2DLoad(FileName):
    CellX = {}
    CellY = {}
    CellSize = {}
    CellUave = {}
    CellDepth = {}
    Info2DFIle = '2D-Info/'+FileName
    df = pd.read_csv(Info2DFIle)
    
    if df.columns.values[0].find('ft') > 0:
        factor = 0.3048
    else:
        factor = 1
        
    print(factor)
    
    Data2D = df.values.tolist()
    
    Xvalue = [row[0] for row in Data2D]
    setX = set(Xvalue)
    if len(Xvalue) == len(setX):
        Foption = 2
    else:
        Foption = 1
        
    if Foption == 2:
        for i in range(len(Data2D)-1):
            key = 'Cell_'+str(i+1)
            CellX[key] = (Data2D[i][0]+Data2D[i+1][0])*0.5*factor
            CellY[key] = (Data2D[i][1]+Data2D[i+1][1])*0.5*factor
            CellSize[key] = factor*((Data2D[i][0]-Data2D[i+1][0])**2+(Data2D[i][1]-Data2D[i+1][1])**2)**0.5
            CellUave[key] = (Data2D[i][2]+Data2D[i+1][2])*0.5*factor
            CellDepth[key] = (Data2D[i][3]+Data2D[i+1][3])*0.5*factor
    elif Foption ==1:
        for i in range(len(Data2D)-1):
            key = 'Cell_'+str(i+1)
            CellX[key] = Data2D[i][0]*factor
            CellY[key] = Data2D[i][1]*factor
            CellSize[key] = -1
            CellUave[key] = Data2D[i][2]*factor
            CellDepth[key] = Data2D[i][3]*factor

    return CellX,CellY,CellSize,CellUave,CellDepth,Foption


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