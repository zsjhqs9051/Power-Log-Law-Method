# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 00:07:50 2022

@author: chen.li

Program Function:
power law and log law method to calculate 2.5D shear stress from SRH2D data
power law and log law are assumed to be suitable for all flow field, including the near-structure region

"""

import os
import sys
import shutil
import numpy as np
import pandas as pd
from scipy.optimize import minimize,minimize_scalar
from Script.Support import _const,Log
from Script.Support import *
import matplotlib.pyplot as plt

#input the roughness height, unit is in meter
d50 = 1e-3 #meter
ks = 2*d50


#log file
Log = Log()
logPath = 'Log.txt'
Log.loginitial(logPath) 
content = ''

#SRH2D information load
#Data order: Coordinate_1, Coordinate_2,U_ave Depth
#Coordinate_1 and Coordinate_2:
#   Surface Field: X,Y
#   Cross-Section: distance,elevation
i=1
for file in os.listdir('2D-Info'):
    FileName = file
    CaseName = file.split('.')[0]

    Log.Log(logPath,CaseName)
    
    CellX = {}
    CellY = {}
    CellSize = {}
    CellUave = {}
    CellDepth = {}

    
    if os.path.exists('2D-Info/'+FileName):
        Log.Log(logPath,'SRH2D Information Loaded')
        CellX,CellY,CellSize,CellUave,CellDepth,Foption = Info2DLoad(FileName)
        if Foption == 1:
            content = 'Function Code 1: Surface Feild: calculate 2.5D shear for each SRH2D cell\n'
            Log.Log(logPath,content)
        elif Foption == 2:
            content = 'Function Code 2: Cross-Section: calculate 2.5D shear with fixed m index value, m=7\n'
            Log.Log(logPath,content)

    else:
        Log.Log(logPath,'SRH2D Information does NOT exist!!!\nProgram Terminate')
        sys.exit()
    
    #Function Code 3 additional requirement
    if os.path.exists('3D_CFD/U_Profiles'):
        shutil.rmtree('3D_CFD/U_Profiles')

    #2.5D shear stress calculation
    Log.Log(logPath,'Calculating 2.5D Shear Stress')
    CellWss = {}
    for key in CellX.keys():
        Uave = CellUave[key]
        h = CellDepth[key]
        m = 1.0/7.0
        Upower = []
        Ulog = []
        Err = 1
        ymin = 0
        ymax = 0.15*h

        Utau0 = _const.viscosity*1000/ymax
        Iternumber = 0
        if m>0:
            while (Err>1e-6) and (Iternumber <= 100):
                Iternumber = Iternumber + 1
                theta = 5*_const.viscosity/Utau0
                yss = []
                ys = np.linspace(ymin,ymax,50)
                for y in ys:
                    yplus=yPluslog(y,Utau0)
                    if yplus >30 and yplus<1000:
                        yss.append(y)
                args  = {'ys':yss,'Uave':Uave,'h':h,'m':m,'ks':ks}
                if theta > ks:
                    res = minimize(UerrorS,Utau0,args,method='SLSQP')
                    UtauNew = res.x[0]
                    Err = np.abs(UtauNew - Utau0)
                    Utau0 = UtauNew
                    ymin = 30*_const.viscosity/Utau0
                    ymax = 1000*_const.viscosity/Utau0
                    print('smooth')
                else:
                    res = minimize(UerrorR,Utau0,args,method='SLSQP')
                    UtauNew = res.x[0]
                    Err = np.abs(UtauNew - Utau0)
                    Utau0 = UtauNew
                    ymin = ks#30*_const.viscosity/Utau0
                    ymax = 0.15*h#1000*_const.viscosity/Utau0
                    print('rough')
            if Iternumber > 100:
                resS = minimize(UerrorS,Utau0,args,method='SLSQP')
                resR = minimize(UerrorR,Utau0,args,method='SLSQP')
                UtauNew = max(resS.x[0], resR.x[0])
                Err = np.abs(UtauNew - Utau0)
                Utau0 = UtauNew
                ymin = ks#30*_const.viscosity/Utau0
                ymax = 0.15*h#1000*_const.viscosity/Utau0
                print('hybrid')
            print(key,Iternumber,UtauNew)
        else:
            UtauNew = 0
            CellSize[key] = 0
        
        CellWss[key] = _const.density*UtauNew**2
    Log.Log(logPath,'Calculation Finished!')

    Log.Log(logPath,'Cell-weight 2.5D shear stress on the Cross-Section done!')
    if Foption == 2:
        ShearForce = 0
        totalA = 0
        for key in CellSize.keys():
            if not np.isnan(CellWss[key]):
                ShearForce = ShearForce + CellSize[key]*CellWss[key]
                totalA = totalA + CellSize[key]
        CellAverageTau=ShearForce/totalA
    
    Log.Log(logPath,'Calculate details are written to the txt file\n')
    
    #write results to txt file
    FileSumamry = CaseName+'.txt'
    if os.path.exists(FileSumamry):
        os.remove(FileSumamry)
    with open(FileSumamry,'a+') as f:
        if Foption == 2: 
            content1 = '2.5D Tau=\t' + str(CellAverageTau) +'\t(Pa)\n'
            content3 = 'Roughness Height ks = \t'+str(ks)+'\tm\n'
            content2 = 'Cell\tCoordinate_1(m)\tCoordinate_2(m)\tCell Size(m)\tUave(m/s)\tDepth(m)\tm-index\tTau(Pa)\n'
            content = content1 +content3 +content2
        elif Foption == 1:
            content3 = 'Roughness Height ks = \t'+str(ks)+'\tm\n'
            content =content3+ 'Cell\tCoordinate_1(m)\tCoordinate_2(m)\tUave(m/s)\tDepth(m)\tm-index\tTau(Pa)\n'
        f.write (content)
    
    for key in CellX.keys():
        datacell = ''
        if Foption == 2: 
            datacell = key+'\t'+str(CellX[key])+'\t'+str(CellY[key])+'\t'+str(CellSize[key])+'\t'+str(CellUave[key])+'\t'\
            +str(CellDepth[key])+'\t' +str(7)+'\t'+str(CellWss[key])+'\n'
        elif Foption == 1:
            datacell = key+'\t'+str(CellX[key])+'\t'+str(CellY[key])+'\t'+str(CellUave[key])+'\t'\
            +str(CellDepth[key])+'\t' +str(7)+'\t'+str(CellWss[key])+'\n'
        with open(FileSumamry,'a+') as f:
            f.write (datacell)
    #2.5D shear stress distribution
    XX = []
    YY = []
    TT = []
    for key in CellX.keys():
        XX.append(CellX[key])
        YY.append(CellY[key])
        TT.append(CellWss[key])
    T1 = sorted(TT)
    MinInd = int(np.floor(0.1 * len(T1)))
    MaxInd = int(np.ceil(0.9 * len(T1)))
    plt.figure(i)
    if Foption == 1:
        ssize = 0.1
    else:
        ssize = 5
    
    plt.scatter(XX,YY,c=TT,s=ssize,cmap='jet',marker = "o", vmin = T1[MinInd], vmax = T1[MaxInd])
    cbar = plt.colorbar()
    cbar.set_label('2.5D Shear Stress (pa)')
    plt.axis("equal")
    if Foption == 1:
        plt.xlabel('X (m)')
        plt.ylabel('Y (m)')
    else:
        plt.xlabel('Distance (m)')
        plt.ylabel('Elevation (m)')

    plt.title(CaseName+' '+'2.5D Shear Stress Distribution')
    figName = CaseName +'-2.5D Shear Stress Distribution'
    figsave = figName + '.jpg'
    plt.savefig(figsave,dpi = 1024)
    i = i+1
