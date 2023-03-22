# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 11:10:21 2023

@author: chen.li
"""

import pandas as pd
import os
import glob
import sys
import shutil
import re

folder_path = 'Data'
data_files = glob.glob(os.path.join(folder_path,'*.dat'))

for datafile in data_files:
    GroupName1 = datafile.split('.')[0]
    GroupName = GroupName1.split('\\')[1]
    if not os.path.exists('input/'+GroupName):
        os.makedirs('input/'+GroupName)
    with open(datafile) as f:
        lines = []
        for line in f.readlines():
            line_replace = line.replace('\t',',').replace('\\',',').replace(' ','')
            line_replace = re.sub('\x00','',line_replace)            #line_replace = line.replace.replace('\\',',')
            if 'Arc' in line_replace:
                linedata=line_replace.split(',')
                lines.append(linedata)
    lines.sort(key=lambda x:(x[1],x[2]))
    for i in range(int(len(lines)/2)):
        j = 2*i
        caseID = lines[j][1]
        name1 = lines[j][2]
        name2 = lines[j+1][2]
        if 'ft' in name1:
            header = 'X_ft,Y_ft,U_ftps,Depth_ft\n'
        else:
            header = 'X_m,Y_m,U_mps,Depth_m\n'
        x=[]
        y=[]
        U=[]
        D=[]
        for k in range(int((len(lines[j])-3)/2)):
            kk = 2*k+3
            if lines[j][kk+1] !='':
                x.append(float(lines[j][kk]))
                y.append(0.0)
                if 'Vel' in name1:
                    U.append(float(lines[j][kk+1]))
                    D.append(float(lines[j+1][kk+1]))
                else:
                    U.append(float(lines[j+1][kk+1]))
                    D.append(float(lines[j][kk+1]))
        caseinputfile = caseID + '.csv'
        if os.path.exists('input/'+GroupName + '/'+caseinputfile):
            os.remove('input/'+GroupName + '/'+caseinputfile)
        with open('input/'+GroupName + '/'+caseinputfile,'w') as f:
            f.write(header)
            for i in range(len(x)):
                content = str(x[i]) + ',' + str(y[i]) + ',' + str(U[i]) + ',' + str(D[i]) + '\n'
                f.write(content)