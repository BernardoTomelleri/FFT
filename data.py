# -*- coding: utf-8 -*-
"""
Created on Sat May  9 17:57:56 2020

@author: berni
Reads in input from data files and defines its source for subsequent analysis
"""
from lab import np, srange, std_unc
DSO = False # Sampling from Digital Oscilloscope
AC = False # Cut constant DC offset from zero (AC coupling)

m1 = 1.040; m5 = 4.70
# Extrazione dei vettori di dati grezzi
Dir = './RLC_data/'
V1, V2 = np.loadtxt(Dir +'long0.47uF_b3_new' +'.txt', unpack=True)#,
                    #skiprows=2*256, max_rows=250)
V1*=1e-6
x_min = 20e-6; x_max = x_min + 0.047

if DSO:
    Dir ='phase_shift/DSO/' 
    V1, V2 = np.genfromtxt(Dir +'1_27' +'.csv', float, delimiter=',',
                     skip_header = 2, usecols=(0, 1), unpack = True)
    x_min = -4.3; x_max = -0.2
    
# Trasformazione dei dati nelle grandezze da analizzare
x = V1
dx = std_unc(x) if DSO else np.full(len(x), 4e-6)
y = V2 - np.mean(V2) if AC else V2
dy = np.full(len(y), (V2[1] - V2[0])/4) if DSO else np.ones_like(V2)
if not DSO: y*=m1; dy*=m1
# Estrazione di un sottointervallo di dati
sx, sdx, sy, sdy = srange(x, dx, y, dy, x_min, x_max)