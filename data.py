# -*- coding: utf-8 -*-
"""
Created on Sat May  9 17:57:56 2020

@author: berni
Reads in input from data files and defines its source for subsequent analysis
"""
from lab import np, srange
DSO = True # Sampling from Digital Oscilloscope
AC = True # Cut constant DC offset from zero (AC coupling)

# Extrazione dei vettori di dati grezzi
Dir = './RLC_data/'
V1, V2 = np.loadtxt(Dir +'long0.1uF_b3' +'.txt', unpack=True)#,
                    #skiprows=2*256, max_rows=250)
V1*=1e-6
x_min = 0; x_max = x_min + 0.047

if DSO:
    Dir ='phase_shift/DSO/' 
    V1, V2 = np.genfromtxt(Dir +'DSO1_41' +'.csv', float, delimiter=',',
                     skip_header = 2, usecols=(0, 1), unpack = True)
    x_min = -0.2; x_max = 2.5
    
# Trasformazione dei dati nelle grandezze da analizzare
x = V1
dx = np.full(len(x), (V1[1]-V1[2])/2)
y = V2 - np.mean(V2) if AC else V2
dy = np.full(len(y), (V2[1]-V2[2])/5) if DSO else np.ones_like(V2)
# Estrazione di un sottointervallo di dati
sx, sdx, sy, sdy = srange(x, dx, y, dy, x_min, x_max)