# -*- coding: utf-8 -*-
"""
Created on Sat Nov  2 17:19:42 2019

@author: berna
Computes Fast Fourier Transforms and fits periodic signals.
"""
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from sys import exit
from lab import (chitest, sine, propfit, sampling, srange, grid, errcor,
                 prnpar, prncor, outlier, pltfitres)

''' Variables that control the script '''
DSO = False # Sampling from Digital Oscilloscope
fit = True # attempt to fit the data
log = False # log-scale axis/es
dB = True # plots response y-axis in deciBels
tick = True # manually choose spacing between axis ticks
tex = True # LaTeX typesetting maths and descriptions
re_im = False # Whether to plot (Re, Im) or (Mag, Phs) of FFT
angle = False # Whether to plot phase of FFT or V(t)

# Extrazione dei vettori di dati grezzi
Dir = './RC_int/'
V1, V2 = np.loadtxt(Dir+'syncsin_50.txt', unpack=True)#,
                    #skiprows=2*256, max_rows=250)
V1*=1e-6
x_min = 0; x_max = 1e6
if DSO:
    Dir ='phase_shift/DSO/' 
    V1, V2 = np.genfromtxt(Dir+'DSO1_41.csv', float, delimiter=',',
                     skip_header = 2, usecols=(0, 1), unpack = True)
    x_min = -0.3; x_max = 2.5

# Trasformazione dei dati nelle grandezze da fittare
x = V1
dx = np.full(len(x),(V1[1]-V1[2])/2)
y = V2 - np.mean(V2)
dy = np.full(len(y), 1)
# Estrazione di un sottointervallo di dati
sx, sdx, sy, sdy = srange(x, dx, y, dy, x_min, x_max)

# Grafico preliminare dati
if tex:
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

init=(60, 50, 0., np.mean(sy))
xx = np.linspace(min(x), max(x), 500)
fig, ax = plt.subplots()
ax = grid(ax, xlab='Time [ms]', ylab='ADC voltage [digit]')
ax.errorbar(x, y, dy, dx, 'ko', ms=1.2, elinewidth=0.8, capsize= 1.1,
        ls='',label='data', zorder=5)
ax.plot(x, y, 'gray', ls='-', lw=0.8, alpha = 0.8)
if fit: ax.plot(xx, sine(xx, *init), 'k--', lw=0.8, zorder=10, alpha =0.6,
                label='initial fit')
legend = ax.legend(loc ='best')
if tick:
    ax.xaxis.set_major_locator(plt.MultipleLocator(5e-2))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(1e-2))
    ax.yaxis.set_major_locator(plt.MultipleLocator(20))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(5))
    if DSO:
        ax.xaxis.set_major_locator(plt.MultipleLocator(5e-2))
        ax.xaxis.set_minor_locator(plt.MultipleLocator(1e-2))
        ax.yaxis.set_major_locator(plt.MultipleLocator(0.1))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(2e-2))

# Costanti per la trasformata discreta di Fourier     
fres, frstd = sampling(space=sx, v=True)
# Lunghezza dell'array su cui calcola la trasformata
fftsize = 100000
tran = np.fft.fftshift(np.fft.fft(sy*np.bartlett(len(sy)), fftsize))
freq = np.fft.fftshift(np.fft.fftfreq(fftsize, d=fres))

# plotfft(freq, tran, dB=False, re_im=False)
fft = tran.real if re_im else 20*np.log10(np.abs(tran)) if dB else np.abs(tran)

fig, (ax1, ax2) = plt.subplots(2,1, True, gridspec_kw={'wspace':0.05,
     'hspace':0.05})
ax1 = grid(ax1)
ax1.plot(freq, fft, c='k', lw='0.9')
ax1.set_ylabel('$\widetilde{V}(f)$ Magnitude [%s]' %('dB' if dB else 'arb. un.'))
if re_im: ax1.set_ylabel('Fourier Transform [Re]')    
if log:
    ax1.set_yscale('log')
    if tick:
        ax1.yaxis.set_major_locator(plt.LogLocator(numticks=16))
        ax1.yaxis.set_minor_locator(plt.LogLocator(subs=np.arange(2, 10)*.1,
                                              numticks = 16))
elif tick:
    ax1.yaxis.set_major_locator(plt.MultipleLocator(1e3))
    ax1.yaxis.set_minor_locator(plt.MultipleLocator(2e2))
    if dB:
        ax1.yaxis.set_major_locator(plt.MultipleLocator(20))
        ax1.yaxis.set_minor_locator(plt.MultipleLocator(4))

fft = tran.imag if re_im else np.angle(tran)
ax2.plot(freq, fft, c='k', lw='0.9')
ax2 = grid(ax2, 'Frequency $f$ [Hz]', '$\widetilde{V}(f)$ Phase [rad]')
if re_im: ax2.set_ylabel('Fourier Transform [Im]')    
ax2.set_xlim(0, 1200)
if tick:
    ax2.xaxis.set_major_locator(plt.MultipleLocator(50))
    ax2.xaxis.set_minor_locator(plt.MultipleLocator(10))
    ax2.yaxis.set_major_locator(plt.MultipleLocator(0.5))
    ax2.yaxis.set_minor_locator(plt.MultipleLocator(0.1))

plt.show()
if not fit: exit()
    
#Fit con sinusoide
pars, covm = curve_fit(sine, sx, sy, init, sdy, absolute_sigma = False)
sdy, pars, covm = propfit(sx, sdx, sy, sdy, sine, pars)
perr, pcor = errcor(covm)
print('Parametri del fit:\n', pars)
prnpar(pars, perr, model=sine)
print('Matrice di correlazione:\n', pcor) 
prncor(pcor, model=sine)

res = sy - sine(sx, *pars)
chisq, ndof, resnorm = chitest(sy, sdy, sine(sx, *pars), ddof=len(pars), v=True)
print('Chi quadro ridotto_v:', chisq/ndof)
   
#Plot DV vs t
if tex:
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

fig1, (ax1, ax2) = pltfitres(sx, sdx, sy, sdy, sine, pars)
ax1.set_ylabel('Differenza di potenziale $\Delta V$ [mV]')
legend = ax1.legend(loc ='lower right', framealpha = 0.3)
legend.set_zorder(100)
if log: ax1.set_xscale('log')
elif tick:
    ax1.xaxis.set_major_locator(plt.MultipleLocator(5e-2))
    ax1.xaxis.set_minor_locator(plt.MultipleLocator(1e-2))
    ax1.yaxis.set_major_locator(plt.MultipleLocator(20))
    ax1.yaxis.set_minor_locator(plt.MultipleLocator(5))

ax2.set_xlabel('Tempo $t$ [ms]', x=0.92)
ax2.set_ylabel('Residui')
if tick:
    ax2.xaxis.set_major_locator(plt.MultipleLocator(10))
    ax2.xaxis.set_minor_locator(plt.MultipleLocator(2))
    ax2.yaxis.set_major_locator(plt.MultipleLocator(5))
    ax2.yaxis.set_minor_locator(plt.MultipleLocator(1))

# Fit V(t) con rimozione degli outliers
TT, dTT, VV, dVV, outT, doutT, outV, doutV = outlier(
    sx, sdx, sy, sdy, sine, pars, thr=3, out=True)
pars, covm = curve_fit(sine, TT, VV, pars, dVV, absolute_sigma = False)
perr, pcor = errcor(covm)
print('Parametri del fit:\n', pars)
prnpar(pars, perr, sine)
print('Matrice di correlazione:\n', pcor)
prncor(pcor, sine)

normout = (outV-sine(outT, *pars))/doutV
chisqin, ndof, normin = chitest(VV, dVV, sine(TT, *pars), ddof=len(pars))
print('Chi quadro ridotto:', chisqin/ndof)

# Plot V(t) con outliers
fig2, (ax1, ax2) = pltfitres(TT, dTT, VV, dVV, sine, pars)
ax1.errorbar(outT, outV, doutV, doutT, 'gx',  ms=3, elinewidth=1.,
             capsize=1.5, ls='', label='outliers')
ax1.set_ylabel('Differenza di potenziale $\Delta V$ [mV]')
if log: ax1.set_xscale('log')
elif tick:
    ax1.xaxis.set_major_locator(plt.MultipleLocator(5e-2))
    ax1.xaxis.set_minor_locator(plt.MultipleLocator(1e-2))
    ax1.yaxis.set_major_locator(plt.MultipleLocator(20))
    ax1.yaxis.set_minor_locator(plt.MultipleLocator(5))
legend = ax1.legend(loc ='lower right', framealpha = 0.3)
legend.set_zorder(100)

ax2.errorbar(outT, normout, None, None, 'gx', elinewidth = 0.7, capsize=0.7,
             ms=3., ls='', zorder=5)
ax2.set_xlabel('Tempo $t$ [ms]', x=0.92)
ax2.set_ylabel('Residui')
if tick:
    ax2.xaxis.set_major_locator(plt.MultipleLocator(10))
    ax2.xaxis.set_minor_locator(plt.MultipleLocator(2))
    ax2.yaxis.set_major_locator(plt.MultipleLocator(2))
    ax2.yaxis.set_minor_locator(plt.MultipleLocator(0.5))
    ax2.set_ylim(np.min(normin)-np.std(normin), np.max(normin)+np.std(normin))
plt.show()