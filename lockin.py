# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 17:46:05 2020
@author: berni
"""
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from lab import chitest, propfit, outlier, errcor, prncor, prnpar, pltfitres, srange

""" Estrazione dei dati corrispondenti all'intervallo di ddp selezionato
tempo in microsecondi, differenza di potenziale in unità arbitrarie """
Dir = './RLC_data/'
t, DV  = np.loadtxt(Dir+'harm0.1uF_b5.txt', unpack=True)
dt = np.ones_like(t); dDV = np.ones_like(DV)
t *= 1e-3; dt *= 1e-3
t_min = 0.; t_max = 45.
st, sdt, sDV, sdDV = srange(t, dt, DV, dDV, t_min, t_max)

# Fit oscillazioni RLC smorzate di ddp nel tempo
def osc(t, A, tau, omega, f, B):
    return A*np.exp(-t/tau)*np.sin(omega*t +f) + B

def damp(t, A, tau, B):
    return A*np.exp(-t/tau) + B

''' Variables that control the script '''
DSO = False # Sampling from Digital Oscilloscope
fit = False # attempt to fit the data
log = False # log-scale axis/es
tix = False # manually choose spacing between axis ticks
tex = True # LaTeX typesetting maths and descriptions

init = (900., 4., 2., 0., 450.)
pars, covm = curve_fit(osc, st, sDV, init, sdDV, absolute_sigma = True)
print('Parametri del fit:\n', pars)
perr, corm = errcor(covm)
prnpar(pars, perr, osc)
res = sDV - osc(st, *pars)
chisq, ndof, resnorm = chitest(sDV, sdDV, osc(st, *pars), len(pars))
print('Chi quadro ridotto:', (chisq/ndof))
# Covarianza tra i parametri
print('Matrice di correlazione:\n', corm)    
prncor(corm, osc)

# Plot DV vs t con esponenziali al contorno
tt = np.linspace(min(st-sdt), max(st+sdt), 2000)
if tex:
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    
fig, (ax1, ax2) = pltfitres(st, sdt, sDV, sdDV, osc, pars)
ax1.set_ylabel('Differenza di potenziale $\Delta V$ [a.u.]')
if fit:   
    ax1.plot(tt, damp(tt, pars[0, 1, 3]), c='b')
    ax1.plot(tt, -damp(tt, *[pars[:2], -pars[3]]), c='b')
if tix:
    ax1.yaxis.set_major_locator(plt.MultipleLocator(100.))
    ax1.yaxis.set_minor_locator(plt.MultipleLocator(20.))
legend = ax1.legend(loc='best')

ax2.set_xlabel('Tempo $t$ [ms]', x=0.92)
ax2.set_ylabel('Residui [a.u.]')
if tix:
    ax2.xaxis.set_major_locator(plt.MultipleLocator(1.))
    ax2.xaxis.set_minor_locator(plt.MultipleLocator(0.2))
    ax2.yaxis.set_major_locator(plt.MultipleLocator(5.))
    ax2.yaxis.set_minor_locator(plt.MultipleLocator(1.))
sdDV, pars, covm = propfit(st, sdt, sDV, sdDV, osc, pars)

# Fit con rimozione outliers
TT, dTT, VV, dVV, outT, doutT, outV, doutV = outlier(
    st, sdt, sDV, sdDV, osc, pars, thr=3.5, out=True)
pars, covm = curve_fit(osc, TT, VV, pars, dVV, absolute_sigma = True)
perr, corm = errcor(covm)
print('Parametri del fit:\n', pars)
prnpar(pars, perr, osc)
normout = (outV-osc(outT, *pars))/doutV
chisqin, ndof, normin = chitest(VV, dVV, osc(TT, *pars), ddof=len(pars))
print('Chi quadro ridotto:', chisqin/ndof)
print('Matrice di correlazione:\n', corm)
prncor(corm, osc)

fig2, (ax1, ax2) = pltfitres(TT, dTT, VV, dVV, osc, pars)
ax1.errorbar(outT, outV, doutV, doutT, 'gx',  ms=3.5, elinewidth=1.,
             capsize=1.5, ls='', label='outliers')
ax1.set_ylabel('Differenza di potenziale $\Delta V$ [a.u.]')
legend = ax1.legend(loc='best')
if tix:
    ax1.yaxis.set_major_locator(plt.MultipleLocator(100.))
    ax1.yaxis.set_minor_locator(plt.MultipleLocator(20.))

ax2.errorbar(outT, normout, None, None, 'gx', elinewidth = 0.7, capsize=0.7,
             ms=3., ls='', zorder=5)
ax2.set_xlabel('Tempo $t$ [ms]', x=0.92)
ax2.set_ylabel('Residui')
if tix:
    ax2.ticklabel_format(axis='both', style='sci', scilimits=None,
                         useMathText=True)
    ax2.xaxis.set_major_locator(plt.MultipleLocator(1))
    ax2.xaxis.set_minor_locator(plt.MultipleLocator(0.2))
    ax2.yaxis.set_major_locator(plt.MultipleLocator(1.))
    ax2.yaxis.set_minor_locator(plt.MultipleLocator(0.2))
    ax2.set_ylim(np.min(normin)-np.std(normin), np.max(normin)+np.std(normin))
plt.show()