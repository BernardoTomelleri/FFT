# -*- coding: utf-8 -*-
"""
Created on Sat Nov  2 17:19:42 2019

@author: berna
Computes Fast Fourier Transforms and fits periodic signals.
"""
from lab import (np, plt, curve_fit, chitest, sine, propfit, grid, errcor,
                 std_unc, prnpar, prncor, outlier, pltfitres, tick, logy,
                 FFT, plotfft, FWHM, optm, cosn, dosc)
from sys import exit

''' Variables that control the script '''
fit = True # attempt to fit the data to function
lock = False # Lock-in amplify data before FFT
tix = False # manually choose spacing between axis ticks
tex = True # LaTeX typesetting maths and descriptions
log = True # log-scale axis/es
dB = False # plots response y-axis in deciBels
function = sine

if lock: from lockin import t as x, Xmag as y, DSO; fit = False
else: from data import sx as x, sdx as dx, sy as y, sdy as dy, DSO

# Grafico preliminare dati
if tex: plt.rc('text', usetex=True); plt.rc('font', family='serif')
if lock: dx = std_unc(x); dy = std_unc(y)
xx = np.linspace(min(x), max(x), 5000)
fig, ax = plt.subplots()
grid(ax, xlab='Time $t$ [s]', ylab='ADC voltage [digit]')
if DSO: ax.set_xlabel('Time $t$ [s]'); ax.set_ylabel('DSO voltage [V]')
ax.errorbar(x, y, dy, dx, 'ko', ms=1.2, elinewidth=0.8, capsize= 1.1,
        ls='',label='data', zorder=5)
ax.plot(x, y, 'gray', ls='-', lw=0.8, alpha = 0.8)
if fit:
    init = (0.7, 168, 0., np.mean(y)) if DSO else (400, 720, 0., np.mean(y), 0.02)
    ax.plot(xx, function(xx, *init), 'k--', lw=0.8, zorder=10, alpha =0.6,
                label='initial fit')
legend = ax.legend(loc ='best')
if tix:
    tick(ax, xmaj=5e-2, ymaj=20, ymin=5)
    if DSO: tick(ax, ymaj=0.1)

# FFT Computation with numpy window functions
freq, tran, fres, frstd = FFT(time=x, signal=y, window=np.bartlett, beta=None)
fmax, fftmax = optm(np.fft.fftshift(freq), np.fft.fftshift(tran), absv=True)
print("Fundamental frequency = %.2f peak magnitude = %.2f" %(fmax, fftmax))
fwhm = FWHM(freq, tran, FT=True); print("FWHM = %.2f Hz" %fwhm)
# FFT and signal plot
fig, (ax1, ax2) = plotfft(freq, tran, signal=(x,dx,y,dy), norm=True,
                          dB=False, re_im=False)
ax2.set_xlim(0, 1000)
if log: logy(ax2, 16)
elif tix:
    tick(ax2, xmaj=100, ymaj=1e3)
    if dB: tick(ax2, ymaj=20, ymin=5)

if tix: 
    tick(ax1, xmaj=1e-2, ymaj=20, ymin=5)
    if DSO: tick(ax1, ymaj=0.2, ymin=5e-2)

if not fit: plt.show(); exit()
#Fit V(t) con sinusoide
pars, covm = curve_fit(function, x, y, init, dy, absolute_sigma = False)
dy, pars, covm = propfit(x, dx, y, dy, function, pars)
perr, pcor = errcor(covm)
print('Parametri del fit:\n', pars)
prnpar(pars, perr, function)
print('Matrice di correlazione:\n', pcor) 
prncor(pcor, function)

res = y - function(x, *pars)
chisq, ndof, resnorm = chitest(y, dy, function(x, *pars), ddof=len(pars), v=True)
print('Chi quadro ridotto: %.1f' %chisq/ndof)
   
#Plot DV vs t
fig1, (ax1, ax2) = pltfitres(x, dx, y, dy, function, pars)
ax1.set_ylabel('Differenza di potenziale $\Delta V$ [digit]')
legend = ax1.legend(loc ='lower right', framealpha = 0.3)
legend.set_zorder(100)
#if log: ax1.set_yscale('log')
if tix: tick(ax1, ymaj=20, ymin=5)

ax2.set_xlabel('Tempo $t$ [ms]', x=0.92); ax2.set_ylabel('Residui')
if tix: tick(ax2, xmaj=10, ymaj=5)

# Fit V(t) con rimozione degli outliers
TT, dTT, VV, dVV, outT, doutT, outV, doutV = outlier(
    x, dx, y, dy, function, pars, thr=3, out=True)
pars, covm = curve_fit(function, TT, VV, pars, dVV, absolute_sigma = False)
perr, pcor = errcor(covm)
print('Parametri del fit:\n', pars)
prnpar(pars, perr, function)
print('Matrice di correlazione:\n', pcor)
prncor(pcor, function)

normout = (outV-function(outT, *pars))/doutV
chisqin, ndof, normin = chitest(VV, dVV, function(TT, *pars), ddof=len(pars))
print('Chi quadro ridotto: %.1f' %chisqin/ndof)

# Plot V(t) con outliers
fig2, (ax1, ax2) = pltfitres(TT, dTT, VV, dVV, function, pars)
ax1.errorbar(outT, outV, doutV, doutT, 'gx',  ms=3, elinewidth=1.,
             capsize=1.5, ls='', label='outliers')
ax1.set_ylabel('Differenza di potenziale $\Delta V$ [digit]')
#if log: ax1.set_yscale('log')
if tix: tick(ax1, ymaj=20, ymin=5)
legend = ax1.legend(loc ='lower right', framealpha = 0.3)
legend.set_zorder(100)

ax2.errorbar(outT, normout, None, None, 'gx', elinewidth = 0.7, capsize=0.7,
             ms=3., ls='', zorder=5)
ax2.set_xlabel('Tempo $t$ [ms]', x=0.92); ax2.set_ylabel('Residui')
if tix:
    tick(ax2, xmaj=10, ymaj=2, ymin=0.5)
    ax2.set_ylim(np.min(normin)-np.std(normin), np.max(normin)+np.std(normin))
plt.show()