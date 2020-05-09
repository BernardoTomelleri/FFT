# -*- coding: utf-8 -*-
"""
Created on Sat Nov  2 17:19:42 2019

@author: berna
Computes Fast Fourier Transforms and fits periodic signals.
"""
from lab import (np, plt, curve_fit, chitest, sine, propfit, grid, errcor,
                 prnpar, prncor, outlier, pltfitres, tick, logy, dosc)
from sys import exit

''' Variables that control the script '''
fit = False # attempt to fit the data
log = True # log-scale axis/es
dB = False # plots response y-axis in deciBels
tix = False # manually choose spacing between axis ticks
tex = True # LaTeX typesetting maths and descriptions
re_im = False # Whether to plot (Re, Im) or (Mag, Phs) of FFT
angle = False # Whether to plot phase of FFT or V(t)
lock = True # Lock-in amplify data before FFT

if lock: from lockin import t as x, Xmag as y, fres, fftsize, DSO
else: from data import x, dx, y, dy, sx, sdx, sy, sdy, init, fres, fftsize, DSO

# Grafico preliminare dati
if tex: plt.rc('text', usetex=True); plt.rc('font', family='serif')
if lock: dx = dy = None
xx = np.linspace(min(x), max(x), 500)
fig, ax = plt.subplots()
grid(ax, xlab='Time $t$ [s]', ylab='ADC voltage [digit]')
if DSO: ax.set_xlabel('Time $t$ [s]'); ax.set_ylabel('DSO voltage [V]')
ax.errorbar(x, y, dy, dx, 'ko', ms=1.2, elinewidth=0.8, capsize= 1.1,
        ls='',label='data', zorder=5)
ax.plot(x, y, 'gray', ls='-', lw=0.8, alpha = 0.8)
if fit: ax.plot(xx, sine(xx, *init), 'k--', lw=0.8, zorder=10, alpha =0.6,
                label='initial fit')
legend = ax.legend(loc ='best')
if tix:
    tick(ax, xmaj=5e-2, ymaj=20, ymin=5)
    if DSO: tick(ax, ymaj=0.1)

tran = np.fft.fftshift(np.fft.fft(y*np.kaiser(len(y), 4), fftsize))
freq = np.fft.fftshift(np.fft.fftfreq(fftsize, d=fres))
# plotfft(freq, tran, dB=False, re_im=False)
fft = tran.real if re_im else 20*np.log10(np.abs(tran)) if dB else np.abs(tran)
fft/=np.max(fft)
if angle or re_im:
    fig, (ax1, ax2) = plt.subplots(2,1, True, gridspec_kw={'wspace':0.05,
         'hspace':0.05})
else: fig, (ax2, ax1) = plt.subplots(2, 1, gridspec_kw={'wspace':0.25,
         'hspace':0.25})

ax1 = grid(ax1, xlab = 'Frequency $f$ [Hz]')
ax1.plot(freq, fft, c='k', lw='0.9')
ax1.set_ylabel('$\widetilde{V}(f)$ Magnitude [%s]' %('dB' if dB else 'arb. un.'))
if re_im: ax1.set_ylabel('Fourier Transform [Re]')    
ax1.set_xlim(0, 700)
if log: logy(ax1, 16)
elif tix:
    tick(ax1, xmaj=100, ymaj=1e3)
    if dB: tick(ax1, ymaj=20, ymin=5)

ax2 = grid(ax2, 'Time $t$ [s]', '$V(t)$ [digit]')
if DSO: ax2.set_ylabel('$V(t)$ [V]')
if angle or re_im: 
    fft = tran.imag if re_im else np.angle(tran)
    ax2.plot(freq, fft, c='k', lw='0.9')
    ax2.set_xlabel('Frequency $f$ [Hz]'); 
    ax2.set_ylabel('$\widetilde{V}(f)$ Phase [rad]')
    if re_im: ax2.set_ylabel('Fourier Transform [Im]')    
    if tix: tick(ax2, xmaj=50, ymaj=0.5)
else:
    ax2.errorbar(x, y, dy, dx, 'ko', ms=1.2, elinewidth=0.8, capsize= 1.1,
                 ls='',label='data', zorder=5)
    # ax2.set_xlim(0, 0.1)
    if tix: 
        tick(ax2, xmaj=1e-2, ymaj=20, ymin=5)
        if DSO: tick(ax2, ymaj=0.2, ymin=5e-2)

if not fit: plt.show(); exit()
#Fit V(t) con sinusoide
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
fig1, (ax1, ax2) = pltfitres(sx, sdx, sy, sdy, sine, pars)
ax1.set_ylabel('Differenza di potenziale $\Delta V$ [digit]')
legend = ax1.legend(loc ='lower right', framealpha = 0.3)
legend.set_zorder(100)
if log: ax1.set_yscale('log')
elif tix: tick(ax1, ymaj=20, ymin=5)

ax2.set_xlabel('Tempo $t$ [ms]', x=0.92); ax2.set_ylabel('Residui')
if tix: tick(ax2, xmaj=10, ymaj=5)

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
ax1.set_ylabel('Differenza di potenziale $\Delta V$ [digit]')
if log: ax1.set_xscale('log')
elif tix: tick(ax1, ymaj=20, ymin=5)
legend = ax1.legend(loc ='lower right', framealpha = 0.3)
legend.set_zorder(100)

ax2.errorbar(outT, normout, None, None, 'gx', elinewidth = 0.7, capsize=0.7,
             ms=3., ls='', zorder=5)
ax2.set_xlabel('Tempo $t$ [ms]', x=0.92); ax2.set_ylabel('Residui')
if tix:
    tick(ax2, xmaj=10, ymaj=2, ymin=0.5)
    ax2.set_ylim(np.min(normin)-np.std(normin), np.max(normin)+np.std(normin))
plt.show()