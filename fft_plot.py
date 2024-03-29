# -*- coding: utf-8 -*-
"""
Created on Sat Nov  2 17:19:42 2019

@author: berni
Computes Fast Fourier Transforms and fits periodic signals.
"""
from lab import (np, plt, curve_fit, chitest, propfit, grid, errcor,
                 std_unc, prnpar, prncor, outlier, pltfitres, tick, logy,
                 FFT, plotfft, FWHM, sine, optm, cosn, dosc, args)

''' Variables that control the script '''
fit = True # attempt to fit the data to function
lock = False # Lock-in amplify data before FFT
tix = False # manually choose spacing between axis ticks
tex = True # LaTeX typesetting maths and descriptions
log = False # log-scale axis/es
dB = True # plots response y-axis in deciBels
if lock and fit or log and dB:
    raise Warning('Mutually exclusive options, choose either one.')

from data import x, dx, y, dy, sx, sdx, sy, sdy, DSO, AC
function = dosc
# Grafico preliminare dati
if tex: plt.rc('text', usetex=True); plt.rc('font', family='serif')
if lock: dx = std_unc(x); dy = std_unc(y)
xx = np.linspace(np.min(0.9*sx), np.max(1.05*sx), 2000)
fig, ax = plt.subplots()
grid(ax, xlab='Time $t$ [s]', ylab='ADC voltage [digit]')
if DSO: ax.set_xlabel('Time $t$ [s]'); ax.set_ylabel('DSO voltage [V]')
ax.errorbar(x, y, dy, dx, 'ko', ms=1.2, elinewidth=0.8, capsize= 1.1,
        ls='',label='data', zorder=5)
ax.plot(x, y, 'gray', ls='-', lw=0.8, alpha = 0.8)
if fit:
    init = (np.std(sy), 300, 1.5, round(np.mean(y)), 0.015)
    if DSO: init = (0.5, 44, 0., round(np.mean(y)))
    ax.plot(xx, function(xx, *init), 'k--', lw=0.8, zorder=10, alpha =0.6,
                label='initial fit')
legend = ax.legend(loc ='best')
if tix:
    tick(ax, xmaj=5e-2, ymaj=100, ymin=None)
    if DSO: tick(ax, ymaj=0.1)

if fit:
    #Fit V(t) con sinusoide
    pars, covm = curve_fit(function, sx, sy, init, sdy, absolute_sigma = False)
    pars, covm, sdy = propfit(function, sx, sy, sdx, sdy, p0=pars,
                             tail=3, rtol=0.3, max_iter=20, v=True)
    perr, pcor = errcor(covm)
    #print('Matrice di correlazione:\n', pcor) 
    #prncor(pcor, function)

    chisq, ndof, resnorm = chitest(function(sx, *pars), data=sy, unc=sdy,
                                   ddof=len(pars), v=True)
    print(f'Chi quadro ridotto: {chisq/ndof:.1f}' )
   
    #Plot DV vs t
    fig1, (ax1, ax2) = pltfitres(function, sx, sy, sdx, sdy, pars)
    ax1.set_ylabel('Differenza di potenziale $\Delta V$ [digit]')
    legend = ax1.legend(loc ='lower right', framealpha = 0.3)
    legend.set_zorder(100)
    #if log: ax1.set_yscale('log')
    if tix: tick(ax1, ymaj=50, ymin=None)

    ax2.set_xlabel('Tempo $t$ [ms]', x=0.92); ax2.set_ylabel('Residui')
    if tix: tick(ax2, xmaj=1e-2, ymaj=20)

    # Fit V(t) con rimozione degli outliers
    tin, vin, dtin, dvin, mask = outlier(function, sx, sy, sdx, sdy, pars,
                                         thr=5, mask=True)
    pars, covm = curve_fit(function, tin, vin, pars, dvin, absolute_sigma=False)
    perr, pcor = errcor(covm)
    print('Parametri del fit:\n', pars)
    prnpar(pars, perr, function)
    print('Matrice di correlazione:\n', pcor)
    prncor(pcor, function)
    chisqin, ndof, normin = chitest(function(tin, *pars), data=vin, unc=dvin,
                                    ddof=len(pars), v=True)
    print(f'Chi quadro ridotto: {chisqin/ndof:.1f}' )
    
    # Plot V(t) con outliers
    fig2, (ax1, ax2) = pltfitres(function, sx, sy, sdx, sdy, pars, in_out=mask)
    ax1.set_ylabel('Differenza di potenziale $\Delta V$ [digit]')
    #if log: ax1.set_yscale('log')
    if tix: tick(ax1, ymaj=50, ymin=None)
    legend = ax1.legend(loc ='lower right', framealpha = 0.3)
    legend.set_zorder(100)

    ax2.set_xlabel('Tempo $t$ [ms]', x=0.92); ax2.set_ylabel('Residui')
    if tix:
        tick(ax2, xmaj=1e-2, ymaj=5, ymin=None)
        ax2.set_ylim(np.min(normin)-np.std(normin), np.max(normin)+np.std(normin))

if lock:
    from lockin import t as time, Xmag as signal, DSO#, DAQ, Ve, Vpu, 
    dt = ds = None
else: time = sx; signal = sy; dt = sdx; ds = sdy
# FFT Computation with numpy window functions
freq, tran, fres, frstd = FFT(time, signal, window=np.kaiser, beta=None)
fmax, fftmax = optm(np.fft.fftshift(freq), np.fft.fftshift(tran), absv=True)
print(f"Fundamental frequency = {fmax:.2e} +- {2*fres/len(sy):.2e} peak at = {fftmax:.2e}" )
if AC and not lock: 
    fwhm = FWHM(freq, tran, FT=True)
    if fwhm: print("FWHM = %.1f Hz" %fwhm)
    Qf = np.sqrt(3)*fmax/fwhm; print(f'Qf = {Qf:.1f}')

# FFT and signal plot
fig, (ax1, ax2) = plotfft(freq, tran, signal=(time, signal, dt, ds),
                          norm=True, dB=dB, re_im=False)

if log: logy(ax2, 16)
elif tix:
    tick(ax2, xmaj=200, ymaj=1e3)
if dB: tick(ax2, ymaj=20, ymin=5); ax2.set_ylim(-100, 1)

if tix:
    tick(ax1, xmaj=5e-2, ymaj=100, ymin=None)
    if DSO:
        tick(ax1, xmaj=0.1, ymaj=0.2, ymin=5e-2)
    ax2.set_xlim(0, 3000)
    if lock:
        tick(ax1, xmaj=2e-2, ymaj=0.1, ymin=None)
    
if fit: 
    ax1.plot(xx, function(xx, *pars), c='gray', lw=1.2,
                 label='fit$\chi^2 = %.e/%d$' %(chisq, ndof))
    out = np.invert(mask)
    ax1.errorbar(sx[out], sy[out], sdy[out], sdx[out], 'gx',  ms=3,
                 elinewidth=1., capsize=1.5, ls='', label='outliers')
    C = 0.47e-6; L = 1./(C*((2*np.pi*pars[1])**2 + (1./pars[-1])**2))
    r = 2*L/pars[-1]
    print(f'Inductance: {L:.2f} H'); print(f'Circuit resistance: {r:.2f} ohm')
legend = ax1.legend(loc ='best', framealpha = 0.3)
plt.show()