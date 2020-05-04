# -*- coding: utf-8 -*-
"""
Created on Sat Apr 11 23:19:50 2020

@author: berni
"""
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from inspect import getfullargspec
from collections import namedtuple
import glob

# Standard argument order for periodical functions:
# A: Amplitude, frq: frequency, phs: phase, ofs: DC offset, tau: damp time
args=namedtuple('pars', 'A frq phs ofs')
 
def sine(t, A, frq, phs=0, ofs=0):
    """ Sinusoidal model function. Standard argument order """
    return A*np.sin(2*np.pi*frq*t + phs) + ofs
def cosn(t, A, frq, phs=0, ofs=0):
    return sine(t, A, frq, phs + np.pi/2., ofs)
def dosc(t, A, frq, phs, ofs, tau):
    """ Dampened oscillation model function. Std argument order"""
    return np.exp(-t/tau)*sine(t, A, frq, phs) + ofs

def parabola(x, A, T, phs=0, ofs=0):
    """ Modular definition of parabolae / integral of a triangle wave."""
    x +=phs
    if (x < T/2.):
        return (A*x*((2*x/T)- 1) + ofs)
    else:
        return (A*(3*x -(2*x**2)/T - T) + ofs)

def mod(vect, A, T, phs=0, ofs=0):
    """ Implementation of series of parabolas from integral of trg(x) """
    y = []
    for x in vect:
        y.append(parabola((x+phs)%T, A, T, ofs))
    return y

def fder(f, x, pars):
    return np.gradient(f(x, *pars), 1)

def LPF(data, gain=1e-3):
    fir = [data[0]]
    for x in data[1:]:
        fir.append(fir[-1] + gain*(x - fir[-1]))
    return np.asarray(fir)

# UTILITIES FOR MANAGING PARAMETER ESTIMATES AND TEST RESULTS
def chitest(data, unc, model, ddof=0, gauss=False, v=False):
    """ Evaluates Chi-square goodness of fit test for a function, model, to
    a set of data. """
    resn = (data - model)/unc
    ndof = len(data) - ddof
    chisq = (resn**2).sum()
    if gauss:
        sigma = (chisq - ndof)/np.sqrt(2*ndof)
        if v: print('Chi quadro/ndof = %.1f/%d [%+.1f]' % (chisq, ndof, sigma))
        return chisq, ndof, resn, sigma
    if v: print('Chi quadro/ndof = %.1f/%d' % (chisq, ndof))
    return chisq, ndof, resn

def errcor(covm):
    perr = np.sqrt(covm.diagonal())
    corm = np.asarray(covm)
    for i in range(np.size(corm, 0)):
        for j in range(np.size(corm, 1)):
            corm[i][j] /= perr[i]*perr[j]
    return perr, corm

def prncor(corm, model):
    pnms = getfullargspec(model)[0][1:]  
    for i in range(np.size(corm, 1)):
        for j in range(1 + i, np.size(corm, 1)):
            print('corr_%s-%s = %.3f' %(pnms[i], pnms[j], corm[i][j]))

def prnpar(pars, perr, model, prec=2):
    """ Prints formatted fit parameters to specified precision %.<prec>f """
    pnms = getfullargspec(model)[0][1:]
    dec = np.abs((np.log10(perr) - prec)).astype(int)
    for nam, par, err, d in zip(pnms, pars, perr, dec):
        print(f'{nam} = {par:.{d}f} +- {err:.{d}f}')        

# UTILITIES FOR MANAGING FIGURE AXES AND PLOTS
def grid(ax, xlab = None, ylab = None):
    """ Adds standard grid and labels for measured data plots to ax. """
    ax.grid(color = 'gray', ls = '--', alpha=0.7)
    ax.set_xlabel('%s' %(xlab if xlab else 'x [a.u.]'), x=0.9)
    ax.set_ylabel('%s' %(ylab if ylab else 'y [a.u.]'))
    ax.minorticks_on()
    ax.tick_params(direction='in', length=5, width=1., top=True, right=True)
    ax.tick_params(which='minor', direction='in', width=1., top=True, right=True)
    return ax
    
def tick(ax, xmaj = None, xmin = None, ymaj = None, ymin = None):
    """ Adds linearly spaced ticks to ax. """
    if not xmin: xmin = xmaj/5. if xmaj else None
    if not ymin: ymin = ymaj/5. if xmaj else None
    if xmaj: ax.xaxis.set_major_locator(plt.MultipleLocator(xmaj))
    if xmin: ax.xaxis.set_minor_locator(plt.MultipleLocator(xmin))
    if ymaj: ax.yaxis.set_major_locator(plt.MultipleLocator(ymaj))
    if ymin: ax.yaxis.set_minor_locator(plt.MultipleLocator(ymin))
    return ax

def logx(ax, tix=None):
    ax.set_xscale('log')
    if tix:
        ax.xaxis.set_major_locator(plt.LogLocator(numticks=tix))
        ax.xaxis.set_minor_locator(plt.LogLocator(subs=np.arange(2, 10)*.1,
                                                  numticks = tix))
def logy(ax, tix=None):
    ax.set_yscale('log')
    if tix:
        ax.yaxis.set_major_locator(plt.LogLocator(numticks=tix))
        ax.yaxis.set_minor_locator(plt.LogLocator(subs=np.arange(2, 10)*.1,
                                                  numticks = tix))

# LEAST SQUARE FITTING ROUTINES AND OUTPUT GRAPHS
# Compatibilità errori
def propfit(xmes, dx, ymes, dy, model, p0=None, thr=5, max_iter=20):
    """ Modified non-linear least squares (curve_fit) algorithm to
    fit model to a set of data, accounting for x-axis uncertainty
    using linear error propagation onto y-axis uncertainty.

    Parameters
    ----------
    xmes : array_like
        Central values of the independent variable where the data is measured.
    dx : array_like
        Uncertainties associated to xmes.
    ymes : array_like
        Central values of the dependent measured data.
    dy : array_like
        Uncertainties associated to ymes.
    model : callable
        The model function to be fitted to the data. 
    p0 : array_like, optional
        Initial guess for the parameters, assumes 1 if p0 is None.
    thr : float, optional
        Arbitrary constant, x-axis uncertainty is considered negligible
        and iterated fit complete whenever dy > thr * |f'(x)|dx. 5 by default.
    max_iter : int, optional
        Arbitrary natural constant, iteration is stopped if the thr condition
        is not met before max_iter iterations. 20 by default.

    Returns
    -------
    dy: array_like
        The propagated uncertainty of ymes after the iterated fit proces,
        the same dy passed as argument if propagation was not necessary.
    popt : array
        Optimal values for the parameters of model that mimimize chi square
        test, found by fitting with propagated uncertainties.
    pcov : 2darray
        The estimated covariance of popt.

    """
    deff = np.asarray(dy)
    for n in range(max_iter):
        popt, pcov = curve_fit(model, xmes, ymes, p0, deff, 
                               absolute_sigma=False)
        deff = np.sqrt(deff**2 + (dx * fder(model, xmes, popt))**2)
        # print(n, np.mean(deff) - thr*np.mean(dx * abs(fder(model, xmes,popt))))
        if np.mean(deff) > thr*np.mean(dx * abs(fder(model, xmes, popt))):
            perr, pcor = errcor(pcov)
            print('Parametri ottimali:')
            prnpar(popt, perr, model)
            # Test Chi quadro
            chisq, ndof, resn = chitest(ymes, deff, model(xmes, *popt),
                                         ddof=len(popt))          
            print('Chi quadro ridotto:', chisq/ndof)
            # Covarianza normalizzata
            print('Matrice di correlazione:\n', pcor)
            break
    return deff, popt, pcov

def outlier(xmes, dx, ymes, dy, model, pars, thr=5, out=False):
    """ Removes outliers from measured data. A sampled point is considered an
    outlier if it has absolute deviation y - model(x, *pars) > thr*dy. """
    isin = [np.abs(y - model(x, *pars)) < thr*sigma
            for x, y, sigma in zip(xmes, ymes, dy)]
    if out:
        isout = np.invert(isin)
        return (xmes[isin], dx[isin], ymes[isin], dy[isin], 
            xmes[isout], dx[isout], ymes[isout], dy[isout])
    
    return xmes[isin], dx[isin], ymes[isin], dy[isin]

def pltfitres(xmes, dx, ymes, dy=None, model=None, pars=None, **kwargs):
# Variables that control the script 
    # kwargs.setdefault(
    #     {
    #     'DSO' : True, # Sampling from Digital Oscilloscope
    #     'fit' : False, # attempt to fit the data
    #     'log' : True, # log-scale axis/es
    #     'dB' : True, # plots response y-axis in deciBels
    #     'tick' : True, # manually choose spacing between axis ticks
    #     'tex' : True, # LaTeX typesetting maths and descriptions
    #     'real_imag' : False, # Whether to plot (Re, Im) or (Mag, Phs) of FFT
    #     'angle' : False, # Whether to plot phase of FFT or V(t)
    #     })
    fig, (ax1, ax2) = plt.subplots(2,1, True, gridspec_kw={
    'wspace':0.05, 'hspace':0.05, 'height_ratios': [3, 1]})
    space = np.linspace(np.min(xmes-dx), np.max(xmes+dx), 2000)
    chisq, ndof, resn = chitest(ymes, dy, model(xmes, *pars), ddof=len(pars))
    ax1 = grid(ax1)
    ax1.errorbar(xmes, ymes, dy, dx, 'ko', ms=1.5, elinewidth=1., capsize=1.5,
             ls='', label='data')
    ax1.plot(space, model(space, *pars), c='gray', 
             label='fit$\chi^2 = %.1f/%d$' %(chisq, ndof))
    
    ax2 = grid(ax2)
    ax2.errorbar(xmes, resn, None , None, 'ko', elinewidth=0.5, capsize=1.,
             ms=1., ls='--', zorder=5)
    ax2.axhline(0, c='r', alpha=0.7, zorder=10)
    return fig, (ax1, ax2)

# UTILITIES FOR MANAGING DATA FILES
def srange(x, dx, y, dy, x_min=0, x_max=1e9):
    """ Extracts and returns the data inside selected range [min, max]. """ 
    xsup = x[x > x_min]; sx = xsup[xsup < x_max];
    ysup = y[x > x_min]; sy = ysup[xsup < x_max];
    dxsup = dx[x > x_min]; sdx = dxsup[xsup < x_max];
    dysup = dy[x > x_min]; sdy = dysup[xsup < x_max];
    return sx, sdx, sy, sdy
    
def sampling(space, v=False):
    """ Evaluates average sampling interval Dx and its standard deviation. """ 
    Deltas=np.zeros(len(space)-1)
    sort=np.sort(space)
    for i in range(len(Deltas)):
        Deltas[i] = sort[i+1] - sort[i]
    Dxavg=np.mean(Deltas)
    Dxstd=np.std(Deltas)
    if v:
        print(f"Delta t average = {Dxavg:.2e} s" )
        print(f"Delta t stdev = {Dxstd:.2e}")
    return Dxavg, Dxstd

# Estrazione e stampa a schermo delle medie delle misure di x
def digitddp():
    ddpfiles = glob.glob('../../plothist_data/avgdevs/*txt')
    for f in ddpfiles:
        #print(f)
        D = np.loadtxt(f, unpack = True, usecols=0)
        n = len(D)
        av = np.mean(D)
        dev = np.std(D, ddof=1)/np.sqrt(n)
        yield [av, dev]