# -*- coding: utf-8 -*-
"""
Created on Sat Apr 11 23:19:50 2020

@author: berni
"""
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from scipy import signal as sg
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
    """ Dampened oscillation model function. Std argument order + tau"""
    return np.exp(-t/tau)*cosn(t, A, frq, phs) + ofs

def trg(t, A, frq, phs=0, ofs=0, duty=0.5):
    """ Triangular wave model function. Std argument order + duty cycle"""
    return A * sg.sawtooth(2*np.pi*frq*t + phs, duty) + ofs
    
def sqw(t, A, frq, phs=0, ofs=0, duty=0.5):
    """ Square wave model function. Std argument order + duty cycle"""
    return A * sg.square(2*np.pi*frq*t + phs, duty) + ofs

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

def LPF(data, gain=2e-4):
    """ Naive implementation of a digital Low Pass Filter (first order). """
    fltr = [data[0]]
    for x in data[1:]:
        fltr.append(fltr[-1] + gain*(x - fltr[-1]))
    return np.asarray(fltr)

def butf(signal, order, fc, ftype='lp', sampf=None):
    """
    Butterworth {lowpass, highpass, bandpass, bandstop} digital filter

    Parameters
    ----------
    signal : array_like
        Central values of the dependent measured signal over time.
    order : int
        Order of the butterworth filter, the higher the slower the response.
    fc : scalar or tuple
        Critical frequencies of the filter: High, Low cutoff (-3dB) frequencies.
    ftype : str, optional
        Type of filter {'LPF', 'HPF', 'BPF', 'BSF'}. The default is 'lp'.
    sampf : float, optional
        Sampling frequency of the system that measured the signal. If not
        given fc must be normalized to fc/fs. The default is None.

    Returns
    -------
    ndarray
        The digitally filtered signal.

    """
    butw = sg.butter(N=order, Wn=fc, btype=ftype, output='sos', fs=sampf)
    return sg.sosfilt(butw, signal)

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
    """ Computes parameter error and correlation matrix from covariance. """
    perr = np.sqrt(covm.diagonal())
    corm = np.array(covm, copy=True)
    for i in range(np.size(corm, 0)):
        for j in range(np.size(corm, 1)):
            corm[i][j] /= perr[i]*perr[j]
    return perr, corm

def prncor(corm, model):
    """ Prints formatted covariance of fit model args to precision %.3f """    
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
        
def RMS(seq):
    """ Evaluates Root Mean Square of data sequence through Numpy module. """
    return np.sqrt(np.mean(np.square(np.abs(seq))))

def RMSE(seq, exp=None):
    """ Associated Root Mean Square Error, expected value = RMS(seq) """
    if not exp: exp = RMS(seq)
    return np.sqrt(np.mean(np.square(seq - exp)))

def FWHM(x, y, FT=None):
    """ Evaluates FWHM of fundamental peak for y over dynamic variable x."""
    if FT: x=np.fft.fftshift(x); y=np.abs(np.fft.fftshift(y))
    d = y - (np.max(y) / 2.)
    indexes = np.where(d > 0)[0]
    if FT: indexes = [i for i in indexes if x[i] > 0]
    if len(indexes) >= 2:   
        return np.abs(x[indexes[-1]] - x[indexes[0]])

def optm(x, y, minm=None, absv=None):
    """ Evaluates minima or maxima (default) of y as a function of x. """
    x = np.asarray(x); y = np.abs(y) if absv else np.asarray(y) 
    yopt = np.min(y) if minm else np.max(y)
    xopt = np.where(y == yopt)[0]
    xopt = xopt[0]
    return x[xopt], yopt

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
    """ Log-scales x-axis, can add tix logarithmically spaced ticks to ax. """
    ax.set_xscale('log')
    if tix:
        ax.xaxis.set_major_locator(plt.LogLocator(numticks=tix))
        ax.xaxis.set_minor_locator(plt.LogLocator(subs=np.arange(2, 10)*.1,
                                                  numticks = tix))
def logy(ax, tix=None):
    """ Log-scales y-axis, can add tix logarithmically spaced ticks to ax. """
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
    dy: ndarray
        The propagated uncertainty of ymes after the iterated fit proces,
        the same dy passed as argument if propagation was not necessary.
    popt : ndarray
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

def pltfitres(xmes, dx, ymes, dy=None, model=None, pars=None, out=None):
# Variables that control the script 
    # kwargs.setdefault(
    #     {
    #     'DSO' : True, # Sampling from Digital Oscilloscope
    #     'fit' : False, # attempt to fit the data
    #     'log' : True, # log-scale axis/es
    #     'dB' : True, # plots response y-axis in deciBels
    #     'tick' : True, # manually choose spacing between axis ticks
    #     'tex' : True, # LaTeX typesetting maths and descriptions
    #     })
    fig, (ax1, ax2) = plt.subplots(2,1, True, gridspec_kw={
    'wspace':0.05, 'hspace':0.05, 'height_ratios': [3, 1]})
    space = np.linspace(np.min(0.9*xmes), np.max(1.1*xmes), 5000)
    if out  is not None: space = np.linspace(np.min(0.9*out), np.max(1.1*out), 5000)
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

# UTILITIES FOR FOURIER TRANSFORMS OF DATA ARRAYS
    
def FFT(time, signal, window=None, beta=0, specres=None):
    """
    Computes Discrete Fourier Transform of signal and its frequency space.

    Parameters
    ----------
    time : array_like
        Time interval over which the signal is sampled.
    signal : array_like
        The ADC sampled signal to be transformed with fft.
    window : function, optional
        Numpy window function with which to filter fft. The default is None.
    beta : int, optional
        Kaiser window's beta parameter. The default is 0, rect filter.
    specres : float, optional
        Sample spacing of the system that acquired signal. The default is None.

    Returns
    -------
    freq : ndarray
        Discrete sample frequency spectrum.
    tran : complex ndarray
        One-dimensional discrete Fourier Transform of the signal.

    """
    # Spectral resolution and number of points over which FFT is computed      
    if specres: fres = specres
    else: fres, frstd = sampling(space=time, dev=True, v=True)
    fftsize = len(time)
    
    if beta: window = lambda M: np.kaiser(M, beta=beta)
    elif window == np.kaiser: window = lambda M: np.kaiser(M, beta=0)
    tran = np.fft.rfft(signal*window(len(signal)), fftsize)
    freq = np.fft.rfftfreq(fftsize, d=fres)
    if not specres: return freq, tran, fres, frstd
    return freq, tran

def plotfft(freq, tran, signal=None, norm=False, dB=False, re_im=False, mod_ph=False):
    fft = tran.real if re_im  else np.abs(tran)
    if norm: fft/=np.max(fft)
    if dB: fft = 20*np.log10(np.abs(fft))
    if mod_ph or re_im:
        fig, (ax1, ax2) = plt.subplots(2,1, True, gridspec_kw={'wspace':0.05,
                                                               'hspace':0.05})
    else: fig, (ax2, ax1) = plt.subplots(2, 1, gridspec_kw={'wspace':0.25,
                                                            'hspace':0.25}) 
    ax1 = grid(ax1, xlab = 'Frequency $f$ [Hz]')
    ax1.plot(freq, fft, c='k', lw='0.9')
    ax1.set_ylabel('$\widetilde{V}(f)$ Magnitude [%s]' %('dB' if dB else 'arb. un.'))
    if re_im: ax1.set_ylabel('Fourier Transform [Re]')    

    ax2 = grid(ax2, 'Time $t$ [s]', '$V(t)$ [V]')
    if mod_ph or re_im: 
        fft = tran.imag if re_im else np.angle(tran)
        ax2.plot(freq, fft, c='k', lw='0.9')
        ax2.set_xlabel('Frequency $f$ [Hz]'); 
        ax2.set_ylabel('$\widetilde{V}(f)$ Phase [rad]')
        if re_im: ax2.set_ylabel('Fourier Transform [Im]')    
    else:
        xmes, dx, ymes, dy = signal
        ax2.errorbar(xmes, ymes, dy, dx, 'ko', ms=1.2, elinewidth=0.8,
                     capsize= 1.1, ls='-', lw=0.7, label='data', zorder=5)
    if signal: ax1, ax2 = [ax2, ax1]
    return fig, (ax1, ax2)
    
# UTILITIES FOR MANAGING DATA FILES
def srange(x, dx, y, dy, x_min=0, x_max=1e9):
    """ Extracts and returns the data inside selected range [min, max]. """ 
    xsup = x[x > x_min]; sx = xsup[xsup < x_max];
    ysup = y[x > x_min]; sy = ysup[xsup < x_max];
    dxsup = dx[x > x_min]; sdx = dxsup[xsup < x_max];
    dysup = dy[x > x_min]; sdy = dysup[xsup < x_max];
    return sx, sdx, sy, sdy
    
def sampling(space, dev=None, v=False):
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
    if dev: return Dxavg, Dxstd
    return Dxavg

def std_unc(measure, ADC=None):
    """ Associates default uncertainty to measured data array."""
    V = np.asarray(measure)
    if ADC: return np.ones_like(measure)
    unc = np.diff(V)/2/np.sqrt(12)
    return np.append(unc, np.mean(unc))

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