# -*- coding: utf-8 -*-
"""
Created on Sat Apr 11 23:19:50 2020

@author: berni
"""
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from scipy.misc import derivative
import glob

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

# Modello sinusoidale e parabolico
def sine(t, A, omega, f, B):
    return A*np.sin(omega*t + f) + B

def parabola(x, A , T, B):
    """ Definizione modulare della parabola integrale di un onda triangolare"""
    if (x < T/2.):
        return (A*x*((2*x/T)- 1) + B)
    else:
        return (A*(3*x -(2*x**2)/T - T) + B)

def mod(vect, A, T, f, B):
    y = []
    for x in vect:
        y.append(parabola((x+f)%T, A, T, B))
    return y

def f_1 (x, pars):
    return derivative(mod, x, dx = 1e-6, n = 1, args = pars);

def chitest(data, unc, model, ddof=0):
    """ Evaluates Chi-square goodness of fit test for a function, model, to
    a set of data """
    res = data - model
    resnorm = res/unc
    ndof = len(data) - ddof
    chisq = (resnorm**2).sum()
    sigma = (chisq - ndof)/np.sqrt(2*ndof)
    return chisq, ndof, sigma    

def fit(model, x, dx, y, dy, p0=None, absigma=False, **kwargs):
# Variables that control the script 
    kwargs.setdefault(
        {
        'DSO' : True, # Sampling from Digital Oscilloscope
        'fit' : False, # attempt to fit the data
        'log' : True, # log-scale axis/es
        'dB' : True, # plots response y-axis in deciBels
        'tick' : True, # manually choose spacing between axis ticks
        'tex' : True, # LaTeX typesetting maths and descriptions
        'real_imag' : False, # Whether to plot (Re, Im) or (Mag, Phs) of FFT
        'angle' : False # Whether to plot phase of FFT or V(t)
        })