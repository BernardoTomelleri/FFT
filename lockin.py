# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 17:46:05 2020
@author: berni
FFT by itself cannot resolve signals buried in noise of more than 30 times
greater amplitude. Simulates working principle of a lock in amplifier, given
enough time it can isolate signal from noise of arbitrarily large scale. 
"""
import numpy as np
from matplotlib import pyplot as plt
import scipy.signal as sg
from lab import srange, sine, cosn, grid, LPF, args, sampling

''' Variables that control the script '''
log = False # log-scale axis/es
dB = False # plots response y-axis in deciBels
tix = False # manually choose spacing between axis ticks
tex = True # LaTeX typesetting maths and descriptions
DAQ = False # Use a real sampled singal or simulated

x_min = 0.; x_max = 2
t = np.linspace(x_min, x_max, 50_000)
cut = round(len(t)*1e-1) # Discarded initial points from LPFed components 

ref = args(A=10, frq=1e2, phs=0, ofs=0)
sig = args(A=1, frq=1e2, phs=0.2, ofs=0)
Vsin = sine(t, *ref); Vcos = cosn(t, *ref)
noise = np.random.normal(loc=0, scale=1000, size=len(t))
Vsig = cosn(t, *sig) + (noise - np.mean(noise))
Vsin = np.multiply(Vsin, Vsig); Vcos = np.multiply(Vcos, Vsig)
X = LPF(Vcos)/(ref.A/2); Y = LPF(Vsin)/(ref.A/2)
Z = X + 1j*Y
Xmag = np.sqrt(X**2 + Y**2); Xphs = np.arctan2(Y, X)
avg_mag = np.mean(Xmag); avg_phs = np.mean(Xphs)
std_mag = np.std(Xmag); std_phs = np.std(Xphs)

# Plot
if tex: plt.rc('text', usetex=True); plt.rc('font', family='serif')
fig_sig, (ax1, ax2) = plt.subplots(2, 1, True, gridspec_kw={'wspace':0.05,'hspace':0.05})
grid(ax1, ylab='$V_{\\rm{sig}}(t)$')
ax1.errorbar(t, Vsig, ls='--', lw=0.9, c='k')
grid(ax2, xlab='Time $t$ [ms]', ylab='$V_{\\rm{ref}}(t)$')
ax2.errorbar(t, Vsin, lw=0.9, label='sin'); ax2.plot(t, Vcos, lw=0.9, label='cos')
leg_sig = ax2.legend(loc='best', framealpha=0.3)

t = t[cut:]; Vsin = Vsin[cut:]; Vcos = Vcos[cut:]; X = X[cut:]; Y = Y[cut:]
Xmag = Xmag[cut:]; Xphs = Xphs[cut:];
fig_mult, (ax1, ax2) = plt.subplots(2, 1, True, gridspec_kw={'wspace':0.05,'hspace':0.05})
grid(ax1, ylab='Magnitude $V(t)$')
ax1.plot(t, Xmag)
ax1.axhline(avg_mag + std_mag, c='orange', ls='--', lw=.9)
ax1.axhline(avg_mag, c='r', label='$V(t)$ Average Mag = %.2f' %avg_mag)
ax1.axhline(avg_mag - std_mag, c='orange', lw=.9, ls='--', label='Mag std. dev = %.2f'%std_mag)

grid(ax2, xlab ='Time $t$ [ms]', ylab='Phase $\phi(t)$')
ax2.plot(t, Xphs)
ax2.axhline(avg_phs + std_phs, c='orange', ls='--', lw=.9)
ax2.axhline(avg_phs, c='r', label='$V(t)$ Average Phase = %.2f' %avg_phs)
ax2.axhline(avg_phs - std_phs, c='orange', lw=.9, ls='--', label='phs std. dev = %.2f'%std_phs)
leg1 = ax1.legend(loc ='best', framealpha = 0.3); leg2 = ax2.legend(loc ='best', framealpha = 0.3)
plt.show()
