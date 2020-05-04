# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 17:46:05 2020
@author: berni
"""
import numpy as np
from matplotlib import pyplot as plt
from lab import srange, sine, cosn, grid, LPF, args

''' Variables that control the script '''
log = False # log-scale axis/es
dB = False # plots response y-axis in deciBels
tix = False # manually choose spacing between axis ticks
tex = True # LaTeX typesetting maths and descriptions

x_min = 0.; x_max = 30000
t = np.linspace(x_min, x_max, 100_000)
ref = args(A=1000, frq=1000, phs=0, ofs=0)
sig = args(A=1, frq=1000, phs=1, ofs=0)
Vsin = sine(t, *ref); Vcos = cosn(t, *ref)
Vsig = sine(t, *sig) + np.random.normal(loc=0, scale=10, size=len(t))
Xsin = LPF(np.multiply(Vsin, Vsig))/ref.A; Xcos = LPF(np.multiply(Vcos, Vsig))/ref.A 
Xmag = np.sqrt(Xsin**2 + Xcos**2); Xphs = np.arctan2(Xcos, Xsin)

# Plot
if tex: plt.rc('text', usetex=True); plt.rc('font', family='serif')
fig, (ax1, ax2) = plt.subplots(2, 1, True, gridspec_kw={'wspace':0.05,'hspace':0.05})
grid(ax1, ylab='$V_{\\rm{sig}}(t)$')
ax1.errorbar(t, Vsig, ls='--', lw=0.9, c='k')
grid(ax2, xlab='Time $t$ [ms]', ylab='$V_{\\rm{ref}}(t)$')
ax2.errorbar(t, Vsin, lw=0.9, label='sin'); ax2.plot(t, Vcos, lw=0.9, label='cos')
legend = ax2.legend(loc='best', framealpha=0.3)

fig, (ax1, ax2) = plt.subplots(2, 1, True, gridspec_kw={'wspace':0.05,'hspace':0.05})
grid(ax1, ylab='Magnitude $V(t)$')
ax1.plot(t, LPF(Xmag, gain=1e-2))
ax1.axhline(np.mean(LPF(Xmag, gain=1e-2)),c='orange', label='$V(t)$ Average Mag')
grid(ax2, xlab ='Time $t$ [ms]', ylab='Phase $\phi(t)$')
ax2.plot(t, LPF(Xphs, gain=1e-2))
ax2.axhline(np.mean(LPF(Xphs, gain=1e-2)), c='orange', label='$V(t)$ Average Phase')
leg1 = ax1.legend(loc ='best', framealpha = 0.3); leg2 = ax2.legend(loc ='best', framealpha = 0.3)
plt.show()