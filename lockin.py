# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 17:46:05 2020
@author: berni
FFT by itself cannot resolve signals buried in noise of more than 30 times
greater amplitude. Simulates working principle of a lock in amplifier, given
enough time it can isolate signal phase and frequency from noise of
arbitrarily large scale through FFT of the demodulated argument, with improving
range from noise as LPF gain approaches 0 and window filters become narrower.
FFT of the modulus seems only to have noticeable harmonics for
gain <1e-3, seems to offer improvement on signal extraction from noisy
environment over regular FFT (SNR < 1/30). Experimentally best results in
extraction are obtained for gain = 1/(noise.scale * 100) and kaiser beta ~ 4. 
"""
from lab import np, plt, sine, cosn, grid, RMS, args, butf, sampling

''' Variables that control the script '''
log = False # log-scale axis/es
dB = False # plots response y-axis in deciBels
tix = False # manually choose spacing between axis ticks
tex = True # LaTeX typesetting maths and descriptions
DAQ = True # Use a real sampled singal or simulate it internally 

x_min = 0.; x_max = 1
t = np.linspace(x_min, x_max, 100_000)
if DAQ: ref = args(A=1, frq=167.72, phs=0, ofs=0)
else: ref = args(100, frq=1e2, phs=0, ofs=0)
sig = args(A=1, frq=1e2, phs=0.2, ofs=0)
cut = round(len(t)*5e-2) # Discarded initial points from LPFed components

Vsig = sine(t, *sig)
if DAQ: from data import x as t, y as Vsig, DSO
else: DSO = False
# fs = 1./sampling(t)
fs = len(t)/x_max
noise = np.random.normal(loc=0, scale=20, size=len(t))
print('Gaussian Noise: avg = %.2f std = %.2f' %(np.mean(noise), np.std(noise)))
if not DAQ: Vsig+=noise
Vsin = sine(t, *ref); Vcos = cosn(t, *ref)

Vsin = np.multiply(Vsin, Vsig); Vcos = np.multiply(Vcos, Vsig)
Y = 2*butf(Vcos, 1, fc=2, sampf=fs); X = 2*butf(Vsin, 1, fc=2, sampf=fs)
X /= ref.A; Y /= ref.A; Z = X + 1j*Y
Xmag = np.sqrt(X**2 + Y**2); Xphs = np.arctan2(Y, X);

# Plot
if tex: plt.rc('text', usetex=True); plt.rc('font', family='serif')
fig_sig, (ax1, ax2) = plt.subplots(2, 1, True, gridspec_kw={'wspace':0.05,'hspace':0.05})
grid(ax1, ylab='$V_{\\rm{sig}}(t)$')
ax1.errorbar(t, Vsig, ls='--', lw=0.9, c='k')
grid(ax2, xlab='Time $t$ [ms]', ylab='$V_{\\rm{ref}}(t)$')
ax2.plot(t, Vsin, lw=0.8, label='sin'); ax2.plot(t, Vcos, lw=0.8, label='cos')
leg_sig = ax2.legend(loc='best', framealpha=0.3)

# Discard initial transient of the digital filter > 1/10 of length
t = t[cut:]; Vsin = Vsin[cut:]; Vcos = Vcos[cut:]; X = X[cut:]; Y = Y[cut:]
Xmag = Xmag[cut:]; Xphs = Xphs[cut:]; Z = Z[cut:]
avg_mag = np.mean(Xmag); avg_phs = np.mean(Xphs); RMS_mag = RMS(Xmag)
std_mag = np.std(Xmag); std_phs = np.std(Xphs)

fig_mult, (ax1, ax2) = plt.subplots(2, 1, True, gridspec_kw={'wspace':0.05,'hspace':0.05})
grid(ax1, ylab='Magnitude $V(t)$')
ax1.plot(t, Xmag)
ax1.axhline(avg_mag + std_mag, c='orange', ls='--', lw=.9)
ax1.axhline(avg_mag, c='r', label='$V(t)$ Average Mag = %.2f' %avg_mag)
ax1.axhline(RMS_mag, c='g', lw=1, label='$V(t)$ RMS = %.2f' %RMS_mag)
ax1.axhline(avg_mag - std_mag, c='orange', lw=.9, ls='--', label='Mag std. dev = %.2f'%std_mag)

grid(ax2, xlab ='Time $t$ [ms]', ylab='Phase $\phi(t)$')
ax2.plot(t, Xphs)
ax2.axhline(avg_phs + std_phs, c='orange', ls='--', lw=.9)
ax2.axhline(avg_phs, c='r', label='$V(t)$ Average Phase = %.2f' %avg_phs)
ax2.axhline(avg_phs - std_phs, c='orange', lw=.9, ls='--', label='phs std. dev = %.2f'%std_phs)
leg1 = ax1.legend(loc ='best', framealpha = 0.3); leg2 = ax2.legend(loc ='best', framealpha = 0.3)
plt.show()