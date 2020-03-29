# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 17:46:05 2020

@author: berni
"""
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from sys import exit
""" Estrazione dei dati corrispondenti all'intervallo di ddp selezionato
tempo in microsecondi, allora quando lo dividi per 10^3 lo trasformi 
in millisecondi, dio cane. differenza di potenziale in unità arbitrarie """
Dir = './RLC_data/'
t, DV  = np.loadtxt(Dir+'harm0.1uF_b5.txt', unpack=True)
dt = np.ones_like(t)
dDV = np.ones_like(DV)
t/=10.**3
dt/=10.**3
t_min = 0.
t_max = 45.
t1 = t[t>t_min]; st = t1[t1<t_max];
dt1 = dt[t>t_min]; sdt = dt1[t1<t_max];
DV1 = DV[t>t_min]; sDV = DV1[t1<t_max];
dDV1 = dDV[t>t_min]; sdDV = dDV1[t1<t_max];

# Fit oscillazioni RLC smorzate di ddp nel tempo
def osc(t, A, tau, omega, f, B):
    return A*np.exp(-t/tau)*np.sin(omega*t +f) + B

def damp(t, A, tau, B):
    return A*np.exp(-t/tau) + B

def chitest(data, unc, model, ddof=0):
    """ Evaluates Chi-square goodness of fit test for a function, model, to
    a set of data """
    res = data - model
    resnorm = res/unc
    ndof = len(data) - ddof
    chisq = (resnorm**2).sum()
    sigma = (chisq - ndof)/np.sqrt(2*ndof)
    return chisq, ndof, sigma    

''' Variables that control the script '''
DSO = False # Sampling from Digital Oscilloscope
fit = False # attempt to fit the data
log = False # log-scale axis/es
tick = True # manually choose spacing between axis ticks
tex = True # LaTeX typesetting maths and descriptions

init = (900., 4., 2., 0., 450.)
pars, covm = curve_fit(osc, st, sDV, init, sdDV, absolute_sigma = True)
A_fit, tau_fit, omega_fit, f_fit, B_fit = pars
dA_fit, dtau_fit, domega_fit, df_fit, dB_fit = np.sqrt(covm.diagonal())
print('Parametri del fit:\n', pars)
print('Matrice di Covarianza:\n', covm)
print('A = %f +- %f a.u.' % (A_fit, dA_fit))
print('tau = %f +- %f s' % (tau_fit, dtau_fit))
print('omega = %f +- %f rad/s' % (omega_fit, domega_fit))
print('f = %f +- %f rad' %(f_fit, df_fit))
print('B = %f +- %f a.u.' % (B_fit, dB_fit))
res = sDV - osc(st, A_fit, tau_fit, omega_fit, f_fit, B_fit)
resnorm = res/sdDV
chisq, ndof, sigma = chitest(sDV, sdDV, osc(st, *pars), len(pars))
sigma = (chisq - ndof)/np.sqrt(2*ndof)
print('Chi quadro/ndof = %f/%d [%+.1f]' % (chisq, ndof, sigma))
print('Chi quadro ridotto:', (chisq/ndof))

# Covarianza tra i parametri
corr_Atau = covm[0][1]/(dA_fit*dtau_fit)
corr_Aomega = covm[0][2]/(dA_fit*domega_fit)
corr_Af = covm[0][3]/(dA_fit*df_fit)
corr_AB = covm[0][4]/(dA_fit*dB_fit)
corr_tauomega = covm[1][2]/(domega_fit*dtau_fit)
corr_tauf = covm[1][3]/(dtau_fit*df_fit)
corr_tauB = covm[1][4]/(dtau_fit*dB_fit)
corr_omegaf = covm[2][3]/(domega_fit*df_fit)
corr_omegaB = covm[2][4]/(domega_fit*dB_fit)
corr_fB = covm[3][4]/(dB_fit*df_fit)
corm = np.zeros((5,5))
for i in range(5):
    for j in range (5):
        corm[i][j] = covm[i][j]/covm[i][i]

print('Matrice di correlazione:\n', corm)    
print('Covarianza normalizzata Atau:', corr_Atau)
print('Covarianza normalizzata Aomega:', corr_Aomega)
print('Covarianza normalizzata Af:', corr_Af)
print('Covarianza normalizzata AB:', corr_AB)
print('Covarianza normalizzata tauomega:', corr_tauomega)
print('Covarianza normalizzata tauf:', corr_tauf)
print('Covarianza normalizzata tauB:', corr_tauB)
print('Covarianza normalizzata omegaf:', corr_omegaf)
print('Covarianza normalizzata omegaB:', corr_omegaB)
print('Covarianza normalizzata fB:', corr_fB)

# Plot DV vs t con esponenziali al contorno
tt = np.linspace(min(st)-0.05, max(st)+0.05, 5000)
if tex:
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    
fig,(ax1, ax2) = plt.subplots(2,1, True, gridspec_kw={
    'wspace':0.05, 'hspace':0.05, 'height_ratios': [3, 1]})
ax1.set_ylabel('Differenza di potenziale $\Delta V$ [a.u.]')
ax1.grid(color ='gray', linestyle = '--', alpha = 0.7)
ax1.errorbar(st, sDV, sdDV, sdt, 'ko', ms=1.5, elinewidth=1., capsize=1.5,
             ls='', label='data')
ax1.plot(tt, osc(tt, A_fit, tau_fit, omega_fit, f_fit, B_fit), c='gray',
         label='fit$\chi^2 = %.1f/%d$' %(chisq, ndof))
if fit:   
    ax1.plot(tt, damp(tt, A_fit, tau_fit, B_fit), c='b')
    ax1.plot(tt, -damp(tt, A_fit, tau_fit, -B_fit), c='b')
if tick:
    ax1.yaxis.set_major_locator(plt.MultipleLocator(100.))
    ax1.yaxis.set_minor_locator(plt.MultipleLocator(20.))
ax1.tick_params(direction='in', length=5, width=1., top=True, right=True)
ax1.tick_params(which='minor', direction='in', width=1., top=True, right=True)
legend = ax1.legend(loc='best')

ax2.plot(st, 0*sDV, 'r', zorder =10)
ax2.errorbar(st, resnorm, None , None, 'ko', elinewidth = 0.5, capsize=1.,
             ms=1., ls='', zorder = 5)
ax2.grid(color ='gray', ls = '--', alpha = 0.7)
ax2.set_xlabel('Tempo $t$ [ms]', x=0.92)
ax2.set_ylabel('Residui [a.u.]')
if tick:
    ax2.xaxis.set_major_locator(plt.MultipleLocator(1.))
    ax2.xaxis.set_minor_locator(plt.MultipleLocator(0.2))
    ax2.yaxis.set_major_locator(plt.MultipleLocator(5.))
    ax2.yaxis.set_minor_locator(plt.MultipleLocator(1.))
ax2.tick_params(direction='in', length=5, width=1., top=True, right=True)
ax2.tick_params(which='minor', direction='in', width=1., top=True, right=True)

# Fit con rimozione degli outliers
TT=np.array([])
VV=np.array([])
dTT=np.array([])
dVV=np.array([])

outT=np.array([])
outV=np.array([])
doutT=np.array([])
doutV=np.array([])
'Contatori in e out liers'
j=0
k=0
# tengo solo i dati che si discostano dalla funzione di fit per meno di 2.5
# deviazioni standard
soglia = 3.5
for i in range (len(sDV)):
    if (np.abs(sDV[i] - osc(st, *pars)[i])< soglia*sdDV[i]):
        TT=np.insert(TT, j, st[i])
        dTT=np.insert(dTT, j, sdt[i])
        VV=np.insert(VV, j, sDV[i])
        dVV=np.insert(dVV, j, sdDV[i])
        j+=1
    else:
        outT=np.insert(outT, k, st[i])
        doutT=np.insert(doutT, k, dt[i])
        outV=np.insert(outV, k, sDV[i])
        doutV=np.insert(doutV, k, dDV[i])
        k+=1

pars, covm = curve_fit(osc, TT, VV, pars, dVV, absolute_sigma = True)
A_fit, tau_fit, omega_fit, f_fit, B_fit = pars
dA_fit, dtau_fit, domega_fit, df_fit, dB_fit = np.sqrt(covm.diagonal())
print('Parametri del fit:\n', pars)
print('Matrice di Covarianza:\n', covm)
print('A = %f +- %f a.u.' % (A_fit, dA_fit))
print('tau = %f +- %f s' % (tau_fit, dtau_fit))
print('omega = %f +- %f rad/s' % (omega_fit, domega_fit))
print('f = %f +- %f rad/s' %(f_fit, df_fit))
print('B = %f +- %f a.u.' % (B_fit, dB_fit))
normin = (VV-osc(TT, *pars))/dVV
normout = (outV-osc(outT, *pars))/doutV
chisqin, ndof, sigmain = chitest(VV, dVV, osc(TT, *pars), ddof=len(pars))
print('Chi quadro/ndof = %f/%d [%+.1f]' % (chisqin, ndof, sigmain))
print('Chi quadro ridotto:', chisqin/ndof)
# Covarianza tra i parametri
corr_Atau = covm[0][1]/(dA_fit*dtau_fit)
corr_Aomega = covm[0][2]/(dA_fit*domega_fit)
corr_Af = covm[0][3]/(dA_fit*df_fit)
corr_AB = covm[0][4]/(dA_fit*dB_fit)
corr_tauomega = covm[1][2]/(domega_fit*dtau_fit)
corr_tauf = covm[1][3]/(dtau_fit*df_fit)
corr_tauB = covm[1][4]/(dtau_fit*dB_fit)
corr_omegaf = covm[2][3]/(domega_fit*df_fit)
corr_omegaB = covm[2][4]/(domega_fit*dB_fit)
corr_fB = covm[3][4]/(dB_fit*df_fit)
corm = np.zeros((5,5))
for i in range(5):
    for j in range (5):
        corm[i][j] = covm[i][j]/covm[i][i]

print('Matrice di correlazione:\n', corm)    
print('Covarianza normalizzata Atau:', corr_Atau)
print('Covarianza normalizzata Aomega:', corr_Aomega)
print('Covarianza normalizzata Af:', corr_Af)
print('Covarianza normalizzata AB:', corr_AB)
print('Covarianza normalizzata tauomega:', corr_tauomega)
print('Covarianza normalizzata tauf:', corr_tauf)
print('Covarianza normalizzata tauB:', corr_tauB)
print('Covarianza normalizzata omegaf:', corr_omegaf)
print('Covarianza normalizzata omegaB:', corr_omegaB)
print('Covarianza normalizzata fB:', corr_fB)

fig2,(ax1, ax2) = plt.subplots(2,1, True, gridspec_kw={
    'wspace':0.05, 'hspace':0.05, 'height_ratios': [3, 1]})
ax1.errorbar(TT, VV, dVV, dTT, 'ko',  ms=1.5, elinewidth=1., capsize=1.5,
             ls='', label='data')
ax1.errorbar(outT, outV, doutV, doutT, 'gx',  ms=3.5, elinewidth=1.,
             capsize=1.5, ls='', label='outliers')
ax1.set_ylabel('Differenza di potenziale $\Delta V$ [a.u.]')
ax1.grid(color = 'gray', ls = '--', alpha=0.7)
ax1.plot(tt, osc(tt, A_fit, tau_fit, omega_fit, f_fit, B_fit), c='gray',
         label='fit $\chi^2 = %.1f/%d$' %(chisqin, ndof), zorder=10, alpha=0.7)
if tick:
    ax1.yaxis.set_major_locator(plt.MultipleLocator(100.))
    ax1.yaxis.set_minor_locator(plt.MultipleLocator(20.))
ax1.tick_params(direction='in', length=5, width=1., top=True, right=True)
ax1.tick_params(which='minor', direction='in', width=1., top=True, right=True)
legend = ax1.legend(loc ='best')

ax2.set_xlabel('Tempo $t$ [ms]', x=0.92)
ax2.set_ylabel('Residui')
ax2.axhline(0, c='r', alpha=0.7, zorder=10)
ax2.errorbar(TT, normin, None, None, 'ko', elinewidth = 0.7, capsize=0.7, ms=1.,
             ls='--', lw=1., zorder=5)
ax2.errorbar(outT, normout, None, None, 'gx', elinewidth = 0.7, capsize=0.7,
             ms=3., ls='', zorder=5)
ax2.grid(color ='gray', ls = '--', alpha=0.7)
if tick:
    ax2.ticklabel_format(axis='both', style='sci', scilimits=None,
                         useMathText=True)
    ax2.xaxis.set_major_locator(plt.MultipleLocator(1))
    ax2.xaxis.set_minor_locator(plt.MultipleLocator(0.2))
    ax2.yaxis.set_major_locator(plt.MultipleLocator(1.))
    ax2.yaxis.set_minor_locator(plt.MultipleLocator(0.2))
    ax2.set_ylim(min(normin)-np.std(normin), max(normin)+np.std(normin))
ax2.tick_params(direction='in', length=5, width=1., top=True, right=True)
ax2.tick_params(which='minor', direction='in', width=1., top=True, right=True)
plt.show()