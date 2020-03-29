# -*- coding: utf-8 -*-
"""
Created on Sat Nov  2 17:19:42 2019

@author: berna
Visualizza dati con opportune bellurie
"""
import numpy as np
import cmath as c
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from scipy.misc import derivative
import glob
from sys import exit

# Estrazione e stampa a schermo delle medie delle misure di x
ddpfiles = glob.glob('../../plothist_data/avgdevs/*txt')
def digitddp():
    x=[]
    dx=[]
    for f in ddpfiles:
        #print(f)
        D = np.loadtxt(f, unpack = True, usecols=0)
        n = len(D)
        av = np.mean(D)
        dev = np.std(D, ddof=1)/np.sqrt(n)
        x.append(av)
        dx.append(dev)
        print('X = %f +- %f ' %(av, dev))
    return x, dx

''' Variables that control the script '''
DSO = True # Sampling from Digital Oscilloscope
fit = False # attempt to fit the data
log = True # log-scale axis/es
dB = True # plots response y-axis in deciBels
tick = True # manually choose spacing between axis ticks
tex = True # LaTeX typesetting maths and descriptions
real_imag = False # Whether to plot (Re, Im) or (Mag, Phs) of FFT
angle = False # Whether to plot phase of FFT or V(t)

# Extrazione dei vettori di dati grezzi
Dir = './phase_shift/'
V1, V2 = np.loadtxt(Dir+'tranb5_1.43.txt', unpack=True)#,
                    #skiprows=2*256, max_rows=250)

# Trasformazione dei dati nelle grandezze da fittare
x = V1/1e3
dx = np.full(len(x),(V1[1]-V1[2])/2)/1e3
y = V2 - np.mean(V2)
dy = np.full(len(y), 1)
if DSO:
    V1, V2 = np.genfromtxt(Dir+'DSO/DSO1_41.csv', float, delimiter=',',
                     skip_header = 2, usecols=(0, 1), unpack = True)
    x = V1
    dx = np.full(len(x),(V1[1]-V1[2])/2)
    y = V2 - np.mean(V2)
    dy = np.full(len(y), (V2[1]-V2[2])/2)

# Estrazione di un sottointervallo di dati
x_min = -0.3
x_max = 2.5
x1 = x[x>x_min]; sx = x1[x1<x_max];
y1 = y[x>x_min]; sy = y1[x1<x_max];
dx1 = dx[x>x_min]; sdx = dx1[x1<x_max];
dy1 = dy[x>x_min]; sdy = dy1[x1<x_max];

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

# Grafico preliminare dati
if tex:
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
init=(200, 2*np.pi*150, 0., 0.)
xx = np.linspace(min(x), max(x), 500)
fig, ax = plt.subplots()
ax.errorbar(sx, sy, sdy, sdx, 'ko', ms=1.2, elinewidth=0.8, capsize= 1.1,
        ls='',label='data', zorder=5)
ax.plot(sx, sy, 'gray', ls='-', lw=0.8, alpha = 0.8)
ax.grid(color = 'gray', ls = '--', alpha=0.7)
ax.set_xlabel('Time [ms]', x=0.9)
ax.set_ylabel('ADC voltage [digit]')
ax.minorticks_on()
ax.tick_params(direction='in', length=5, width=1., top=True, right=True)
ax.tick_params(which='minor', direction='in', width=1., top=True, right=True)

if tick:
    ax.xaxis.set_major_locator(plt.MultipleLocator(10))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(2))
    ax.yaxis.set_major_locator(plt.MultipleLocator(25))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(5))
    if DSO:
        ax.xaxis.set_major_locator(plt.MultipleLocator(5e-2))
        ax.xaxis.set_minor_locator(plt.MultipleLocator(1e-2))
        ax.yaxis.set_major_locator(plt.MultipleLocator(0.1))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(2e-2))
if fit:
    ax.plot(xx, sine(xx, *init), 'k--', lw=0.8, zorder=10, alpha =0.6,
            label='initial fit')
legend = ax.legend(loc ='best')

fig, (ax1, ax2) = plt.subplots(2,1, True, gridspec_kw={'wspace':0.05,
     'hspace':0.05})

""" Costanti per la trasformata discreta di Fourier """    
# elabora la media e deviazione standard degli intervalli di campionamento Dt
Dx=np.zeros(len(sx)-1)
xsort=np.sort(sx)

for i in range (0,len(Dx)):
    Dx[i]=xsort[i+1]-xsort[i]

Dxavg=np.average(Dx)
Dxstd=np.std(Dx)
print("Delta t average = %.e s" %Dxavg)
print("Delta t stdev = %.e s" %Dxstd)

# Lunghezza dell'array su cui calcola la trasformata
fftsize = 100000
fres = Dxavg
trans = np.fft.fft(sy*np.bartlett(len(sy)), fftsize)
freq = np.fft.fftfreq(fftsize, d=fres)
trans = np.fft.fftshift(trans)
freq = np.fft.fftshift(freq)

if real_imag:
    ax1.plot(freq, trans.real,  c='k', lw='0.9')
    ax1.set_ylabel('Fourier Transform [Re]')    
elif dB:
    log = False
    ax1.plot(freq, 20*np.log10(np.abs(trans)),  c='k', lw='0.9')
    ax1.set_ylabel('$\widetilde{V}(f)$ Magnitude [dB]')
else:
    ax1.plot(freq, np.abs(trans),  c='k', lw='0.9')
    ax1.set_ylabel('$\widetilde{V}(f)$ Magnitude [arb. un.]')
    

ax1.grid(color = 'gray', ls = '--', alpha=0.7)
ax1.minorticks_on()
ax1.tick_params(direction='in', length=5, width=1., top=True, right=True)
ax1.tick_params(which='minor', direction='in', width=1., top=True, right=True)

if log:
    ax1.set_yscale('log')
    if tick:
        ax1.yaxis.set_major_locator(plt.LogLocator(numticks=16))
        ax1.yaxis.set_minor_locator(plt.LogLocator(subs=np.arange(2, 10)*.1,
                                              numticks = 16))
elif tick:
    ax1.yaxis.set_major_locator(plt.MultipleLocator(1e3))
    ax1.yaxis.set_minor_locator(plt.MultipleLocator(2e2))
    if dB:
        ax1.yaxis.set_major_locator(plt.MultipleLocator(20))
        ax1.yaxis.set_minor_locator(plt.MultipleLocator(4))

if real_imag:
    ax2.plot(freq, trans.imag,  c='k', lw='0.9')
    ax2.set_ylabel('Fourier Transform [Im]')
else: 
    ax2.plot(freq, np.angle(trans), c='k', lw='0.9')
    ax2.set_ylabel('$\widetilde{V}(f)$ Phase [rad]')
    
ax2.grid(color = 'gray', ls = '--', alpha=0.7)
ax2.minorticks_on()
ax2.tick_params(direction='in', length=5, width=1., top=True, right=True)
ax2.tick_params(which='minor', direction='in', width=1., top=True, right=True)
ax2.set_xlabel('Frequency $f$ [Hz]', x=0.9)

if tick:
    ax2.xaxis.set_major_locator(plt.MultipleLocator(50))
    ax2.xaxis.set_minor_locator(plt.MultipleLocator(10))
    ax2.yaxis.set_major_locator(plt.MultipleLocator(0.5))
    ax2.yaxis.set_minor_locator(plt.MultipleLocator(0.1))
    
ax2.set_xlim(0, 1200)
plt.show()

if not fit:
    exit()
#Fit con parabola
pars, covm = curve_fit(mod, sx, sy, init, sdy, absolute_sigma = False)
A_fit, T_fit, f_fit, B_fit = pars
dA_fit, dT_fit, df_fit, dB_fit = np.sqrt(covm.diagonal())
print('Parametri del fit:\n', pars)
print('Matrice di Covarianza:\n', covm)
print('A = %f +- %f mV' % (A_fit, dA_fit))
print('T = %f +- %f ms' % (T_fit, dT_fit))
print('fase = %f +- %f ms' % (f_fit, df_fit))
print('B = %f +- %f mV' % (B_fit, dB_fit))
#Test Chi quadro per par
res = sy - mod(sx, *pars)
resnorm = res/sdy
chisq, ndof, sigma = chitest(sy, sdy, mod(sx, *pars), ddof=len(pars))
print('Chi quadro/ndof = %f/%d [%+.1f]' % (chisq, ndof, sigma))
print('Chi quadro ridotto_v:', chisq/ndof)
# Covarianza tra i parametri parabola
corr_AT = covm[0][1]/(dA_fit*dT_fit)
corr_Af = covm[0][2]/(dA_fit*df_fit)
corr_AB = covm[0][3]/(dA_fit*dB_fit)
corr_Tf = covm[1][2]/(dT_fit*df_fit)
corr_BT = covm[1][3]/(dB_fit*dT_fit)
corr_Bf = covm[2][3]/(dB_fit*df_fit)
corm = np.zeros((4,4))
for i in range(4):
    for j in range (4):
        corm[i][j] = covm[i][j]/covm[i][i]

print('Matrice di correlazione:\n', corm)    
print('Covarianza normalizzata AT:', corr_AT)
print('Covarianza normalizzata Af:', corr_Af)
print('Covarianza normalizzata AB:', corr_AB)
print('Covarianza normalizzata Tf:', corr_Tf)
print('Covarianza normalizzata BT:', corr_BT)
print('Covarianza normalizzata Bf:', corr_Bf)

#Plot DV vs t per par
if tex:
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

# Plot y vs x
xx = np.linspace(min(sx)-0.3, max(sx)+0.3, 500)
fig1,(ax1, ax2) = plt.subplots(2,1, True, gridspec_kw={'wspace':0.05,
     'hspace':0.05, 'height_ratios': [3, 1]})
if log:
    tick = False
    xx = np.logspace(np.log10(min(sx)), np.log10(max(sx)), len(xx))
    #ax1.set_yscale('log')
    ax1.set_xscale('log')
    ax1.minorticks_on()
ax1.set_ylabel('Differenza di potenziale $\Delta V$ [mV]')
ax1.grid(color = 'gray', ls = '--', alpha=0.7)
ax1.errorbar(sx, sy, sdy, sdx, 'ko', ms=1.5, elinewidth=1., capsize=1.5,
             ls='',label='data', zorder=5)
ax1.plot(xx, mod(xx, *pars), c='gray',
        label='parabola \n $\chi^2 = %.f/%d$' %(chisq, ndof), zorder=10, alpha=0.7)
if tick:
    ax1.yaxis.set_major_locator(plt.MultipleLocator(100.))
    ax1.yaxis.set_minor_locator(plt.MultipleLocator(20))
ax1.tick_params(direction='in', length=5, width=1., top=True, right=True)
ax1.tick_params(which='minor', direction='in', width=1., top=True, right=True)
legend = ax1.legend(loc ='lower right', framealpha = 0.3)
legend.set_zorder(100)

ax2.set_xlabel('Tempo $t$ [ms]', x=0.92)
ax2.set_ylabel('Residui')
ax2.axhline(0, c='r', zorder=10)
ax2.errorbar(sx, resnorm, None, None, 'ko', elinewidth = 0.7,
             capsize=0.7, ms=1., ls='', zorder=0)
ax2.grid(color ='gray', ls = '--', alpha=0.7)
ax2.ticklabel_format(axis='both', style='sci', scilimits=None,
                     useMathText=True)
if tick:
    ax2.xaxis.set_major_locator(plt.MultipleLocator(10))
    ax2.xaxis.set_minor_locator(plt.MultipleLocator(2))
    ax2.yaxis.set_major_locator(plt.MultipleLocator(2))
    ax2.yaxis.set_minor_locator(plt.MultipleLocator(0.5))
ax2.tick_params(direction='in', length=5, width=1., top=True, right=True)
ax2.tick_params(which='minor', direction='in', width=1., top=True, right=True)

# Compatibilità errori
soglia = 4
max_iter=100
Compatibili = True
if(np.mean(sdy) > soglia*(np.mean(sdx)*f_1(sdx, pars))):
    Compatibili = False
if(Compatibili):
    deff = sdy
    for n in range(max_iter):
        popt, pcov = curve_fit(mod, sx, sy, init, deff, absolute_sigma=False)
        m, q = popt
        dm, dq = np.sqrt(pcov.diagonal())
        deff=np.sqrt(deff**2 + (f_1(sx, popt) * sdx)**2)
        print(n, np.mean(deff) - soglia*np.mean(f_1(sx,popt)))
        if (np.mean(deff) > soglia*np.mean(f_1(sx,popt))):
            print(deff)
            print('Parametri ottimali:')
            print('m = %f +- %f' % (m, dm))
            print('q = %f +- %f compatibile con 0?' % (q, dq))
            # Tesx Chi quadro
            res = sy - mod(sx, *pars)
            resnorm = res/sdy
            chisq, ndof, sigma = chitest(sy, deff, mod(sx, *pars),
                                         ddof=len(pars))
            print('Chi quadro/ndof = %f/%d [%+.1f]' % (chisq, ndof, sigma))
            print('Chi quadro ridotto:', chisq/ndof)
            # Covarianza tra q e m
            corr = covm[0][1]/(dm*dq)
            corm = np.zeros((2,2))
            for i in range(2):
                for j in range (2):
                    corm[i][j] = covm[i][j]/covm[i][i]
    
            print('Covarianza normalizzata:', corr)
            print('Matrice di correlazione:\n', corm)
            break

# Fit parabole con rimozione degli outliers
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
#tengo solo i dati che si discostano dal modello per meno
# di 3 deviazioni standard
soglia = 3
for i in range (len(sy)):
    if (np.abs(sy[i] - mod(sx, *pars)[i])< soglia*sdy[i]): 
        TT=np.insert(TT, j, sx[i])
        dTT=np.insert(dTT, j, sdx[i])
        VV=np.insert(VV, j, sy[i])
        dVV=np.insert(dVV, j, sdy[i])
        j+=1
    else:
        outT=np.insert(outT, k, sx[i])
        doutT=np.insert(doutT, k, sdx[i])
        outV=np.insert(outV, k, sy[i])
        doutV=np.insert(doutV, k, sdy[i])
        k+=1

pars, covm = curve_fit(mod, TT, VV, init, dVV, absolute_sigma = False)
A_fit, omega_fit, f_fit, B_fit = pars
dA_fit, domega_fit, df_fit, dB_fit = np.sqrt(covm.diagonal())
print('Parametri del fit:\n', pars)
print('Matrice di Covarianza:\n', covm)
print('A = %f +- %f mV' % (A_fit, dA_fit))
print('T = %f +- %f ms' % (T_fit, dT_fit))
print('phi = %f +- %f rad' % (f_fit, df_fit))
print('B = %f +- %f mV' % (B_fit, dB_fit))

normin = (VV-mod(TT, *pars))/dVV
normout = (outV-mod(outT, *pars))/doutV
chisqin, ndof, sigmain = chitest(VV, dVV, mod(TT, *pars), ddof=len(pars))
print('Chi quadro/ndof = %f/%d [%+.1f]' % (chisqin, ndof, sigmain))
print('Chi quadro ridotto:', chisqin/ndof)
# Covarianza tra i parametri della parabola
corr_AT = covm[0][1]/(dA_fit*T_fit)
corr_Af = covm[0][2]/(dA_fit*df_fit)
corr_AB = covm[0][3]/(dA_fit*dB_fit)
corr_Tf = covm[1][2]/(dT_fit*df_fit)
corr_BT = covm[1][3]/(dB_fit*T_fit)
corr_Bf = covm[2][3]/(dB_fit*df_fit)
corm = np.zeros((4,4))
for i in range(4):
    for j in range (4):
        corm[i][j] = covm[i][j]/covm[i][i]
    
print('Matrice di correlazione:\n', corm)
print('Covarianza normalizzata AT:', corr_AT)
print('Covarianza normalizzata Af:', corr_Af)
print('Covarianza normalizzata AB:', corr_AB)
print('Covarianza normalizzata omegaf:', corr_Tf)
print('Covarianza normalizzata Bomega:', corr_BT)
print('Covarianza normalizzata Bf:', corr_Bf)

# Plot DV vs t con outliers di par
fig2,(ax1, ax2) = plt.subplots(2,1, True, gridspec_kw={'wspace':0.05,
     'hspace':0.05, 'height_ratios': [3, 1]})
if(log):
    ax1.set_xscale('log')
    ax1.minorticks_on()
ax1.errorbar(TT, VV, dVV, dTT, 'ko',  ms=1.5, elinewidth=1.,
             capsize=1.5, ls='', label='data')
ax1.errorbar(outT, outV, doutV, doutT, 'gx',  ms=3, elinewidth=1.,
             capsize=1.5, ls='', label='outliers')
ax1.set_ylabel('Differenza di potenziale $\Delta V$ [mV]')
ax1.grid(which = 'major', color = 'gray', ls = '--', alpha=0.7)
ax1.plot(xx, mod(xx, *pars), c='gray',
         label='fit \n $\chi^2 = %.f/%d$' %(chisqin, ndof), zorder=10, alpha=0.7)
if tick:
    ax1.yaxis.set_major_locator(plt.MultipleLocator(100))
    ax1.yaxis.set_minor_locator(plt.MultipleLocator(20))

ax1.tick_params(direction='in', length=5, width=1., top=True, right=True)
ax1.tick_params(which='minor', direction='in', width=1., top=True, right=True)
legend = ax1.legend(loc ='lower right', framealpha = 0.3)
legend.set_zorder(100)

ax2.set_xlabel('Tempo $t$ [ms]', x=0.92)
ax2.set_ylabel('Residui')
ax2.axhline(0, c='r', alpha=0.7, zorder=10)
ax2.errorbar(TT, normin, None, None, 'ko', elinewidth = 0.5, capsize=0.5,
             ms=1., ls='', lw=1., zorder=5)
ax2.errorbar(outT, normout, None, None, 'gx', elinewidth = 0.7, capsize=0.7,
             ms=3., ls='', zorder=5)
ax2.grid(color ='gray', ls = '--', alpha=0.7)
ax2.ticklabel_format(axis='both', style='sci', scilimits=None,
                     useMathText=True)
if tick:
    ax2.xaxis.set_major_locator(plt.MultipleLocator(10))
    ax2.xaxis.set_minor_locator(plt.MultipleLocator(2))
    ax2.yaxis.set_major_locator(plt.MultipleLocator(1))
    ax2.yaxis.set_minor_locator(plt.MultipleLocator(0.2))
    ax2.set_ylim(min(normin)-np.std(normin), max(normin)+np.std(normin))
ax2.tick_params(direction='in', length=5, width=1., top=True, right=True)
ax2.tick_params(which='minor', direction='in', width=1., top=True, right=True)


ax2.set_xlabel('Lettura digitalizzata [digit]', x=0.83)
ax2.set_ylabel('Residui')
ax2.axhline(0, c='r', zorder=10)
ax2.errorbar(sx, resnorm, None, None, 'ko', elinewidth=1, capsize=2, ms=2.5,
             ls='--', lw=1., zorder=0)
ax2.grid(color ='gray', ls = '--', alpha=0.7)
if tick:
    ax2.xaxis.set_major_locator(plt.MultipleLocator(abs(x[0]-x[len(x)-1])/8.))
    ax2.xaxis.set_minor_locator(plt.MultipleLocator(abs(x[0]-x[len(x)-1])/40.))
    ax2.yaxis.set_major_locator(plt.MultipleLocator(0.5))
    ax2.yaxis.set_minor_locator(plt.MultipleLocator(0.1))
ax2.tick_params(direction='in', length=5, width=1., top=True, right=True)
ax2.tick_params(which='minor', direction='in', width=1., top=True, right=True)

ax2.set_xlabel('Lettura digitalizzata [digit]', x=0.83)
ax2.set_ylabel('Residui')
ax2.axhline(0, c='r', zorder=10)
ax2.errorbar(sx, resnorm, None, None, 'ko', elinewidth=1, capsize=2, ms=2.5,
             ls='--', lw=1., zorder=0)
ax2.grid(color ='gray', ls = '--', alpha=0.7)
if tick:
    ax2.xaxis.set_major_locator(plt.MultipleLocator(abs(x[0]-x[len(x)-1])/8.))
    ax2.xaxis.set_minor_locator(plt.MultipleLocator(abs(x[0]-x[len(x)-1])/40.))
    ax2.yaxis.set_major_locator(plt.MultipleLocator(0.5))
    ax2.yaxis.set_minor_locator(plt.MultipleLocator(0.1))
ax2.tick_params(direction='in', length=5, width=1., top=True, right=True)
ax2.tick_params(which='minor', direction='in', width=1., top=True, right=True)