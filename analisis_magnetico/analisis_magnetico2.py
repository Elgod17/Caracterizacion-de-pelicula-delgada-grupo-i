import numpy as np
import matplotlib
matplotlib.use('TkAgg')  # o prueba 'Agg' si no vas a mostrar ventanas
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
from scipy.constants import Boltzmann
from scipy.signal import medfilt
with open('magnetismo.txt', 'r') as file:
    content = file.read()
    lineas = content.split('\n')

pos = lineas.index('[Data]')
headers = lineas[pos+1].split(',')
dic = {x:[] for x in headers}

for i in range(pos+2, len(lineas)):
    data = lineas[i].split(',')
    for j in range(len(data)):
        dic[headers[j]].append(data[j])


### 1. Graficar Moment (emu) vs Magnetic Field (Oe)

## extraer y limpiar datos 
mu = np.array(dic['Moment (emu)'])
H = np.array(dic['Magnetic Field (Oe)'])
pos_empty = np.where(mu == '')[0]
mu = np.delete(mu, pos_empty).astype(float)
H = np.delete(H, pos_empty).astype(float)
mask_pos = H > 0
mask_neg = H < 0

mu_f = np.zeros_like(mu)
mu_f[mask_pos] = medfilt(mu[mask_pos], 3)
mu_f[mask_neg] = medfilt(mu[mask_neg], 3)




#plt.figure(figsize = (10,8))
#plt.plot(H, mu_f, marker='o', linestyle='-', color='b', markersize=5, label = 'Filtrado')
#plt.plot(H, mu, marker='o', linestyle='-', color='r', markersize=5, label = 'No filtrado')
#plt.legend()
#plt.show()







temp_ = np.delete(np.array(dic['Temperature (K)']), pos_empty).astype(float)

x_ = H / temp_


def linear(x, m, b):
    return m*x + b
## X es la susceptibilidad magnetica
# susceptibilidad magnetica en la rama positiva de H
popt1, pcov1 = curve_fit(linear, H[mask_pos], mu_f[mask_pos])
X1, b_fit1 = popt1
errX1, errb1 = np.sqrt(np.diag(pcov1))
# en la rama negativa de H
popt, pcov = curve_fit(linear, H[mask_neg], mu_f[mask_neg])
X, b_fit = popt
errX, errb = np.sqrt(np.diag(pcov))

#print(f'Susceptibilidad magnetica rama positiva: X1 = {X1:.20f} ± {errX1:.20f}')
#print(f'Susceptibilidad magnetica rama negativa: X = {X:.20f} ± {errX:.20f}')



### momento magnetico del aluminio

mu_Al = np.zeros_like(mu)
mu_Al[mask_neg] = mu_f[mask_neg] - X*H[mask_neg]
mu_Al[mask_pos] = mu_f[mask_pos] - X1*H[mask_pos] 


#mu_Al = np.array(dic['Al_moment (emu)'])
#plt.figure(figsize = (10,8))
#plt.plot(x_, mu_Al, '.', markersize=3)
#plt.title('Magnetic Moment vs Magnetic Field')
#plt.xlabel('Magnetic Field (Oe)')
#plt.ylabel('Magnetic Moment Al (emu)')
#plt.grid()
#plt.savefig('moment_vs_field_Al.png')
#plt.show()
#plt.close()
## Comparación: Modelo Langevin o Brillouin

def coth(x):
    val = np.sinh(x)
    for j in val: 
        if j == 0:
            return ValueError("coth undefined")
    else: 
        return np.cosh(x)/np.sinh(x)

def ajuste_Brillouin(x,s,a,b,c):
    T1 = 1 + 1 / (2*s)
    T2 = coth(T1*a*x)
    T3 = -1 /(2*s) * coth(a*x / (2*s))
    return b * (T1*T2+T3) + c



def ajuste_Langevin(x,a,b,c):
    return b*(coth(a*x) - 1/(a*x)) + c

###################### TODO BIEN HASTA AQUI ########################################################
####################################################################################################
####################################################################################################
x_fit = np.linspace(min(x_), max(x_), 100)
popt_B, pcov_B = curve_fit(ajuste_Brillouin, x_, mu_Al)
s, mod_B, amp_B, corte_B = popt_B
serr, moderr_B, amperr_B, corterr_B = np.sqrt(np.diag(pcov_B))
brillouin_fit = ajuste_Brillouin(x_fit, s, mod_B, amp_B, corte_B)



popt_L, pcov_L = curve_fit(ajuste_Langevin, x_, mu_Al)
mod_L, amp_L, corte_L = popt_L
moderr_L, amperr_L, corterr_L = np.sqrt(np.diag(pcov_L))
langevin_fit = ajuste_Langevin(x_fit, moderr_L, amp_L, corte_L)

#plt.figure(figsize=(10,8))
#plt.plot(x_, mu_Al, 'o', label='Datos experimentales filtrados', markersize=5)
##plt.plot(x_fit, langevin_fit, '-', label=f'Langevin fit: amp={amp_L:.4}±{amperr_L:.4}, corte={corte_L:.4}±{corterr_L:.4}, , mod={mod_L:.4}±{moderr_L:.4}', linewidth=2)
##plt.plot(x_fit, brillouin_fit, '-', label=f'Brillouin fit: s={s:.2}±{serr:.2}, amp={amp_B:.4}±{amperr_B:.4}, corte={corte_B:.4}±{corterr_B:.4}', linewidth=2)
#plt.plot(x_fit, langevin_fit, '-', label=f'Langevin fit: amp=({93}±{2})e-8, corte=({30}±{2})e-8, , mod={3.9}±{1.4}', linewidth=2)
#plt.plot(x_fit, brillouin_fit, '-', label=f'Brillouin fit: s={0.5}±{serr:.2}, amp=({91}±{2})e-8, corte=({30}±{2})e-8,mod={mod_L:.2}±{moderr_L:.2}', linewidth=2)
#plt.title('Ajustes de los modelos Langevin y Brillouin', fontsize = 22)
#plt.xlabel('H/T (Oe/K)', fontsize = 22)
#plt.ylabel('Momento magnetico del Aluminio (emu)', fontsize = 22)
#plt.legend()
#plt.grid()
#plt.savefig('langevin_fit.png')
#plt.show()
#plt.close()
#
#
#plt.figure(figsize=(10,8))
#plt.plot(x_, mu_Al - ajuste_Langevin(x_, moderr_L, amp_L, corte_L), '.', label='Residuos modelo de Langevin')
#plt.plot(x_, mu_Al - ajuste_Brillouin(x_, s, mod_B, amp_B, corte_B), '.', label='Residuos modelo de Brillouin s = 0.51')
#plt.axhline(0, color='k')
#plt.legend()
#plt.show()
dimasc = (H < 450) & (H > -450)
mualin = mu_Al[dimasc]
hlin = H[dimasc]

popt1, pcov1 = curve_fit(linear, hlin, mualin)
Xal, b = popt1
errXal, errb = np.sqrt(np.diag(pcov1))


plt.figure()
plt.plot(hlin, mualin, 'o', label='Datos del régimen lineal', markersize=5)
plt.plot(hlin, linear(hlin, Xal, b), '-', label=f'Ajuste lineal: Xal={Xal:.2e}±{errXal:.2e}', linewidth=2)
plt.title('Ajuste lineal del momento magnetico del Aluminio en campo debil', fontsize=16)
plt.xlabel('Campo  magnético (Oe)', fontsize=14)
plt.ylabel('Magnetización (emu)', fontsize=14)
plt.legend()
plt.grid()
plt.show()
m = 0.2035
errm = 0.0002
Xm = Xal / m
errXm = np.sqrt((errXal / Xal)**2 + (errm * Xal / m)**2)
print(f'Susceptibilidad magnetica del Aluminio: Xal = {Xal:.3e} ± {errXal:.5e}')
print(f'Susceptibilidad magnetica masica del Aluminio: Xm = {Xm:.1} ± {errXm:.11}')

v = 0.07744

plt.figure()
plt.plot(H, mu_Al/v, 'o', markersize=5)
plt.title('Magnetización por unuidad de volumen del Aluminio', fontsize=16)
plt.xlabel('Campo magnético (Oe)', fontsize=14)
plt.ylabel('Magnetización (emu/cm³)', fontsize=14)
plt.grid()
plt.show()