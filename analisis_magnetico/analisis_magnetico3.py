import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import binned_statistic

dH = 60.0

with open('magnetismo.txt', 'r') as file:
    content = file.read()
    lineas = content.split('\n')

pos = lineas.index('[Data]')
headers = lineas[pos+1].split(',')
dic = {x: [] for x in headers}

for i in range(pos+2, len(lineas)):
    if lineas[i].strip() == '':
        continue
    data = lineas[i].split(',')
    for j in range(len(headers)):
        dic[headers[j]].append(data[j])

mu = np.array(dic['Moment (emu)'], dtype=float)
H  = np.array(dic['Magnetic Field (Oe)'], dtype=float)
T  = np.array(dic['Temperature (K)'], dtype=float)

x = H / T

def linear(x, m, b):
    return m*x + b

mask_pos = H > 0
mask_neg = H < 0

popt_p, _ = curve_fit(linear, H[mask_pos], mu[mask_pos])
popt_n, _ = curve_fit(linear, H[mask_neg], mu[mask_neg])

Xp, bp = popt_p
Xn, bn = popt_n

mu_Al = np.zeros_like(mu)
mu_Al[mask_pos] = mu[mask_pos] - Xp*H[mask_pos]
mu_Al[mask_neg] = mu[mask_neg] - Xn*H[mask_neg]

bins_H = np.arange(H.min(), H.max() + dH, dH)

mu_mean, _, _ = binned_statistic(H, mu_Al, statistic='mean', bins=bins_H)
mu_std,  _, _ = binned_statistic(H, mu_Al, statistic='std',  bins=bins_H)
x_mean,  _, _ = binned_statistic(H, x,      statistic='mean', bins=bins_H)
counts,  _, _ = binned_statistic(H, mu_Al,  statistic='count', bins=bins_H)

mask = (~np.isnan(mu_mean)) & (counts > 2)

x_b   = x_mean[mask]
mu_b  = mu_mean[mask]
err_b = mu_std[mask] / np.sqrt(counts[mask])

def coth(x):
    return np.cosh(x)/np.sinh(x)

def ajuste_Langevin(x, a, b, c):
    return b*(coth(a*x) - 1/(a*x)) + c

def ajuste_Brillouin(x, s, a, b, c):
    t1 = (2*s + 1)/(2*s)
    t2 = 1/(2*s)
    return b*(t1*coth(t1*a*x) - t2*coth(a*x/(2*s))) + c

p0_L = [1.0, 1e-6, 0.0]
popt_L, pcov_L = curve_fit(
    ajuste_Langevin,
    x_b, mu_b,
    p0=p0_L,
    sigma=err_b,
    absolute_sigma=True,
    maxfev=20000
)

p0_B = [0.5, 1.0, 1e-6, 0.0]
popt_B, pcov_B = curve_fit(
    ajuste_Brillouin,
    x_b, mu_b,
    p0=p0_B,
    sigma=err_b,
    absolute_sigma=True,
    maxfev=20000
)

def chi2_red(y, yfit, sigma, npar):
    return np.sum(((y - yfit)/sigma)**2)/(len(y) - npar)

chi_L = chi2_red(mu_b, ajuste_Langevin(x_b, *popt_L), err_b, len(popt_L))
chi_B = chi2_red(mu_b, ajuste_Brillouin(x_b, *popt_B), err_b, len(popt_B))

x_fit = np.linspace(x_b.min(), x_b.max(), 500)

plt.figure(figsize=(10,8))
plt.errorbar(x_b, mu_b, yerr=err_b, fmt='o', ms=4, label='Datos binneados')
plt.plot(x_fit, ajuste_Langevin(x_fit, *popt_L), '-', lw=2,
         label=f'Langevin χ²ᵣ={chi_L:.2f}')
plt.plot(x_fit, ajuste_Brillouin(x_fit, *popt_B), '--', lw=2,
         label=f'Brillouin χ²ᵣ={chi_B:.2f}')
plt.title('Ajustes de los modelos Langevin y Brillouin', fontsize=22)
plt.xlabel('H/T (Oe/K)', fontsize=22)
plt.ylabel('Momento magnético del Aluminio (emu)', fontsize=22)
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()

plt.figure(figsize=(10,8))
plt.plot(x_b, mu_b - ajuste_Langevin(x_b, *popt_L), '.', label='Residuos Langevin')
plt.plot(x_b, mu_b - ajuste_Brillouin(x_b, *popt_B), '.', label='Residuos Brillouin')
plt.axhline(0, color='k')
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()

