import numpy as np
import matplotlib
matplotlib.use('TkAgg')  # o prueba 'Agg' si no vas a mostrar ventanas
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
from scipy.constants import Boltzmann

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

mu = np.array(dic['Moment (emu)'])
H = np.array(dic['Magnetic Field (Oe)'])
pos_empty = np.where(mu == '')[0]
mu = np.delete(mu, pos_empty).astype(float)
H = np.delete(H, pos_empty).astype(float)
temp = np.delete(np.array(dic['Temperature (K)']), pos_empty).astype(float)
H_2 = H[H<=0]
mu_2 = mu[H<=0]
plt.figure(figsize = (10,8))
plt.plot(H, mu, marker='o', linestyle='-', color='b', markersize=5)
plt.title('Magnetic Moment vs Magnetic Field')
plt.xlabel('Magnetic Field (Oe)')
plt.ylabel('Magnetic Moment (emu)')
plt.grid()
plt.savefig('moment_vs_field.png')
plt.show()
plt.close()

### 2. Ajuste de minimos cuadrados 

def linear(x, m, b):
    return m*x + b

popt, pcov = curve_fit(linear, H_2, mu_2)
X, b_fit = popt
delta_m, delta_b = np.sqrt(np.diag(pcov))
H_fit = np.linspace(min(H_2), max(H_2), 100)
mu_fit = linear(H_fit, X, b_fit)
plt.figure(figsize=(10,8))
plt.plot(H_2, mu_2, 'o', label='Datos experimentales')
plt.plot(H_fit, mu_fit, '-', label=f'y = {X:.2}x+{b_fit:.2}')
plt.title('Ajuste Lineal de Momento Magnético vs Campo Magnético')
plt.xlabel('Magnetic Field (Oe)')
plt.ylabel('Magnetic Moment (emu)')
plt.legend()
plt.grid()
plt.savefig('fit_moment_vs_field.png')
plt.show()
plt.close()

print(f'Pendiente (m): {X:.2} ± {delta_m:.2}; Intersección (b): {b_fit:.2} ± {delta_b:.2}')

### Susceptibilidad magnética del aluminio 

dic['Al_moment (emu)'] = mu - X*H

### Comportamiento real del aluminio vs campo magnético 

mu_Al = np.array(dic['Al_moment (emu)'])

plt.figure(figsize = (10,8))
plt.scatter(H, mu_Al, marker='o', color='r', s = 25)
plt.plot(H, mu_Al, linestyle='--', color='black', alpha = 0.7)
plt.title('Magnetic Moment vs Magnetic Field')
plt.xlabel('Magnetic Field (Oe)')
plt.ylabel('Magnetic Moment Al (emu)')
plt.grid()
plt.savefig('moment_vs_field_Al.png')
plt.show()
plt.close()

### Comparación: Modelo Langevan o Brillouin

def coth(x):
    val = np.sinh(x)
    for j in val: 
        if j == 0:
            return ValueError("coth undefined")
    else: 
        return np.cosh(x)/np.sinh(x)

def Brillouin(x,s):
    T1 = 1 + 1 / (2*s)
    T2 = coth(T1*x)
    T3 = -1 /(2*s) * coth(x / (2*s))
    return T1*T2+T3

def ajuste_Brillouin(x,a,b):
    return b*Brillouin(x,a)

def ajuste_Langevin(x,a,b):
    return b*(coth(a*x) - 1/(a*x))

def Langevin(x):
    return coth(x) - 1/x 



s = np.arange(0.5,5.5, 0.5)
norm_Al = mu_Al/max(abs(mu_Al))
Brillouin_values = {str(x):[] for x in s}
e_Brillouin = {str(x):0 for x in s}
H_si = H * 1e3 / (4 * np.pi)  # Convertir Oe a A/m
x = (H_si / temp)*0.001  # x = μB*H/(kB*T), con μB en J/T
Langevin_values = Langevin(x)
for i in range(len(s)):
    Brillouin_values[str(s[i])] = Brillouin(x, s[i])
    e_Brillouin[str(s[i])] = np.sum((norm_Al - Brillouin_values[str(s[i])])**2)

e_Langevin = np.sum((norm_Al - Langevin_values)**2)

plt.figure(figsize=(10,8))
plt.plot(H, norm_Al, 'o', label='Datos experimentales', markersize=5)
plt.plot(H, Langevin_values, '-', label='Langevin', linewidth=2)
for i in range(len(s)):
    plt.plot(H, Brillouin_values[str(s[i])], '--', label=f'Brillouin s={s[i]}', linewidth=2)
plt.title('Comparación Modelo Langevin y Brillouin')
plt.xlabel('Magnetic Field (Oe)')
plt.ylabel('Normalized Magnetic Moment Al')
plt.legend()
plt.grid()
plt.savefig('model_comparison.png')
plt.show()
plt.close()
print(f'Error cuadrático Langevin: {e_Langevin:.4f}')
for i in range(len(s)):
    print(f'Error cuadrático Brillouin s={s[i]}: {e_Brillouin[str(s[i])]:.4f}')

### Ajuste Langevin y Brillouin

popt_L, pcov_L = curve_fit(ajuste_Langevin, x, norm_Al)
a_L, b_L = popt_L
delta_a_L, delta_b_L = np.sqrt(np.diag(pcov_L))
x_fit = np.linspace(min(x), max(x), 100)
langevin_fit = ajuste_Langevin(x_fit, a_L, b_L)
plt.figure(figsize=(10,8))
plt.plot(x, norm_Al, 'o', label='Datos experimentales', markersize=5)
plt.plot(x_fit, langevin_fit, '-', label=f'Langevin fit: a={a_L:.2}±{delta_a_L:.2}, b={b_L:.2}±{delta_b_L:.2}', linewidth=2)
plt.title('Ajuste Modelo Langevin')
plt.xlabel('x = μB*H/(kB*T)')
plt.ylabel('Normalized Magnetic Moment Al')
plt.legend()
plt.grid()
plt.savefig('langevin_fit.png')
plt.show()
plt.close()
popt_B, pcov_B = curve_fit(ajuste_Brillouin, x, norm_Al)
a_B, b_B = popt_B
delta_a_B, delta_b_B = np.sqrt(np.diag(pcov_B))
brillouin_fit = ajuste_Brillouin(x_fit, a_B, b_B)
plt.figure(figsize=(10,8))
plt.plot(x, norm_Al, 'o', label='Datos experimentales', markersize=5)
plt.plot(x_fit, brillouin_fit, '-', label=f'Brillouin fit: a={a_B:.2}±{delta_a_B:.2}, b={b_B:.2}±{delta_b_B:.2}', linewidth=2)
plt.title('Ajuste Modelo Brillouin')
plt.xlabel('x = μB*H/(kB*T)')
plt.ylabel('Normalized Magnetic Moment Al')
plt.legend()
plt.grid()
plt.savefig('brillouin_fit.png')
plt.show()
plt.close()
print(f'Ajuste Langevin: a={a_L:.2} ± {delta_a_L:.2}, b={b_L:.2} ± {delta_b_L:.2}')
print(f'Ajuste Brillouin: a={a_B:.2} ± {delta_a_B:.2}, b={b_B:.2} ± {delta_b_B:.2}')

