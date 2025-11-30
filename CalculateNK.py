




"""
fit_nk.py
Fit n(λ) and k(λ) from specular reflectance measurements at 30°, 45°, 60°.
Input CSV columns: wavelength_nm, R30, R45, R60
Outputs a CSV with wavelength_nm, n, k, and plots.

Dependencies: numpy, scipy, pandas, matplotlib
"""

import numpy as np
import pandas as pd
from scipy.optimize import least_squares
import matplotlib
matplotlib.use('TkAgg')  # o prueba 'Agg' si no vas a mostrar ventanas
import matplotlib.pyplot as plt

# --- Fresnel helper functions ---
def fresnel_RsRp(N, theta_inc_rad):
    """
    Compute Rs and Rp for a bulk material with complex refractive index N = n + 1j*k,
    at incidence angle theta_inc_rad (radians) in air (n0=1).
    Returns (Rs, Rp) as floats.
    """
    n0 = 1.0
    # Snell: kz in medium (complex)
    sin_t = n0 * np.sin(theta_inc_rad) / N  # complex
    # sqrt argument for cos_t: sqrt(1 - sin^2)
    cos_t = np.sqrt(1.0 - sin_t**2)
    # cos incidence
    cos_i = np.cos(theta_inc_rad)
    # Fresnel coefficients (complex)
    # i 
    rs = (n0 * cos_i - N * cos_t) / (n0 * cos_i + N * cos_t)
    rp = (N * cos_i - n0 * cos_t) / (N * cos_i + n0 * cos_t)
    Rs = np.abs(rs)**2
    Rp = np.abs(rp)**2
    return Rs.real, Rp.real

# --- Objective for least-squares (unpolarized average) ---
def residuals_nk(x, R_meas, angles_rad):     
    """
    x = [n, k]
    R_meas: array-like measured reflectances at the supplied angles (same order)
    angles_rad: array-like angles in radians
    returns residual vector (len = len(angles_rad))
    """
    n, k = x
    if n < 0 or k < 0:
        # return large residuals to keep solver in physical domain
        return 1e3 * np.ones_like(R_meas)
    N = n + 1j * k
    res = []
    for idx, th in enumerate(angles_rad):
        Rs, Rp = fresnel_RsRp(N, th)
        Rcalc = 0.5 * (Rs + Rp)   # unpolarized average
        res.append(Rcalc - R_meas[idx])
    return np.asarray(res)

# --- Main fit routine ---
def fit_nk_per_wavelength(wvl, Rvals, angles_deg, n0_guess=1.0, k0_guess=1.0):
    """
    Fit n,k for a single wavelength.
    Rvals: measured reflectances list-like, same order as angles_deg
    angles_deg: list like [30,45,60]
    Returns (n_fit, k_fit) or (np.nan, np.nan) if fit fails.
    """
    angles_rad = np.deg2rad(angles_deg)
    # initial guess
    x0 = np.array([n0_guess, k0_guess])
    # bounds: n in [0, 10], k in [0, 50] (adjust if needed)
    lb = [0.0, 0.0]   ############## CORREGIR PARA QUE DE n>1
    ub = [20.0, 200.0]
    try:
        res = least_squares(residuals_nk, x0, args=(Rvals, angles_rad), bounds=(lb, ub),
                            ftol=1e-10, xtol=1e-10, gtol=1e-10, max_nfev=2000)
        n_fit, k_fit = res.x
        return float(n_fit), float(k_fit), res.cost
    except Exception as e:
        print(f"Fit failed at {wvl} nm: {e}")
        return np.nan, np.nan, np.nan

# --- Example usage: load data, loop wavelengths ---
def main(input_csv='filt_reflectance_data.csv', output_csv='filt_nk_fitted.csv'):
    df = pd.read_csv(input_csv)
    required = ['wavelength_nm', 'R30', 'R45', 'R60']
    for col in required:
        if col not in df.columns:
            raise ValueError(f"Input CSV must contain column '{col}'")
    wavelengths = df['wavelength_nm'].values
    R30 = df['R30'].values
    R45 = df['R45'].values
    R60 = df['R60'].values
    angles = [30.0, 45.0, 60.0]

    n_list, k_list, cost_list = [], [], []
    # Good initial guess strategy: use previous wavelength result for continuity
    n_guess, k_guess = 1.0, 5.0

    for i, wl in enumerate(wavelengths):
        Rvals = [R30[i], R45[i], R60[i]]
        # sanity clamp
        Rvals = np.clip(Rvals, 1e-8, 0.999999)
        n_fit, k_fit, cost = fit_nk_per_wavelength(wl, Rvals, angles, n_guess, k_guess)
        # fallback if failed: try different initial guess
        if np.isnan(n_fit):
            n_fit, k_fit, cost = fit_nk_per_wavelength(wl, Rvals, angles, 1.0, 5.0)
        n_list.append(n_fit)
        k_list.append(k_fit)
        cost_list.append(cost)
        # update guess for next wavelength (smoothness)
        if not np.isnan(n_fit):
            n_guess, k_guess = n_fit, k_fit

    out = pd.DataFrame({
        'wavelength_nm': wavelengths,
        'n': n_list,
        'k': k_list,
        'cost': cost_list,
        'R30': R30, 'R45': R45, 'R60': R60
    })
    out.to_csv(output_csv, index=False)
    print(f"Saved fitted n,k to {output_csv}")

    # quick plots
    fig, axes = plt.subplots(2,1, figsize=(7,6), sharex=True)
    axes[0].plot(wavelengths, out['n'], marker='o', label='n')
    axes[0].set_ylabel('n')
    axes[0].grid(True)
    axes[1].plot(wavelengths, out['k'], marker='o', label='k')
    axes[1].set_ylabel('k')
    axes[1].set_xlabel('Wavelength (nm)')
    axes[1].grid(True)
    plt.suptitle('Fitted n and k')
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    main("C:/Users/Asus/caracterizacion_pelicula_delgada/filt_reflectance_data.csv")