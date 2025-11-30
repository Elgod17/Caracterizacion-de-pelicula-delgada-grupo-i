# Copilot Instructions: Optical Thin Film Analysis

## Project Overview

This project analyzes optical spectra from Thorlabs FTS (Fourier Transform Spectrometer) measurements of thin film samples (30nm thickness designation). The workflow processes spectral intensity data across wavelengths to characterize optical properties (reflectance, transmission, etc.).

**Key data sources**: CSV files with Thorlabs FTS metadata headers followed by wavelength-intensity data pairs.

## Data Processing Patterns

### CSV Format & Parsing
- **Format**: Thorlabs FTS header (lines 1-34, metadata with `#Key;Value` syntax) + data rows (columns: wavelength in nm_air, intensity)
- **Parsing pattern**: Skip header with `iloc[34:, 0:2]` to extract wavelength column (0) and intensity column (1)
- **Files**: 
  - `DR-G1-30G.csv`, `dr30.csv` - Direct reflectance measurements
  - `R-30G.csv`, `r30.csv` - Additional reflectance datasets

### NumPy/Pandas Workflow
- Convert columns to NumPy arrays for curve fitting: `.to_numpy()`
- Use `iloc` for header skipping (row ranges) and column selection (0-indexed)
- `scipy.optimize.curve_fit()` for spectral feature modeling (e.g., resonances, absorption)

## Common Tasks

### Plotting & Visualization
- **Backend**: Use `matplotlib.use('TkAgg')` explicitly before importing pyplot (set before interactive operations)
- **Pattern**: Plot intensity vs. wavelength; use `plt.show()` for GUI display on Windows
- **Example**: `plt.plot(intensity_array, wavelength_array)` 

### Adding New Analyses
1. Load CSV → skip header → extract columns
2. Optional: Apply smoothing/normalization to intensity data
3. Fit physical models (Fresnel equations, optical constants, resonance models)
4. Visualize results with wavelength on x-axis, intensity/reflectance on y-axis

## File Organization

```
optico_30.py          # Main analysis script
DR-G1-30G.csv         # Measurement datasets (naming: MEASUREMENT_ID-SAMPLE_ID.csv)
dr30.csv              # 
R-30G.csv
r30.csv
```

## Development Notes

- **Absolute paths required**: CSV paths are hardcoded with full Windows paths (`C:/Users/Asus/...`); update path references when relocating files or running on different systems
- **Header structure consistency**: All CSVs follow Thorlabs FTS format; verify header line count (34 rows) before extracting data
- **Numerical focus**: Emphasis on NumPy arrays and scipy optimization for spectral fitting
