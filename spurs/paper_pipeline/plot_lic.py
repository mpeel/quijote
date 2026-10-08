''' plot_lic.py:
Bottom panel of Fig. 3: P of the full weighted map with a line integral convolution
(LIC) overlay showing the polarisation angles.

Example usage:
python plot_lic.py

Version 1.0 [Oct 2026]
Roke Cepeda-Arroita
roke.cepeda@iac.es
'''

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import healpy as hp
from scipy import stats
import os
from cmcrameri import cm as cmc

from config import PRODUCTS, FIGURES

# Configuration
fits_file = f'{PRODUCTS}/combined/full_60arcmin_n512.fits'
OUTDIR = f'{FIGURES}/fig3/'
os.makedirs(OUTDIR, exist_ok=True)
LIC_transparency = 0.35  # Higher value = more visible LIC (reduced from 0.7 for less washout)
foreground_cmap = cmc.roma_r  # Colormap for P intensity — matches the loops figure
limits = [0, 0.7]  # Linear normalization limits in mK
# Smoothing parameters
fwhm_initial = 1.0  # Current FWHM in degrees
fwhm_final = 1.0    # Target FWHM for LIC in degrees (no additional smoothing)
# Calculate effective smoothing FWHM
fwhm_smooth = np.sqrt(max(0, fwhm_final**2 - fwhm_initial**2))

print("=" * 80)
print(f"Reading FITS file: {fits_file}")
print("=" * 80)

# Read Q and U maps (fields 1 and 2 based on the header)
# Field 0 = P, Field 1 = Q, Field 2 = U
P = hp.read_map(fits_file, field=2, dtype=np.float64)   # columns Q, U, P
Q = hp.read_map(fits_file, field=0, dtype=np.float64)
U = hp.read_map(fits_file, field=1, dtype=np.float64)

print(f"\nMap statistics:")
print(f"P - Min: {np.nanmin(P):.3f}, Max: {np.nanmax(P):.3f}, Mean: {np.nanmean(P):.3f}")
print(f"Q - Min: {np.nanmin(Q):.3f}, Max: {np.nanmax(Q):.3f}, Mean: {np.nanmean(Q):.3f}")
print(f"U - Min: {np.nanmin(U):.3f}, Max: {np.nanmax(U):.3f}, Mean: {np.nanmean(U):.3f}")

# Store original for masking
P_original = np.copy(P)
Q_original = np.copy(Q)
U_original = np.copy(U)

# Smooth Q and U if needed
if fwhm_smooth > 0:
    print(f"\nSmoothing Q and U from {fwhm_initial}° to {fwhm_final}° FWHM...")
    print(f"Effective smoothing: {fwhm_smooth:.3f}°")
    Q_smooth = hp.smoothing(Q, fwhm=np.deg2rad(fwhm_smooth))
    U_smooth = hp.smoothing(U, fwhm=np.deg2rad(fwhm_smooth))
    print("Smoothing done!")
else:
    print(f"\nNo additional smoothing of Q and U (already at {fwhm_initial}° FWHM)")
    Q_smooth = np.copy(Q)
    U_smooth = np.copy(U)

# LIC does NOT like np.nans AT ALL!
P[np.isnan(P_original)] = 0
Q_smooth[np.isnan(Q_original)] = 0
U_smooth[np.isnan(U_original)] = 0

# Compute LIC with correct sign combination: -Q, -U
print("\nComputing LIC with -Q, -U...")
lic = hp.line_integral_convolution(-Q_smooth, -U_smooth)
print("LIC computation done!")

# Apply minimal smoothing to LIC for sharper result
print("Applying minimal smoothing to LIC...")
lic = hp.smoothing(lic, np.deg2rad(0.3))  # Reduced from 1.0 to 0.3 degrees
print("Done!")

# Enhance LIC contrast using histogram equalization
print("Enhancing LIC contrast...")
lic_valid = lic[~np.isnan(P_original)]
lic_normalized = stats.rankdata(lic_valid) / len(lic_valid)
lic_enhanced = np.empty_like(lic)
lic_enhanced[~np.isnan(P_original)] = lic_normalized
lic_enhanced[np.isnan(P_original)] = np.nan
lic = lic_enhanced
print("Done!")

# Restore NaNs
P[np.isnan(P_original)] = np.nan
lic[np.isnan(P_original)] = np.nan

print(f"\nUsing P limits: [{limits[0]:.3f}, {limits[1]:.3f}] mK")
print(f"LIC transparency: {LIC_transparency}")

# Prepare LIC colormap with transparency
# Use 'gray' for sharper contrast, and limit range for better definition
cmap_colors = plt.colormaps['gray'](np.linspace(0.2, 0.8, 256))  # Constrained range for sharper contrast
cmap_colors[..., 3] = LIC_transparency
lic_cmap = ListedColormap(cmap_colors)

# Create final plots
print("\nCreating plots...")

# Plot 1: Mollweide projection
plt.figure(10, figsize=(12, 7))
hp.mollview(P, fig=10, cmap=foreground_cmap,
            xsize=2400, norm='hist',
            title='Polarization Intensity with Sharp LIC overlay',
            cbar=False)
hp.mollview(lic, cmap=lic_cmap, cbar=False, reuse_axes=True,
            xsize=2400, title='')

plt.tight_layout()
plt.savefig(OUTDIR + 'P_with_LIC_sharp_mollweide.pdf', bbox_inches='tight', pad_inches=0)
plt.savefig(OUTDIR + 'P_with_LIC_sharp_mollweide.png', dpi=300, bbox_inches='tight', pad_inches=0)
print("Saved: P_with_LIC_sharp_mollweide.pdf and .png")

# Plot 2: Orthographic projection
plt.figure(11, figsize=(12, 7))
hp.orthview(P, rot=[0, 45], fig=11, cmap=foreground_cmap,
            min=limits[0], max=limits[1], xsize=2400,
            title='Polarization Intensity with Sharp LIC overlay (Orthographic)',
            unit='mK CMB')
hp.orthview(lic, rot=[0, 45], cmap=lic_cmap, cbar=False,
            reuse_axes=True, xsize=2400, title='')

plt.tight_layout()
plt.savefig(OUTDIR + 'P_with_LIC_sharp_orthographic.pdf', bbox_inches='tight', pad_inches=0)
plt.savefig(OUTDIR + 'P_with_LIC_sharp_orthographic.png', dpi=300, bbox_inches='tight', pad_inches=0)
print("Saved: P_with_LIC_sharp_orthographic.pdf and .png")

# Plot 3: Just the LIC for reference
plt.figure(12, figsize=(11, 6))
hp.mollview(lic, cmap='gray', title='LIC (enhanced, minimal smoothing, -Q, -U)',
            xsize=2400, cbar=True)
plt.tight_layout()
plt.savefig(OUTDIR + 'LIC_only_sharp.png', dpi=300, bbox_inches='tight', pad_inches=0)
print("Saved: LIC_only_sharp.png")

print("\n" + "=" * 80)
print("Done! All plots have been saved.")
print("=" * 80)
print("\nKey improvements:")
print(f"  - LIC transparency reduced to {LIC_transparency} (was 0.7)")
print(f"  - LIC smoothing reduced to 0.3° (was 1.0°)")
print("  - Added histogram equalization for enhanced contrast")
print("  - Using constrained 'gray' colormap range (0.2-0.8) for sharper features")
print("=" * 80)
