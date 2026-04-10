#!/usr/bin/env python3
"""
Script to plot polarization intensity P with LIC overlay showing polarization angles.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.colors import ListedColormap
import healpy as hp

# Configuration
fits_file = "512_60.0smoothed_quijotecombwei10_tqu_v1.5_noise_v1.0_-3.1_wmapplanck_qtmask_combine.fits"
LIC_transparency = 0.7  # Higher value = more visible LIC
foreground_cmap = 'jet'  # Colormap for P intensity
limits = [0, 0.7]  # Linear normalization limits in mK
# Smoothing parameters
fwhm_initial = 1.0  # Current FWHM in degrees
fwhm_final = 1.0    # Target FWHM for LIC in degrees
# Calculate effective smoothing FWHM
fwhm_smooth = np.sqrt(fwhm_final**2 - fwhm_initial**2)  # = sqrt(4 - 1) = sqrt(3) ≈ 1.73 deg

print("=" * 80)
print(f"Reading FITS file: {fits_file}")
print("=" * 80)

# Read Q and U maps (fields 1 and 2 based on the header)
# Field 0 = P, Field 1 = Q, Field 2 = U
P = hp.read_map(fits_file, field=0, verbose=False, dtype=np.float64)
Q = hp.read_map(fits_file, field=1, verbose=False, dtype=np.float64)
U = hp.read_map(fits_file, field=2, verbose=False, dtype=np.float64)

print(f"\nMap statistics:")
print(f"P - Min: {np.nanmin(P):.3f}, Max: {np.nanmax(P):.3f}, Mean: {np.nanmean(P):.3f}")
print(f"Q - Min: {np.nanmin(Q):.3f}, Max: {np.nanmax(Q):.3f}, Mean: {np.nanmean(Q):.3f}")
print(f"U - Min: {np.nanmin(U):.3f}, Max: {np.nanmax(U):.3f}, Mean: {np.nanmean(U):.3f}")

# Store original for masking
P_original = np.copy(P)
Q_original = np.copy(Q)
U_original = np.copy(U)

# Smooth Q and U to 2 degrees FWHM for smoother LIC
print(f"\nSmoothing Q and U from {fwhm_initial}° to {fwhm_final}° FWHM...")
print(f"Effective smoothing: {fwhm_smooth:.3f}°")
Q_smooth = hp.smoothing(Q, fwhm=np.deg2rad(fwhm_smooth))
U_smooth = hp.smoothing(U, fwhm=np.deg2rad(fwhm_smooth))
print("Smoothing done!")

# LIC does NOT like np.nans AT ALL!
P[np.isnan(P_original)] = 0
Q_smooth[np.isnan(Q_original)] = 0
U_smooth[np.isnan(U_original)] = 0

# Compute LIC with correct sign combination: -Q, -U
print("\nComputing LIC with -Q, -U...")
lic = hp.line_integral_convolution(-Q_smooth, -U_smooth)
print("LIC computation done!")

# Smooth the LIC slightly
print("Smoothing LIC...")
lic = hp.smoothing(lic, np.deg2rad(0.5))
print("Done!")

# Restore NaNs
P[np.isnan(P_original)] = np.nan
lic[np.isnan(P_original)] = np.nan

# Restore NaNs for P
P[np.isnan(P_original)] = np.nan

print(f"\nUsing P limits: [{limits[0]:.3f}, {limits[1]:.3f}] mK")

# Prepare LIC colormap with transparency
cmap_colors = plt.colormaps['binary'](np.linspace(0, 1, 256))
cmap_colors[..., 3] = LIC_transparency
lic_cmap = ListedColormap(cmap_colors)

# Create final plots
print("\nCreating plots...")

# Plot 1: Mollweide projection
plt.figure(10, figsize=(12, 7))
hp.mollview(P, fig=10, cmap=foreground_cmap,
            xsize=2400, norm='hist',
            title='Polarization Intensity with LIC overlay',
            cbar=False)
hp.mollview(lic, cmap=lic_cmap, cbar=False, reuse_axes=True,
            xsize=2400, title='')

plt.tight_layout()
plt.savefig('P_with_LIC_mollweide.pdf', bbox_inches='tight', pad_inches=0)
plt.savefig('P_with_LIC_mollweide.png', dpi=300, bbox_inches='tight', pad_inches=0)
print("Saved: P_with_LIC_mollweide.pdf and .png")

# Plot 2: Orthographic projection
plt.figure(11, figsize=(12, 7))
hp.orthview(P, rot=[0, 45], fig=11, cmap=foreground_cmap,
            min=limits[0], max=limits[1], xsize=2400,
            title='Polarization Intensity with LIC overlay (Orthographic)',
            unit='mK CMB')
hp.orthview(lic, rot=[0, 45], cmap=lic_cmap, cbar=False,
            reuse_axes=True, xsize=2400, title='')

plt.tight_layout()
plt.savefig('P_with_LIC_orthographic.pdf', bbox_inches='tight', pad_inches=0)
plt.savefig('P_with_LIC_orthographic.png', dpi=300, bbox_inches='tight', pad_inches=0)
print("Saved: P_with_LIC_orthographic.pdf and .png")

# Plot 3: Just the LIC for reference
plt.figure(12, figsize=(11, 6))
hp.mollview(lic, cmap='binary', title='LIC (smoothed, -Q, -U)',
            xsize=2400, cbar=True)
plt.tight_layout()
plt.savefig('LIC_only.png', dpi=300, bbox_inches='tight', pad_inches=0)
print("Saved: LIC_only.png")

print("\n" + "=" * 80)
print("Done! All plots have been saved.")
print("=" * 80)

plt.show()
