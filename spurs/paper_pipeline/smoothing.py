''' smoothing.py:
Smooths I, Q, U maps to a common Gaussian resolution using the full beam
transfer function of each map (not a Gaussian approximation to the beam).

The smoothing is done in harmonic space, with Q and U as a spin-2 field:

    a_lm(out) = a_lm(in) * G_l / (b_l * p_l)

where G_l is the target Gaussian, b_l the beam of the input map and p_l the
HEALPix pixel window of the input Nside. The map stays at its native Nside.
This is a linear operation, so the same call is used for the data and for
every noise realisation.

Example usage:
I, Q, U = smooth_iqu(I, Q, U, 'K', 60.0)

Version 1.0 [Oct 2026]
Roke Cepeda-Arroita
roke.cepeda@iac.es
'''

import os
import numpy as np
import healpy as hp
from functools import lru_cache
from astropy.io import fits

from config import BEAM_DIR


# Beam Files

# WMAP: official 9-year b_l of each differencing assembly
WMAP_DA = {'K': 'K1', 'Ka': 'Ka1', 'Q1': 'Q1', 'Q2': 'Q2'}

# Planck LFI: NPIPE T, E and B beams (the ones used in the Cosmoglobe processing)
LFI_FREQ = {'030': '30', '044': '44'}

# QUIJOTE MFI DR1: b_l of each horn and frequency from the RIMO
MFI_HORN = {'MFI11': '311', 'MFI13': '313', 'MFI17H2': '217', 'MFI17H4': '417', 'MFI19H2': '219', 'MFI19H4': '419'}

# Nominal FWHM in arcmin, only used for the Gaussian tail of the kernel (see transfer_functions)
NOMINAL_FWHM = {'K': 52.8, 'Ka': 39.6, 'Q1': 30.6, 'Q2': 30.6,
                '030': 32.29, '044': 27.94,
                'MFI11': 55.38, 'MFI13': 55.84, 'MFI17H2': 38.3, 'MFI17H4': 38.8, 'MFI19H2': 38.0, 'MFI19H4': 38.8}

# Multipoles where the target Gaussian is below this value are set to zero
GAUSS_FLOOR = 1e-10


def load_beam(band):
    ''' Returns the beam transfer functions (bT, bE, bB) of a band '''

    if band in WMAP_DA:
        d = np.loadtxt(os.path.join(BEAM_DIR, f'wmap_ampl_bl_{WMAP_DA[band]}_9yr_v5p1.txt'))
        bl = np.zeros(int(d[-1, 0])+1)
        bl[d[:, 0].astype(int)] = d[:, 1]
        return bl, bl, bl

    elif band in LFI_FREQ:
        f = LFI_FREQ[band]
        d = fits.getdata(os.path.join(BEAM_DIR, f'Bl_TEB_npipe6v19_{f}GHzx{f}GHz.fits'), 1)
        return d['T'].astype(float), d['E'].astype(float), d['B'].astype(float)

    elif band in MFI_HORN:
        d = fits.getdata(os.path.join(BEAM_DIR, 'quijote_mfi', 'rimo_quijote_mfi_beamtf_dr1.fits'), 1)
        bl = d[f'Bl_{MFI_HORN[band]}'][0].astype(float)
        return bl, bl, bl

    else:
        raise ValueError(f'Unknown band {band}')


def gaussian_beam(fwhm_arcmin, lmax, pol=False):
    ''' Gaussian beam transfer function, (T, E, B) if pol=True '''

    g = hp.gauss_beam(np.deg2rad(fwhm_arcmin/60.), lmax=lmax, pol=True)
    if pol:
        return g[:, 0], g[:, 1], g[:, 2]
    return g[:, 0]


def default_lmax(band, nside_in, fwhm_final):
    ''' Highest multipole needed: where the target Gaussian drops below GAUSS_FLOOR,
        limited by 3*Nside-1 and by the length of the tabulated beam '''

    l = np.arange(4*nside_in)
    g = gaussian_beam(fwhm_final, lmax=len(l)-1)
    if np.any(g < GAUSS_FLOOR):
        l_floor = int(l[np.argmax(g < GAUSS_FLOOR)])
    else:
        l_floor = len(l)-1
    return min(l_floor, 3*nside_in-1, len(load_beam(band)[0])-1)


@lru_cache(maxsize=None)
def transfer_functions(band, nside_in, fwhm_final=60.0, lmax=None, deconvolve_pixwin=True):
    ''' Smoothing kernels (fT, fE, fB) from the native beam of a band to a Gaussian
        of FWHM fwhm_final (arcmin). Cached, so they are only computed once per band '''

    if lmax is None:
        lmax = default_lmax(band, nside_in, fwhm_final)

    # Beams (set to zero beyond the tabulated multipoles), target Gaussian and pixel window
    bT, bE, bB = (np.pad(b, (0, max(0, lmax+1-len(b))))[:lmax+1] for b in load_beam(band))
    gT, gE, gB = gaussian_beam(fwhm_final, lmax=lmax, pol=True)

    if deconvolve_pixwin:
        pT, pP = hp.pixwin(nside_in, pol=True, lmax=lmax)
    else:
        pT = pP = np.ones(lmax+1)

    def ratio(g, b, p):
        f = np.zeros(lmax+1)
        ok = (b*p) > 0
        f[ok] = g[ok]/(b[ok]*p[ok])
        return f

    # Gaussian tail: where the native beam falls faster than the target Gaussian
    # (QUIJOTE horn 3 above l~350), G/b starts to grow with l and would amplify the
    # noise. Above the first l where G/b stops decreasing (searched only where the
    # target is below 0.5), the kernel continues as the Gaussian of FWHM
    # sqrt(fwhm_final^2 - fwhm_nominal^2), matched at that l
    l_search = int(np.argmax(gT < 0.5))
    k = gaussian_beam(np.sqrt(fwhm_final**2 - NOMINAL_FWHM[band]**2), lmax=lmax)

    def gaussian_tail(f):
        rising = np.where(np.diff(f[l_search:]) >= 0)[0]
        if len(rising) == 0:
            return f, None
        l_cross = l_search + int(rising[0])
        f = f.copy()
        f[l_cross:] = f[l_cross] * k[l_cross:] / k[l_cross]
        return f, l_cross

    fT, l_cross_T = gaussian_tail(ratio(gT, bT, pT))
    fE, l_cross_E = gaussian_tail(ratio(gE, bE, pP))
    fB, _ = gaussian_tail(ratio(gB, bB, pP))
    if l_cross_T is not None or l_cross_E is not None:
        print(f'{band}: Gaussian tail used above l = {l_cross_T} (T), {l_cross_E} (E/B)')

    # No monopole or dipole in polarisation
    fE[:2] = 0.
    fB[:2] = 0.

    return fT, fE, fB, lmax


def smooth_iqu(I, Q, U, band, fwhm_final=60.0, lmax=None, iter=0, deconvolve_pixwin=True):
    ''' Smooths I, Q, U (RING) of a band to a Gaussian of FWHM fwhm_final (arcmin).
        Pixels with no data are filled with the median before smoothing and set
        back to NaN afterwards '''

    nside = hp.get_nside(I)
    fT, fE, fB, lmax = transfer_functions(band, nside, fwhm_final, lmax, deconvolve_pixwin)

    # Fill the pixels with no data
    maps = np.array([I, Q, U], dtype=float)
    maps[maps == hp.UNSEEN] = np.nan
    bad = np.isnan(maps)
    if np.any(bad):
        med = np.nanmedian(maps, axis=1)
        for i in range(3):
            maps[i, bad[i]] = med[i]

    # Smooth in harmonic space
    alm = hp.map2alm(maps, lmax=lmax, iter=iter, pol=True)
    alm = [hp.almxfl(a, f) for a, f in zip(alm, (fT, fE, fB))]
    out = hp.alm2map(alm, nside, lmax=lmax, pol=True)

    out[bad] = np.nan
    return out[0], out[1], out[2]
