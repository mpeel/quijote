''' noise_mc.py:
Smoothed maps and Monte Carlo noise for one input map.

The map is smoothed to 60 arcmin at its native Nside (smoothing.py) and then
degraded to each Nside in NSIDES. The noise is propagated with NSIM white-noise
realisations drawn from the per-pixel noise of the map (independent in I, and
correlated in Q and U through the 2x2 QU covariance), each one smoothed and
degraded exactly like the data. The variances and the QU covariance are then
computed from the realisations at each Nside.

Noise inputs:
    Cosmoglobe DR1: I_RMS, Q_RMS, U_RMS and QU_RMS (the QU covariance)
    QUIJOTE MFI DR1: 1/WEI and COV_QU (the DR1 weights already include the
                     noise rescaling factors of Table 12 of Rubino-Martin et al. 2023)

The Planck 2018 solar dipole is subtracted from the Cosmoglobe intensity maps,
and the DR1 satellite-band mask is applied to the QUIJOTE maps (an output pixel
is kept only if all of its native pixels are usable).

Outputs, in PRODUCTS/<name>/:
    <name>_60arcmin_n<nside>.fits     I, Q, U, I_ERR, Q_ERR, U_ERR, QU_COV
    <name>_60arcmin_n64_noise_sims.npy    the Nside-64 realisations (NSIM, 3, 49152)

Example usage:
python noise_mc.py K
python noise_mc.py MFI11 --nsim 1000

Version 1.0 [Oct 2026]
Roke Cepeda-Arroita
roke.cepeda@iac.es
'''

import os
import time
import argparse
import numpy as np
import healpy as hp
from astropy.io import fits

from config import CG_DIR, MFI_DIR, MFI_MASK, PRODUCTS, FWHM, NSIDES, NSIM
from smoothing import smooth_iqu


# Input Maps

# label: (file, kind, unit, output name, frequency in GHz)
DATASETS = {
    'K':       (f'{CG_DIR}/CG_023-WMAP_K_IQU_n0512_v1.fits',  'CG', 'mK', 'CG_WMAP_K', 22.8),
    'Ka':      (f'{CG_DIR}/CG_030-WMAP_Ka_IQU_n0512_v1.fits', 'CG', 'mK', 'CG_WMAP_Ka', 33.0),
    'Q1':      (f'{CG_DIR}/CG_040-WMAP_Q1_IQU_n0512_v1.fits', 'CG', 'mK', 'CG_WMAP_Q1', 40.7),
    'Q2':      (f'{CG_DIR}/CG_040-WMAP_Q2_IQU_n0512_v1.fits', 'CG', 'mK', 'CG_WMAP_Q2', 40.7),
    '030':     (f'{CG_DIR}/CG_030_IQU_n0512_v1.fits',         'CG', 'uK', 'CG_LFI_030', 28.4),
    '044':     (f'{CG_DIR}/CG_044_IQU_n0512_v1.fits',         'CG', 'uK', 'CG_LFI_044', 44.1),
    'MFI11':   (f'{MFI_DIR}/quijote_mfi_skymap_11ghz_512_dr1.fits',              'MFI', 'mK', 'QUIJOTE_MFI_311', 11.1),
    'MFI13':   (f'{MFI_DIR}/quijote_mfi_skymap_13ghz_512_dr1.fits',              'MFI', 'mK', 'QUIJOTE_MFI_313', 12.9),
    'MFI17H2': (f'{MFI_DIR}/horns/quijote_mfi_skymap_17ghz_horn2_512_dr1.fits',  'MFI', 'mK', 'QUIJOTE_MFI_217', 16.7),
    'MFI17H4': (f'{MFI_DIR}/horns/quijote_mfi_skymap_17ghz_horn4_512_dr1.fits',  'MFI', 'mK', 'QUIJOTE_MFI_417', 17.0),
    'MFI19H2': (f'{MFI_DIR}/horns/quijote_mfi_skymap_19ghz_horn2_512_dr1.fits',  'MFI', 'mK', 'QUIJOTE_MFI_219', 18.7),
    'MFI19H4': (f'{MFI_DIR}/horns/quijote_mfi_skymap_19ghz_horn4_512_dr1.fits',  'MFI', 'mK', 'QUIJOTE_MFI_419', 19.0),
}

# Random seed of each map
SEEDS = {'K': 60600, 'Ka': 60601, 'Q1': 60602, 'Q2': 60603, '030': 60604, '044': 60605,
         'MFI11': 60606, 'MFI13': 60607, 'MFI17H2': 60608, 'MFI17H4': 60609, 'MFI19H2': 60610, 'MFI19H4': 60611}

SOLAR_DIPOLE = (3362.08, 264.021, 48.253)        # Planck 2018: amplitude (uK_CMB), l, b (deg)


def solar_dipole_map(nside, unit):
    ''' Planck 2018 solar dipole in mK or uK '''

    amp, l, b = SOLAR_DIPOLE
    amp = {'mK': 1e-3, 'K': 1e-6, 'uK': 1.0}[unit]*amp
    return amp*np.dot(hp.ang2vec(l, b, lonlat=True), np.array(hp.pix2vec(nside, np.arange(12*nside**2))))


def load(band):
    ''' Reads a map and its per-pixel noise (standard deviations and QU covariance) '''

    fname, kind, unit, name, freq = DATASETS[band]

    if kind == 'CG':
        I, Q, U, sI, sQ, sU, cQU = (np.asarray(c, float) for c in hp.read_map(fname, field=range(7)))
        I = I - solar_dipole_map(hp.get_nside(I), unit)
        good = np.ones(len(I), bool)

    else:
        # Columns: I_STOKES, Q_STOKES, U_STOKES, NHITS_I, NHITS_QU, WEI_I, WEI_Q, WEI_U, COV_QU
        I, Q, U, wI, wQ, wU, cQU = (np.asarray(c, float) for c in hp.read_map(fname, field=(0, 1, 2, 5, 6, 7, 8)))
        good = (hp.read_map(MFI_MASK) > 0.5) & np.all(np.isfinite([I, Q, U, wI, wQ, wU, cQU]), axis=0) \
               & np.all(np.array([I, Q, U, cQU]) > -1e29, axis=0) & (wI > 0) & (wQ > 0) & (wU > 0)
        with np.errstate(divide='ignore', invalid='ignore'):
            sI, sQ, sU = 1/np.sqrt(wI), 1/np.sqrt(wQ), 1/np.sqrt(wU)

    for m in (I, Q, U, sI, sQ, sU, cQU):
        m[~good] = np.nan

    return dict(I=I, Q=Q, U=U, sI=sI, sQ=sQ, sU=sU, cQU=cQU, good=good, unit=unit, name=name, kind=kind, freq=freq)


def operator(I, Q, U, band, iter=0):
    ''' Smooths at the native resolution, then degrades to each Nside '''

    I, Q, U = smooth_iqu(I, Q, U, band, FWHM, iter=iter)
    return {n: np.array([hp.ud_grade(m, n) for m in (I, Q, U)]) for n in NSIDES}


def run(band, nsim, seed=None):
    ''' Smoothed map and noise of one band, written to PRODUCTS/<name>/ '''

    d = load(band)
    name, unit = d['name'], d['unit']
    out = os.path.join(PRODUCTS, name)
    os.makedirs(out, exist_ok=True)

    if seed is None:
        seed = SEEDS[band]
    rng = np.random.default_rng(seed)

    good = d['good']
    ng = good.sum()
    npix = len(good)
    keep = {n: hp.ud_grade(good.astype(float), n) > 1 - 1e-9 for n in NSIDES}   # output pixels with all native pixels usable

    # Cholesky factor of the QU covariance in each pixel
    L11 = d['sQ'][good]
    L21 = d['cQU'][good]/L11
    L22 = np.sqrt(np.clip(d['sU'][good]**2 - L21**2, 0, None))

    # Smoothed data
    data = operator(d['I'], d['Q'], d['U'], band, iter=3)

    # Noise realisations: running sums for the variances and the QU covariance
    S1 = {n: np.zeros((3, 12*n**2)) for n in NSIDES}
    S2 = {n: np.zeros((3, 12*n**2)) for n in NSIDES}
    SQU = {n: np.zeros(12*n**2) for n in NSIDES}
    sims64 = np.zeros((nsim, 3, 12*64**2), np.float32)

    print(f'{name}: {nsim} realisations, seed {seed}, {ng/npix:.1%} of the sky used', flush=True)
    t0 = time.time()
    for r in range(nsim):
        z = rng.standard_normal((3, ng))
        n = np.full((3, npix), np.nan)
        n[0, good] = d['sI'][good]*z[0]
        n[1, good] = L11*z[1]
        n[2, good] = L21*z[1] + L22*z[2]
        s = operator(n[0], n[1], n[2], band)
        for k in NSIDES:
            S1[k] += s[k]
            S2[k] += s[k]**2
            SQU[k] += s[k][1]*s[k][2]
        sims64[r] = s[64]
        if (r+1) % 50 == 0 or r == 0:
            elapsed = time.time()-t0
            print(f'  {r+1}/{nsim}  {elapsed:.0f} s, about {elapsed/(r+1)*(nsim-r-1)/60:.0f} min left', flush=True)

    # Write the maps with their errors
    for k in NSIDES:
        mean = S1[k]/nsim
        var = (S2[k] - nsim*mean**2)/(nsim-1)
        cov = (SQU[k] - nsim*mean[1]*mean[2])/(nsim-1)
        err = np.sqrt(np.clip(var, 0, None))
        dk = data[k].copy()
        dk[:, ~keep[k]] = hp.UNSEEN
        err[:, ~keep[k]] = hp.UNSEEN
        cov[~keep[k]] = hp.UNSEEN

        cols = [fits.Column(name=c, format='D', unit=u, array=a) for c, u, a in
                (('I', unit, dk[0]), ('Q', unit, dk[1]), ('U', unit, dk[2]),
                 ('I_ERR', unit, err[0]), ('Q_ERR', unit, err[1]), ('U_ERR', unit, err[2]),
                 ('QU_COV', f'{unit}^2', cov))]
        hdu = fits.BinTableHDU.from_columns(cols)
        for key, v, c in (('PIXTYPE', 'HEALPIX', ''), ('ORDERING', 'RING', ''), ('COORDSYS', 'G', 'Galactic'),
                          ('NSIDE', k, ''), ('INDXSCHM', 'IMPLICIT', ''), ('FWHM', FWHM, 'arcmin'),
                          ('FREQ', d['freq'], 'GHz'), ('BAND', band, ''),
                          ('SRCFILE', os.path.basename(DATASETS[band][0]), ''), ('NSIM', nsim, ''), ('SEED', seed, ''),
                          ('DIPOLE', d['kind'] == 'CG', 'Planck 2018 solar dipole subtracted from I'),
                          ('FDEC', False, 'no function of declination subtracted')):
            hdu.header[key] = (v, c)
        hdu.header['COMMENT'] = 'Smoothed at the native Nside with G_l/(b_l p_l) (spin 2 for Q, U), then ud_grade.'
        hdu.header['COMMENT'] = 'Errors and QU covariance from NSIM white-noise realisations through the same steps.'
        fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(f'{out}/{name}_60arcmin_n{k}.fits', overwrite=True)

    sims64[:, :, ~keep[64]] = np.nan
    np.save(f'{out}/{name}_60arcmin_n64_noise_sims.npy', sims64)

    # Check: the two halves of the realisations should give the same sigma within sqrt(2/(nsim/2))
    a = np.nanstd(sims64[:nsim//2, 1], axis=0)
    b = np.nanstd(sims64[nsim//2:, 1], axis=0)
    ok = keep[64] & (b > 0)
    print(f'{name}: done in {(time.time()-t0)/60:.1f} min. Split-half Q sigma ratio: mean '
          f'{np.mean(a[ok]/b[ok]):.4f}, std {np.std(a[ok]/b[ok]):.4f} (expected {np.sqrt(2/nsim):.4f})', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('band', choices=list(DATASETS))
    parser.add_argument('--nsim', type=int, default=NSIM)
    parser.add_argument('--seed', type=int, default=None)
    args = parser.parse_args()
    run(args.band, args.nsim, args.seed)
