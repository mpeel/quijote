''' combine.py:
Weighted polarisation maps at 10 GHz from the smoothed maps and noise of noise_mc.py.

Each input map i is scaled to NU0 with a power law of index BETA and its colour
correction (fastcc, alpha = beta + 2):

    d_i = s_i (Q_i, U_i),    C_i = s_i^2 [[var_Q, cov_QU], [cov_QU, var_U]],    s_i = cc_i (NU0/nu_i)^beta

and the maps are combined pixel by pixel by generalised least squares:

    m = (sum W_i)^-1 sum W_i d_i,    W_i = C_i^-1,    cov(m) = (sum W_i)^-1

Weighted maps (PRODUCTS/combined/<name>_60arcmin_n<nside>.fits):
    quijote            311, 313, 217, 417, 219, 419
    wmapplanck         WMAP K, Ka, Q1, Q2 and LFI 28.4, 44.1 GHz
    wmapplanck_fdec    same, with FDEC subtracted from each input
    full               all six QUIJOTE maps (with the satellite-band mask) + WMAP/Planck with FDEC
    full_no1719        311, 313 + WMAP/Planck with FDEC
    *_nofdec           full combinations with no FDEC subtracted
    *_1f               QUIJOTE variances scaled by the 1/f factors at 1 deg (ONEF)
Also written: the WMAP Q band (GLS of Q1 and Q2 at 40.7 GHz, with and without FDEC)
and each WMAP/Planck map with FDEC subtracted, at its own frequency.

Example usage:
python combine.py              (Nside 512, 256 and 64)
python combine.py 64

Version 1.0 [Oct 2026]
Roke Cepeda-Arroita
roke.cepeda@iac.es
'''

import os
import sys
import numpy as np
import healpy as hp
from astropy.io import fits
from fastcc.fastcc import fastcc

from config import PRODUCTS, BETA, NU0, FWHM
from fdec import apply_fdec

UNS = hp.UNSEEN


# Input Maps

# key: (product name, frequency in GHz, fastcc band, fastcc detector, factor to mK)
INPUTS = {
    '311': ('QUIJOTE_MFI_311', 11.1, 'Q11', 'Q311', 1.0),
    '313': ('QUIJOTE_MFI_313', 12.9, 'Q13', 'Q313', 1.0),
    '217': ('QUIJOTE_MFI_217', 16.7, 'Q17', 'Q217', 1.0),
    '417': ('QUIJOTE_MFI_417', 17.0, 'Q17', 'Q417', 1.0),
    '219': ('QUIJOTE_MFI_219', 18.7, 'Q19', 'Q219', 1.0),
    '419': ('QUIJOTE_MFI_419', 19.0, 'Q19', 'Q419', 1.0),
    'K':   ('CG_WMAP_K',  22.8, 'WK',  None, 1.0),
    'Ka':  ('CG_WMAP_Ka', 33.0, 'WKa', None, 1.0),
    'Q1':  ('CG_WMAP_Q1', 40.7, 'WQ',  None, 1.0),
    'Q2':  ('CG_WMAP_Q2', 40.7, 'WQ',  None, 1.0),
    '030': ('CG_LFI_030', 28.4, 'P30', None, 1e-3),
    '044': ('CG_LFI_044', 44.1, 'P44', None, 1e-3),
}
QT = ['311', '313', '217', '417', '219', '419']       # QUIJOTE
WP = ['K', 'Ka', 'Q1', 'Q2', '030', '044']            # WMAP and Planck

# name: (inputs, FDEC subtracted from WMAP/Planck, QUIJOTE variances scaled by ONEF)
COMBOS = {
    'quijote':               (QT, False, False),
    'quijote_1f':            (QT, False, True),
    'wmapplanck':            (WP, False, False),
    'wmapplanck_fdec':       (WP, True, False),
    'full':                  (QT + WP, True, False),
    'full_1f':               (QT + WP, True, True),
    'full_no1719':           (['311', '313'] + WP, True, False),
    'full_no1719_1f':        (['311', '313'] + WP, True, True),
    'full_nofdec':           (QT + WP, False, False),
    'full_nofdec_1f':        (QT + WP, False, True),
    'full_no1719_nofdec':    (['311', '313'] + WP, False, False),
    'full_no1719_nofdec_1f': (['311', '313'] + WP, False, True),
}

# Excess noise at 1 deg from the correlated 1/f noise, sigma(real)/sigma(white) in (Q, U),
# measured in the DR1 half-mission null maps (H1-H2)/2 smoothed to 60 arcmin, Nside 64, |b| > 30 deg.
# On top of the DR1 Table-12 factors, which are already in the weights.
ONEF = {'311': (1.446, 1.506), '313': (1.435, 1.458), '217': (1.739, 1.712),
        '417': (1.412, 1.463), '219': (1.820, 1.814), '419': (1.760, 1.677)}


def scale_factor(key):
    ''' Colour correction times (NU0/nu)^beta '''

    name, nu, band, det, _ = INPUTS[key]
    if det:
        cc = fastcc(band, BETA + 2.0, detector=det)
    else:
        cc = fastcc(band, BETA + 2.0)
    return cc*(NU0/nu)**BETA


def read(key, nside):
    ''' Smoothed Q, U and their noise covariance, in mK '''

    name, nu, band, det, u = INPUTS[key]
    d = fits.getdata(f'{PRODUCTS}/{name}/{name}_60arcmin_n{nside}.fits', 1)
    Q, U, sQ, sU, c = (np.asarray(d[k], float) for k in ('Q', 'U', 'Q_ERR', 'U_ERR', 'QU_COV'))
    good = (Q > -1e29) & (sQ > 0) & (sU > 0) & (sQ < 1e29)
    return dict(Q=Q*u, U=U*u, vQ=sQ**2*u**2, vU=sU**2*u**2, c=c*u**2, good=good)


_fdec_cache = {}
def fdec_maps(key, nside):
    ''' FDEC templates (Q, U) of a map: computed at Nside 512, then degraded '''

    if (key, nside) not in _fdec_cache:
        m = read(key, 512)
        Qf, Uf, tq, tu = apply_fdec(np.where(m['good'], m['Q'], UNS), np.where(m['good'], m['U'], UNS))
        tq, tu = (hp.ud_grade(t, nside) for t in (tq, tu))
        _fdec_cache[(key, nside)] = (tq, tu)
    return _fdec_cache[(key, nside)]


def combine(keys, nside, use_fdec=False, scale=True, return_inputs=False, onef=False):
    ''' Per-pixel GLS combination of the maps in keys '''

    npix = 12*nside**2
    Wsum = np.zeros((npix, 2, 2))
    b = np.zeros((npix, 2))
    n = np.zeros(npix, int)
    ins = {}

    for k in keys:
        m = read(k, nside)
        Q, U = m['Q'].copy(), m['U'].copy()
        if use_fdec and k in WP:
            tq, tu = fdec_maps(k, nside)
            Q -= tq
            U -= tu
        s = scale_factor(k) if scale else 1.0
        g = m['good']

        # Noise covariance of the scaled map
        C = np.zeros((npix, 2, 2))
        C[:, 0, 0], C[:, 1, 1], C[:, 0, 1], C[:, 1, 0] = m['vQ'], m['vU'], m['c'], m['c']
        C *= s**2
        if onef and k in ONEF:
            rq, ru = ONEF[k]
            C[:, 0, 0] *= rq**2
            C[:, 1, 1] *= ru**2
            C[:, 0, 1] *= rq*ru
            C[:, 1, 0] *= rq*ru

        # Weights and weighted data
        W = np.zeros_like(C)
        W[g] = np.linalg.inv(C[g])
        d = np.stack([s*Q, s*U], axis=1)
        d[~g] = 0
        Wsum += W
        b += np.einsum('pij,pj->pi', W, d)
        n += g
        ins[k] = (d, C, g, s)

    # Solution and its covariance
    ok = n > 0
    cov = np.full((npix, 2, 2), np.nan)
    cov[ok] = np.linalg.inv(Wsum[ok])
    mm = np.full((npix, 2), np.nan)
    mm[ok] = np.einsum('pij,pj->pi', cov[ok], b[ok])
    out = dict(Q=mm[:, 0], U=mm[:, 1], vQ=cov[:, 0, 0], vU=cov[:, 1, 1], c=cov[:, 0, 1], n=n, ok=ok)

    if return_inputs:
        return out, ins
    return out


def write(out, path, keys, use_fdec, nside, nu_ref, onef=False):
    ''' Writes a weighted map: Q, U, P, Q_ERR, U_ERR, QU_COV and the number of input maps '''

    ok = out['ok']
    f = lambda a: np.where(ok, a, UNS)
    cols = [fits.Column(name=nm, format='D', unit=un, array=a) for nm, un, a in
            (('Q', 'mK_CMB', f(out['Q'])), ('U', 'mK_CMB', f(out['U'])),
             ('P', 'mK_CMB', f(np.hypot(out['Q'], out['U']))),
             ('Q_ERR', 'mK_CMB', f(np.sqrt(out['vQ']))), ('U_ERR', 'mK_CMB', f(np.sqrt(out['vU']))),
             ('QU_COV', 'mK_CMB^2', f(out['c'])), ('NMAPS', '', out['n'].astype(float)))]
    hdu = fits.BinTableHDU.from_columns(cols)
    for key, v, c in (('PIXTYPE', 'HEALPIX', ''), ('ORDERING', 'RING', ''), ('COORDSYS', 'G', ''),
                      ('NSIDE', nside, ''), ('INDXSCHM', 'IMPLICIT', ''), ('FWHM', FWHM, 'arcmin'),
                      ('REFFREQ', nu_ref, 'GHz'), ('BETA', BETA, 'spectral index used for the scaling'),
                      ('INPUTS', ','.join(keys), ''), ('FDEC', use_fdec, 'FDEC subtracted from the WMAP/Planck inputs'),
                      ('ONEF', onef, 'QUIJOTE variances x (1/f factor)^2')):
        hdu.header[key] = (v, c)
    hdu.header['COMMENT'] = 'Per-pixel GLS combination with the Monte Carlo QU covariances; colour corrections from fastcc.'
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(path, overwrite=True)


if __name__ == '__main__':
    nsides = [int(x) for x in sys.argv[1:]] or [512, 256, 64]
    od = f'{PRODUCTS}/combined'
    os.makedirs(od, exist_ok=True)

    for nside in nsides:

        # Weighted maps at 10 GHz
        for name, (keys, fd, of) in COMBOS.items():
            out = combine(keys, nside, use_fdec=fd, onef=of)
            write(out, f'{od}/{name}_60arcmin_n{nside}.fits', keys, fd, nside, NU0, onef=of)
            print(f'{name} n{nside}: {out["ok"].mean():.1%} of the sky', flush=True)

        # WMAP Q band (Q1 + Q2) at 40.7 GHz, for Table 1
        out = combine(['Q1', 'Q2'], nside, scale=False)
        write(out, f'{od}/wmapQ_60arcmin_n{nside}.fits', ['Q1', 'Q2'], False, nside, 40.7)

        # WMAP and Planck maps with FDEC subtracted, at their own frequency
        for k in ('K', 'Ka', 'Q1', 'Q2', '030', '044'):
            m = read(k, nside)
            tq, tu = fdec_maps(k, nside)
            o = dict(Q=m['Q']-tq, U=m['U']-tu, vQ=m['vQ'], vU=m['vU'], c=m['c'], n=m['good'].astype(int), ok=m['good'])
            write(o, f'{od}/{INPUTS[k][0]}_fdec_60arcmin_n{nside}.fits', [k], True, nside, INPUTS[k][1])

        # WMAP Q band with FDEC: each DA filtered, then combined
        out = combine(['Q1', 'Q2'], nside, use_fdec=True, scale=False)
        write(out, f'{od}/wmapQ_fdec_60arcmin_n{nside}.fits', ['Q1', 'Q2'], True, nside, 40.7)
