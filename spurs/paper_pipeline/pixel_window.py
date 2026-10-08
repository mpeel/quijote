''' pixel_window.py:
Effect of deconvolving the native HEALPix pixel window in the smoothing kernel
(Sect. 3.1 of the paper). Each map is smoothed with and without the pixel window,
together with five noise realisations, and compared at Nside 64.

Example usage:
python pixel_window.py

Version 1.0 [Oct 2026]
Roke Cepeda-Arroita
roke.cepeda@iac.es
'''

import numpy as np
import healpy as hp

from config import FWHM
from noise_mc import load
from smoothing import smooth_iqu

rng = np.random.default_rng(5)

for band in ('MFI11', 'K', '030'):
    d = load(band)
    g = d['good']
    res = {}

    for dp in (True, False):

        # Signal
        I, Q, U = smooth_iqu(d['I'], d['Q'], d['U'], band, FWHM, iter=1, deconvolve_pixwin=dp)
        Q64, U64 = (hp.ud_grade(np.where(np.isfinite(m), m, hp.UNSEEN), 64) for m in (Q, U))

        # Noise in Q, five realisations
        v = []
        for _ in range(5):
            n = np.where(g, d['sQ']*rng.standard_normal(len(g)), np.nan)
            _, nq, _ = smooth_iqu(np.zeros_like(n), n, np.zeros_like(n), band, FWHM, deconvolve_pixwin=dp)
            v.append(hp.ud_grade(np.where(np.isfinite(nq), nq, hp.UNSEEN), 64))
        res[dp] = (Q64, U64, np.var(np.array(v), axis=0))

    ok = (res[True][0] > -1e29) & (res[False][0] > -1e29)
    P1 = np.hypot(res[True][0], res[True][1])
    P0 = np.hypot(res[False][0], res[False][1])
    hi = ok & (P1 > np.percentile(P1[ok], 50))
    print(f'{band}: P(with/without pixel window) median {np.median(P1[hi]/P0[hi]):.4f};  '
          f'noise variance ratio median {np.median(res[True][2][ok]/res[False][2][ok]):.4f}', flush=True)
