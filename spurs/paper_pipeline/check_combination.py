''' check_combination.py:
Checks of the weighted maps at Nside 64:
  1) the GLS noise covariance (sum W)^-1 against the scatter of the combined
     noise realisations (each map's realisations pushed through the same weights)
  2) the difference between the full GLS solution and the diagonal formula
     (sum D W)(sum W)^-1 used in the original version of the paper

Writes PRODUCTS/check_combination.json.

Example usage:
python check_combination.py

Version 1.0 [Oct 2026]
Roke Cepeda-Arroita
roke.cepeda@iac.es
'''

import json
import numpy as np

from config import PRODUCTS
from combine import combine, COMBOS, INPUTS, scale_factor

NSIM = 1000
res = {}

for name in ('quijote', 'wmapplanck', 'full', 'full_1f'):
    keys, fd, of = COMBOS[name]
    out, ins = combine(keys, 64, use_fdec=fd, onef=of, return_inputs=True)
    npix = 12*64**2

    # Weights of each map, and the GLS covariance
    W = {}
    Wsum = np.zeros((npix, 2, 2))
    for k, (d, C, g, s) in ins.items():
        Wk = np.zeros_like(C)
        Wk[g] = np.linalg.inv(C[g])
        W[k] = Wk
        Wsum += Wk
    ok = out['ok']
    cov = np.zeros_like(Wsum)
    cov[ok] = np.linalg.inv(Wsum[ok])

    # Combine the noise realisations with the same weights (the 1/f option only changes the
    # weights; the realisations are white noise)
    acc = np.zeros((NSIM, npix, 2))
    for k in keys:
        nm, nu, band, det, u = INPUTS[k]
        sims = np.load(f'{PRODUCTS}/{nm}/{nm}_60arcmin_n64_noise_sims.npy', mmap_mode='r')
        s = scale_factor(k)*u
        n = np.nan_to_num(np.asarray(sims[:, 1:3, :], float).transpose(0, 2, 1))*s
        acc += np.einsum('pij,spj->spi', W[k], n)
    m = np.einsum('pij,spj->spi', cov, acc)
    vq_sim = m[:, :, 0].var(axis=0, ddof=1)
    sel = ok & (out['vQ'] > 0)
    r = vq_sim[sel]/out['vQ'][sel]

    # Diagonal formula of the original code, on the data
    b_diag = np.zeros((npix, 2, 2))
    for k, (d, C, g, s) in ins.items():
        D = np.zeros((npix, 2, 2))
        D[:, 0, 0], D[:, 1, 1] = d[:, 0], d[:, 1]
        b_diag += D @ W[k]
    md = np.einsum('pij,pjk->pik', b_diag[ok], cov[ok])
    dq = md[:, 0, 0] - out['Q'][ok]
    du = md[:, 1, 1] - out['U'][ok]
    rel = np.median(np.hypot(dq, du)/np.sqrt(out['vQ'][ok]))

    res[name] = dict(var_sim_over_analytic_median=float(np.median(r)), p16=float(np.percentile(r, 16)), p84=float(np.percentile(r, 84)),
                     mike_formula_dQU_over_sigma_median=float(rel),
                     mike_formula_dQU_over_sigma_p99=float(np.percentile(np.hypot(dq, du)/np.sqrt(out['vQ'][ok]), 99)))
    print(name, res[name], flush=True)

json.dump(res, open(f'{PRODUCTS}/check_combination.json', 'w'), indent=1)
