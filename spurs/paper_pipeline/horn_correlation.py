''' horn_correlation.py:
Effect on the weighted maps of the noise correlation between the two frequencies
of the same QUIJOTE horn (Sect. 2.1 of the paper). The per-pixel GLS of combine.py
is repeated with the full covariance of all the input maps, where the two maps of
a horn have correlated Q and U noise,

    C_ab = rho * C_a^(1/2) C_b^(1/2)

with rho = 0.3 (311/313), 0.4 (217/219) and 0.2 (417/419), and all other pairs
uncorrelated. With rho = 0 it gives exactly the maps of combine.py.
The S/N is computed as in Table 1 (sn_table.py). Writes PRODUCTS/horn_correlation.json.

Example usage:
python horn_correlation.py

Version 1.0 [Oct 2026]
Roke Cepeda-Arroita
roke.cepeda@iac.es
'''

import json
import numpy as np
from scipy.linalg import sqrtm

import combine as cb
import sn_table as st
from config import PRODUCTS

N = 64
RHO = {('311', '313'): 0.30, ('217', '219'): 0.40, ('417', '419'): 0.20}


def gls(keys, use_fdec, rho_scale):
    ''' GLS combination of the maps in keys with correlated noise between the two maps of each horn '''

    out, ins = cb.combine(keys, N, use_fdec=use_fdec, return_inputs=True)
    npix = 12*N**2
    K = len(keys)

    # Data vector and block covariance of all the inputs, pixel by pixel
    d = np.zeros((npix, 2*K))
    C = np.zeros((npix, 2*K, 2*K))
    g = np.zeros((npix, K), bool)
    for i, k in enumerate(keys):
        di, Ci, gi, _ = ins[k]
        d[:, 2*i:2*i+2] = di
        C[:, 2*i:2*i+2, 2*i:2*i+2] = Ci
        g[:, i] = gi

    # Cross terms between the two frequencies of each horn
    sq = {}
    for (a, b), rho in RHO.items():
        if a in keys and b in keys and rho_scale > 0:
            ia, ib = keys.index(a), keys.index(b)
            both = g[:, ia] & g[:, ib]
            for p in np.where(both)[0]:
                Sa = sq.setdefault((a, p), np.real(sqrtm(C[p, 2*ia:2*ia+2, 2*ia:2*ia+2])))
                Sb = sq.setdefault((b, p), np.real(sqrtm(C[p, 2*ib:2*ib+2, 2*ib:2*ib+2])))
                X = rho*rho_scale*Sa @ Sb
                C[p, 2*ia:2*ia+2, 2*ib:2*ib+2] = X
                C[p, 2*ib:2*ib+2, 2*ia:2*ia+2] = X.T

    # GLS solution in each pixel, using the maps with data there
    Q = np.full(npix, np.nan)
    U = Q.copy()
    vQ = Q.copy()
    vU = Q.copy()
    for p in np.where(g.any(axis=1))[0]:
        sel = np.repeat(g[p], 2)
        Ci = np.linalg.inv(C[p][np.ix_(sel, sel)])
        A = np.tile(np.eye(2), (sel.sum()//2, 1))
        cov = np.linalg.inv(A.T @ Ci @ A)
        m = cov @ A.T @ Ci @ d[p, sel]
        Q[p], U[p], vQ[p], vU[p] = m[0], m[1], cov[0, 0], cov[1, 1]

    return Q, U, np.sqrt(vQ), np.sqrt(vU), out


if __name__ == '__main__':
    r = st.region()
    res = {}
    for name, keys, fd in (('QUIJOTE combination', cb.QT, False), ('Full combination', cb.QT + cb.WP, True)):
        q0, u0, sq0, su0, out = gls(list(keys), fd, 0.0)
        assert np.nanmax(np.abs(q0[out['ok']] - out['Q'][out['ok']])) < 1e-9       # rho = 0 gives the combine.py maps
        q1, u1, sq1, su1, _ = gls(list(keys), fd, 1.0)
        s0 = st.sn(q0, u0, sq0, su0, r)
        s1 = st.sn(q1, u1, sq1, su1, r)
        err = float(np.nanmedian(sq1[r]/sq0[r])), float(np.nanmedian(su1[r]/su0[r]))
        res[name] = dict(sn_uncorrelated=s0, sn_correlated=s1, sigma_ratio_median_QU=err)
        print(f'{name}: S/N Q {s0[0]:.2f} ({s0[1]:.2f}) -> {s1[0]:.2f} ({s1[1]:.2f}),  U {s0[2]:.2f} ({s0[3]:.2f}) -> {s1[2]:.2f} ({s1[3]:.2f});'
              f'  median sigma ratio Q {err[0]:.3f}, U {err[1]:.3f}', flush=True)

    json.dump(dict(rho={f'{a}/{b}': v for (a, b), v in RHO.items()}, results=res),
              open(f'{PRODUCTS}/horn_correlation.json', 'w'), indent=1)
