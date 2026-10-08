''' sn_table.py:
Signal-to-noise ratios of Table 1: median and mean of |Q|/sigma_Q and |U|/sigma_U
per pixel at Nside 64, over the QUIJOTE footprint (the pixels with data in all
six QUIJOTE maps, i.e. inside the satellite-band mask), for every map.

WMAP and Planck values are given with FDEC subtracted (as they enter the full
combination) and without. Writes PRODUCTS/sn_table.json and prints the table rows.

Example usage:
python sn_table.py

Version 1.0 [Oct 2026]
Roke Cepeda-Arroita
roke.cepeda@iac.es
'''

import json
import numpy as np
from astropy.io import fits

from config import PRODUCTS
from combine import ONEF

N = 64
HORNS = ['311', '313', '217', '417', '219', '419']
NAMES = {k: f'QUIJOTE_MFI_{k}' for k in HORNS}


def get(path):
    ''' Q, U, Q_ERR, U_ERR of a map '''

    d = fits.getdata(path, 1)
    return [np.asarray(d[c], float) for c in ('Q', 'U', 'Q_ERR', 'U_ERR')]


def region():
    ''' Pixels with data in all six QUIJOTE maps '''

    r = np.ones(12*N**2, bool)
    for k in HORNS:
        r &= get(f'{PRODUCTS}/{NAMES[k]}/{NAMES[k]}_60arcmin_n{N}.fits')[2] > 0
    return r


def sn(q, u, sq, su, r):
    ''' Median and mean of |Q|/sigma_Q and |U|/sigma_U in region r '''

    out = []
    for x, s in ((q, sq), (u, su)):
        v = np.abs(x[r])/s[r]
        out += [float(np.median(v)), float(np.mean(v))]
    return out


if __name__ == '__main__':
    r = region()
    print(f'Region: {r.sum()} Nside-64 pixels ({r.mean():.1%} of the sky)')
    rows = {}

    # QUIJOTE, with and without the 1/f factors
    for k in HORNS:
        q, u, sq, su = get(f'{PRODUCTS}/{NAMES[k]}/{NAMES[k]}_60arcmin_n{N}.fits')
        rows[f'QUIJOTE {k}'] = sn(q, u, sq, su, r)
        rows[f'QUIJOTE {k} (1/f)'] = sn(q, u, sq*ONEF[k][0], su*ONEF[k][1], r)

    # WMAP and Planck, with and without FDEC
    for lab, name in (('WMAP K', 'CG_WMAP_K'), ('Planck LFI 28.4', 'CG_LFI_030'),
                      ('WMAP Ka', 'CG_WMAP_Ka'), ('Planck LFI 44.1', 'CG_LFI_044')):
        rows[lab] = sn(*get(f'{PRODUCTS}/combined/{name}_fdec_60arcmin_n{N}.fits'), r)
        rows[lab + ' (no FDEC)'] = sn(*get(f'{PRODUCTS}/{name}/{name}_60arcmin_n{N}.fits'), r)
    rows['WMAP Q'] = sn(*get(f'{PRODUCTS}/combined/wmapQ_fdec_60arcmin_n{N}.fits'), r)
    rows['WMAP Q (no FDEC)'] = sn(*get(f'{PRODUCTS}/combined/wmapQ_60arcmin_n{N}.fits'), r)

    # Weighted maps
    for lab, c in (('QUIJOTE combination', 'quijote'), ('QUIJOTE combination (1/f)', 'quijote_1f'),
                   ('WMAP+Planck combination', 'wmapplanck_fdec'),
                   ('WMAP+Planck combination (no FDEC)', 'wmapplanck'), ('Full combination', 'full'),
                   ('Full combination (1/f)', 'full_1f'), ('Full combination, 11/13 GHz only', 'full_no1719'),
                   ('Full combination, 11/13 GHz only (1/f)', 'full_no1719_1f')):
        rows[lab] = sn(*get(f'{PRODUCTS}/combined/{c}_60arcmin_n{N}.fits'), r)

    json.dump(dict(npix=int(r.sum()), rows=rows), open(f'{PRODUCTS}/sn_table.json', 'w'), indent=1)
    for lab, v in rows.items():
        print(f'{lab:40s} & {v[0]:.2f} ({v[1]:.2f}) & {v[2]:.2f} ({v[3]:.2f})')
