''' fdec.py:
Function-of-declination (FDEC) filtering of the WMAP and Planck maps, to remove
the same large-scale power that is removed from the QUIJOTE maps
(Rubino-Martin et al. 2023, Sect. 2.4.2):

  - Q and U are rotated to equatorial coordinates
  - f(dec) is the median of each iso-declination ring, excluding |b| < 10 deg
    and the 10 per cent brightest pixels in P
  - the f(dec) template is rotated back to Galactic and subtracted (only the
    template is rotated, so the data are not interpolated)

It is applied to the smoothed maps at Nside 512, and to the data only.

Example usage:
Q_filtered, U_filtered, tq, tu = apply_fdec(Q, U)

Version 1.0 [Oct 2026]
Roke Cepeda-Arroita
roke.cepeda@iac.es
'''

import numpy as np
import healpy as hp

R_GC = hp.Rotator(coord=['G', 'C'])
R_CG = hp.Rotator(coord=['C', 'G'])


def fdec_template(Q, U, bcut=10.0, bright_frac=0.10, extra_mask=None):
    ''' Galactic (Q, U) FDEC template for Galactic maps Q, U (UNSEEN or NaN = no data) '''

    nside = hp.get_nside(Q)
    npix = len(Q)

    # Rotate the maps and their bad pixels to equatorial coordinates
    bad = ~np.isfinite(Q) | ~np.isfinite(U) | (Q < -1e29) | (U < -1e29)
    q = np.where(bad, 0., Q)
    u = np.where(bad, 0., U)
    _, qe, ue = R_GC.rotate_map_pixel([np.zeros(npix), q, u])
    bad_e = R_GC.rotate_map_pixel(bad.astype(float)) > 0.01

    # Pixels used: away from the Galactic plane, and not among the brightest in P
    th, ph = hp.pix2ang(nside, np.arange(npix))
    thg, _ = R_CG(th, ph)
    use = (np.abs(90 - np.degrees(thg)) > bcut) & ~bad_e
    if extra_mask is not None:
        use &= R_GC.rotate_map_pixel(extra_mask.astype(float)) > 0.5
    P = np.hypot(qe, ue)
    use &= P < np.percentile(P[use], 100*(1 - bright_frac))

    # Median of each iso-declination ring
    fq = np.zeros(npix)
    fu = np.zeros(npix)
    for t in np.unique(th):
        ring = th == t
        sel = ring & use
        if sel.sum() >= 4:
            fq[ring] = np.median(qe[sel])
            fu[ring] = np.median(ue[sel])

    # Back to Galactic coordinates
    _, tq, tu = R_CG.rotate_map_pixel([np.zeros(npix), fq, fu])
    return tq, tu


def apply_fdec(Q, U, **kw):
    ''' Returns the filtered Q, U and the subtracted templates '''

    tq, tu = fdec_template(Q, U, **kw)
    return Q - tq, U - tu, tq, tu
