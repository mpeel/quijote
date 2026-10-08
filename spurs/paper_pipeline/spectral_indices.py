''' spectral_indices.py:
Polarised spectral indices of the regions in Table 2 (Sect. 5.2, Fig. 5 and Appendix D).

For each region (region mask x QUIJOTE satellite-band and NCP mask, at Nside 64):
  - beta from the Fuskeland et al. (2014) method between QUIJOTE 11.1 GHz (horn 3)
    and Planck 28.4 GHz, and between QUIJOTE 11.1 GHz and WMAP 22.8 GHz
  - T-T plots in Q and U between 11.1 and 28.4 GHz, as a cross-check
All fits are colour corrected iteratively (ttplot_fuskeland.linear_fit_cc).

Writes FIGURES/<name>_tt_311_30_{q,u}.png, FIGURES/<label>_fusk.png,
FIGURES/<name>_wmapk_fusk.png and PRODUCTS/spectral_indices.json, and prints the
rows of Table 2.

Example usage:
python spectral_indices.py

Version 1.0 [Oct 2026]
Roke Cepeda-Arroita
roke.cepeda@iac.es
'''

import os
import json
import numpy as np
import healpy as hp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits
from cmcrameri import cm as cmc

from config import PRODUCTS, FIGURES, QT_MASK, REGION_DIR
from ttplot_fuskeland import p_fusk14_cc, linear_fit_cc, LAST

NSIDE = 64
M = REGION_DIR

# Table 2 row: (region mask, name of the T-T plots, name of the Fuskeland plot)
REGIONS = {
    'Loop I (NPS)':         (f'{M}/loop1_only.fits', 'loop1', 'loop1_planck'),
    'Fan region':           (f'{M}/fan_mask.fits', 'fan', 'fan_planck'),
    'Loop III':             (f'{M}/loop3_mask.fits', 'loop3', 'loop3_planck'),
    'Loop III South':       (f'{M}/Loop3_south_mask.fits', 'loop3south', 'loop3south_planck'),
    'Shell above Loop III': (f'{M}/loop3_shell_correct.fits', 'loop3_shell', 'loop3_shell_planck'),
    'Shell near Cygnus':    (f'{M}/shell_south_of_cygnus.fits', 'shell_south_of_cygnus', 'shell_south_of_cygnus'),
    'Loop II':              (f'{M}/loop_2.fits', 'loop_2', 'loop_2_planck'),
    'Loop IX':              (f'{M}/loop9_only.fits', 'loop9', 'loop9_planck'),
    'Loop between IX/GCS':  (f'{M}/loop_between_9_and_GCS_only.fits', 'loop_between_9_and_cgs', 'loop_between_9_and_cgs_planck'),
    'GCS':                  (f'{M}/loop_GCS_only.fits', 'cgs', 'cgs_planck'),
    'Loop X':               (f'{M}/regionx.fits', 'regionx', 'regionx_planck'),
    'Loop LB':              (f'{M}/LoopLB_mask.fits', 'looplb', 'looplb_planck'),
    'Loop XI':              (f'{M}/Loop11_mask.fits', 'loop11', 'loop11_planck'),
}

# Plot style, matching the map figures
C_PT, C_FIT, INK, MUTED, GRID = cmc.roma_r(0.12), cmc.roma_r(0.88), '#1f1f1f', '#5f5f5f', '#e4e4e4'
STYLE = {'font.size': 12, 'font.family': 'serif', 'mathtext.fontset': 'cm', 'axes.formatter.use_mathtext': False,
         'axes.edgecolor': MUTED, 'axes.linewidth': 0.8, 'axes.labelcolor': INK,
         'xtick.color': MUTED, 'ytick.color': MUTED, 'xtick.direction': 'in', 'ytick.direction': 'in',
         'xtick.top': True, 'ytick.right': True, 'legend.frameon': False, 'savefig.dpi': 300}


def load(name, unit=1.0):
    ''' I, Q, U and their variances at Nside 64 (no data = UNSEEN, variance 0) '''

    d = fits.getdata(f'{PRODUCTS}/{name}/{name}_60arcmin_n{NSIDE}.fits', 1)
    m = np.array([np.asarray(d[c], float) for c in ('I', 'Q', 'U')])
    v = np.array([np.asarray(d[c], float)**2 for c in ('I_ERR', 'Q_ERR', 'U_ERR')])
    bad = (m[1] < -1e29) | (v[1] <= 0) | (v[1] > 1e29)
    m *= unit
    v *= unit**2
    m[:, bad] = hp.UNSEEN
    v[:, bad] = 0.0
    return m, v


def region_mask(f):
    ''' Region mask times the QUIJOTE mask at Nside 64 (pixels fully inside both) '''

    r = hp.ud_grade(hp.read_map(f).astype(float), NSIDE)
    q = hp.ud_grade(hp.read_map(QT_MASK).astype(float), NSIDE)
    return ((r > 1 - 1e-6) & (q > 1 - 1e-6)).astype(float)


def _axes():
    fig, ax = plt.subplots(figsize=(5.2, 4.0))
    ax.grid(color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    return fig, ax


def plot_tt(x, y, sx, sy, f1, f2, xlabel, ylabel, path, band2, det2=False):
    ''' T-T plot of the colour-corrected data with the fitted line. Returns beta and its error '''

    beta, ebeta, c, ec, xc, yc, sxc, syc = linear_fit_cc(x, y, sx, sy, f1, f2, 'Q11', band2, 'Q311', det2, 0.001)
    m = (f2/f1)**beta

    with matplotlib.rc_context(STYLE):
        fig, ax = _axes()
        ax.errorbar(xc, yc, xerr=sxc, yerr=syc, fmt='o', ms=2.2, color=C_PT, ecolor=C_PT, elinewidth=0.5,
                    capsize=1.2, capthick=0.5, alpha=0.55, zorder=2, lw=0)
        t = np.linspace(xc.min(), xc.max(), 10)
        ax.plot(t, m*t + c, color=C_FIT, lw=2.0, zorder=3, label=rf'$\beta = {beta:.2f} \pm {abs(ebeta):.2f}$')
        ax.axhline(0, color=MUTED, lw=0.5, zorder=1)
        ax.axvline(0, color=MUTED, lw=0.5, zorder=1)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.legend(loc='upper left')
        fig.tight_layout()
        fig.savefig(path)
        plt.close(fig)

    return float(beta), float(abs(ebeta))


def plot_fusk(path, f2label):
    ''' beta as a function of the angle alpha from the last Fuskeland fit, with the weighted
        mean and its two uncertainties (minimum single-angle error and scatter over angles) '''

    r = LAST
    with matplotlib.rc_context(STYLE):
        fig, ax = _axes()
        ok = np.isfinite(r['sigma']) & (r['beta'] != 0)
        ax.axhspan(r['beta_tot'] - r['std'], r['beta_tot'] + r['std'], color=C_FIT, alpha=0.12, lw=0,
                   label=r'$1\sigma$ (scatter over $\alpha$)')
        ax.axhspan(r['beta_tot'] - r['smin'], r['beta_tot'] + r['smin'], color=C_FIT, alpha=0.28, lw=0,
                   label=r'$1\sigma$ (min. error)')
        ax.axhline(r['beta_tot'], color=C_FIT, lw=1.8,
                   label=rf'$\langle\beta\rangle = {r["beta_tot"]:.2f} \pm {max(r["std"], r["smin"]):.2f}$')
        ax.errorbar(r['alpha'][ok], r['beta'][ok], yerr=r['sigma'][ok], fmt='o', ms=4, color=C_PT,
                    elinewidth=1.2, capsize=2.5, capthick=1.2, zorder=3, lw=1.6, ls='-')
        ax.set_xlabel(r'$\alpha$ [deg]')
        ax.set_ylabel(r'$\beta_{11.1-' + f2label + r'\,\mathrm{GHz}}$')
        ax.set_xlim(-3, 88)
        ax.legend(loc='lower left', bbox_to_anchor=(0, 1.01, 1, 0.2), mode='expand', ncol=3, fontsize=8.5, handlelength=1.4)
        fig.tight_layout()
        fig.savefig(path)
        plt.close(fig)


if __name__ == '__main__':
    matplotlib.rcParams.update(STYLE)
    os.makedirs(FIGURES, exist_ok=True)

    # Maps at Nside 64, in mK
    m311, v311 = load('QUIJOTE_MFI_311')
    p30, vp30 = load('CG_LFI_030', unit=1e-3)
    wk, vwk = load('CG_WMAP_K')

    res = {}
    for reg, (mf, name, flab) in REGIONS.items():
        mask = region_mask(mf)
        sel = (mask == 1) & (m311[1] != hp.UNSEEN) & (p30[1] != hp.UNSEEN)
        print(f'=== {reg}: {sel.sum()} Nside-64 pixels', flush=True)

        # T-T plots in Q and U, 11.1 vs 28.4 GHz
        tt = {}
        for X, j in (('q', 1), ('u', 2)):
            tt[X] = plot_tt(m311[j][sel], p30[j][sel], np.sqrt(v311[j][sel]), np.sqrt(vp30[j][sel]), 11.1, 28.4,
                            rf'QUIJOTE 11.1 GHz ${X.upper()}$ [mK$_{{\rm CMB}}$]', rf'Planck 28.4 GHz ${X.upper()}$ [mK$_{{\rm CMB}}$]',
                            f'{FIGURES}/{name}_tt_311_30_{X}.png', 'P30')

        # Fuskeland method, 11.1-28.4 and 11.1-22.8 GHz
        b30 = p_fusk14_cc(mask, m311, v311, 11.1, 'Q11', 'Q311', p30, vp30, 28.4, 'P30', False, 0.0, FIGURES, flab, make_plot=False)
        plot_fusk(f'{FIGURES}/{flab}_fusk.png', '28.4')
        bk = p_fusk14_cc(mask, m311, v311, 11.1, 'Q11', 'Q311', wk, vwk, 22.8, 'WK', False, 0.0, FIGURES, name + '_wmapk', make_plot=False)
        plot_fusk(f'{FIGURES}/{name}_wmapk_fusk.png', '22.8')

        # Quoted error: the larger of the minimum single-angle error and the scatter over angles
        res[reg] = dict(npix=int(sel.sum()), beta_30=float(b30[0]), err_30=float(b30[2]),
                        beta_K=float(bk[0]), err_K=float(bk[2]), tt_q=tt['q'], tt_u=tt['u'])
        print(f'    11-23: {bk[0]:.2f} +- {bk[2]:.2f}    11-28.4: {b30[0]:.2f} +- {b30[2]:.2f}', flush=True)

    json.dump(res, open(f'{PRODUCTS}/spectral_indices.json', 'w'), indent=1)
    for reg, r in res.items():
        print(f"  {reg:22s} & ${r['beta_K']:.2f} \\pm {r['err_K']:.2f}$ & ${r['beta_30']:.2f} \\pm {r['err_30']:.2f}$")
