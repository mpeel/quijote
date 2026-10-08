''' make_figures.py:
Map figures of the paper, from the products of noise_mc.py and combine.py.

    fig1      QUIJOTE maps of the six horns and the QUIJOTE weighted map (P and angle)
    fig2      WMAP K, QUIJOTE, WMAP+Planck and full weighted maps on the same scale
    fig4      polarisation fraction with respect to Haslam 408 MHz (beta = -3.1), Nside 64
    figA1     WMAP+Planck weighted map, effect of FDEC, and differences with the previous version
    figB1     full weighted map: P, angle, Q and U
    figC1     QUIJOTE minus WMAP+Planck difference maps
    figmasks  region masks of Table 2 over the full weighted map

All single maps are scaled to 10 GHz with the same colour corrections used in the
weighted maps, so every panel is on the same scale. Figure 3 is made with
plot_regions.py and plot_lic.py.

Example usage:
python make_figures.py              (all figures)
python make_figures.py fig1 fig2

Version 1.0 [Oct 2026]
Roke Cepeda-Arroita
roke.cepeda@iac.es
'''

import os
import sys
import numpy as np
import healpy as hp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from astropy.io import fits
from cmcrameri import cm as cmc

from config import PRODUCTS, FIGURES, QT_MASK_PLOT, HASLAM, REGION_DIR, ORIGINAL_MAPS
from combine import scale_factor, INPUTS

NSIDE = 256
VMAX_P, VMAX_D = 1.5, 0.5                            # mK, colour scale limits of P and of the differences
CMAP_P, CMAP_ANG, CMAP_D = cmc.roma_r, cmc.romaO, cmc.vik
MK = r'$\mathrm{mK}_{\mathrm{CMB}}$'


# Reading the Maps

def nanify(a):
    ''' No data -> NaN '''

    a = np.asarray(a, float).copy()
    a[(a < -1e29) | (np.abs(a) > 1e20)] = np.nan
    return a


def single(key, nside=NSIDE, fdec=False):
    ''' Q, U, sigma_Q, sigma_U of one input map, scaled to 10 GHz '''

    name = INPUTS[key][0]
    if fdec:
        f = f'{PRODUCTS}/combined/{name}_fdec_60arcmin_n{nside}.fits'
    else:
        f = f'{PRODUCTS}/{name}/{name}_60arcmin_n{nside}.fits'
    d = fits.getdata(f, 1)
    u = 1.0 if fdec else INPUTS[key][4]
    s = scale_factor(key)*u
    return [nanify(d[c])*s for c in ('Q', 'U', 'Q_ERR', 'U_ERR')]


def combo(name, nside=NSIDE):
    ''' Q, U, sigma_Q, sigma_U of a weighted map '''

    d = fits.getdata(f'{PRODUCTS}/combined/{name}_60arcmin_n{nside}.fits', 1)
    return [nanify(d[c]) for c in ('Q', 'U', 'Q_ERR', 'U_ERR')]


def original_P(key):
    ''' P of a weighted map of the previous version of the paper '''

    m = hp.read_map(ORIGINAL_MAPS[key], field=None)
    m = np.array([nanify(x) for x in np.atleast_2d(m)])
    p = P(m[1], m[2]) if len(m) >= 3 else m[0]
    p[p == 0] = np.nan
    return hp.ud_grade(p, NSIDE)


def mask_qt(m):
    ''' Blanks the pixels outside the QUIJOTE mask '''

    mk = hp.ud_grade(hp.read_map(QT_MASK_PLOT), hp.npix2nside(len(m)))
    m = m.copy()
    m[mk < 0.5] = np.nan
    return m


def P(q, u):
    return np.hypot(q, u)


def ANG(q, u):
    return 0.5*np.degrees(np.arctan2(u, q))


def mas(q, u, sq, su):
    ''' Modified asymptotic estimator of P (Plaszczynski et al. 2014) '''

    p = np.hypot(q, u)
    b2 = ((q*su)**2 + (u*sq)**2)/p**2
    return p - b2*(1 - np.exp(-p**2/b2))/(2*p)


# Plotting

def save_map(m, path, cmap, vmin, vmax, **kw):
    ''' Mollweide map with no colour bar (the colour bars are separate files) '''

    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig = plt.figure(figsize=(6.5, 3.5))
    hp.mollview(m, fig=fig.number, cmap=cmap, min=vmin, max=vmax, cbar=False, title=None,
                xsize=2400, notext=True, badcolor='lightgray', **kw)
    plt.gca().set_frame_on(False)
    plt.savefig(path, dpi=600, bbox_inches='tight', pad_inches=0, transparent=True)
    plt.close(fig)
    print('  saved', path)


def save_cbar(path, cmap, vmin, vmax, label, ticks=None, ticklabels=None):
    ''' Horizontal colour bar '''

    fig, ax = plt.subplots(figsize=(4.5, 0.45))
    fig.subplots_adjust(left=0.05, right=0.95, bottom=0.55, top=0.95)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=mcolors.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cb = fig.colorbar(sm, cax=ax, orientation='horizontal')
    cb.set_label(label, fontsize=12, labelpad=-4)
    cb.set_ticks(ticks or [vmin, vmax])
    cb.set_ticklabels(ticklabels or [f'{vmin:g}', f'{vmax:g}'])
    cb.ax.tick_params(labelsize=11, length=0)
    fig.savefig(path, dpi=600, bbox_inches='tight', pad_inches=0.02, transparent=True)
    plt.close(fig)
    print('  saved', path)


def save_appendix(m, path, cmap, vmin, vmax, unit=MK):
    ''' Appendix maps: Mollweide map with the healpy colour bar '''

    os.makedirs(os.path.dirname(path), exist_ok=True)
    with matplotlib.rc_context({'font.size': 16}):
        hp.mollview(np.where(np.isfinite(m), m, hp.UNSEEN), min=vmin, max=vmax, cmap=cmap, unit=unit, title='')
        for ax in plt.gcf().axes[1:]:
            for t in ax.texts:
                t.set_fontsize(16)
            ax.xaxis.label.set_fontsize(16)
            ax.tick_params(labelsize=14)
        plt.savefig(path)
    plt.close('all')
    print('  saved', path)


# Figures

def fig1():
    od = f'{FIGURES}/fig1'
    for key, lab in (('311', 'fig1_mfi311_pol'), ('313', 'fig1_mfi313_pol'), ('417', 'fig1_mfi417_pol'),
                     ('419', 'fig1_mfi419_pol'), ('217', 'fig1_mfi217_pol'), ('219', 'fig1_mfi219_pol')):
        q, u, _, _ = single(key)
        save_map(mask_qt(P(q, u)), f'{od}/{lab}.pdf', CMAP_P, 0, VMAX_P)
    q, u, _, _ = combo('quijote')
    save_map(mask_qt(P(q, u)), f'{od}/fig2_mfi_combine.pdf', CMAP_P, 0, VMAX_P)
    save_map(mask_qt(ANG(q, u)), f'{od}/fig2_mfi_combine_polang.pdf', CMAP_ANG, -90, 90)
    save_cbar(f'{od}/colorbar_mkcmb.pdf', CMAP_P, 0, VMAX_P, MK, ticklabels=['0.0', f'{VMAX_P:.1f}'])
    save_cbar(f'{od}/colorbar_polang.pdf', CMAP_ANG, -90, 90, r'$\mathrm{deg}$', ticklabels=['-90', '+90'])


def fig2():
    od = f'{FIGURES}/fig2'
    q, u, _, _ = single('K')
    save_map(P(q, u), f'{od}/wmap9_x10.pdf', CMAP_P, 0, VMAX_P)
    q, u, _, _ = combo('quijote')
    save_map(mask_qt(P(q, u)), f'{od}/quijote.pdf', CMAP_P, 0, VMAX_P)
    q, u, _, _ = combo('wmapplanck')
    save_map(P(q, u), f'{od}/wmapplanck.pdf', CMAP_P, 0, VMAX_P)
    q, u, _, _ = combo('full')
    save_map(P(q, u), f'{od}/wmapplanck_quijote.pdf', CMAP_P, 0, VMAX_P)
    save_cbar(f'{od}/colorbar.pdf', CMAP_P, 0, VMAX_P, MK, ticklabels=['0.0', f'{VMAX_P:.1f}'])


def fig4():
    od = f'{FIGURES}/polfrac'
    has = hp.ud_grade(hp.read_map(HASLAM), NSIDE)*1e3*(10/0.408)**(-3.1)      # K -> mK, scaled to 10 GHz

    def frac(q, u, sq, su, qt=False):
        f = 100*mas(q, u, sq, su)/has
        f[~np.isfinite(f) | (f < 0)] = np.nan
        f = hp.ud_grade(np.where(np.isfinite(f), f, hp.UNSEEN), 64)
        f = nanify(f)
        return mask_qt(f) if qt else f

    save_map(frac(*single('311'), qt=True), f'{od}/mfi11.pdf', CMAP_P, 0, 50)
    save_map(frac(*combo('quijote'), qt=True), f'{od}/mfiweighted.pdf', CMAP_P, 0, 50)
    save_map(frac(*single('K')), f'{od}/wmapk9.pdf', CMAP_P, 0, 50)
    save_map(frac(*single('030')), f'{od}/p30.pdf', CMAP_P, 0, 50)
    save_cbar(f'{od}/colorbar.pdf', CMAP_P, 0, 50, r'$\%$')


def figA1():
    od = f'{FIGURES}/appendix'
    q, u, _, _ = combo('wmapplanck')
    pw = P(q, u)
    save_appendix(pw, f'{od}/figA1_cg_wmapplanck.pdf', CMAP_P, 0, VMAX_P)
    qf, uf, _, _ = combo('wmapplanck_fdec')
    save_appendix(P(qf, uf) - pw, f'{od}/figA1_cg_wmapplanck_fdec_minus_nofdec.pdf', CMAP_D, -VMAX_D, VMAX_D)
    save_appendix(P(qf, uf) - original_P('wmapplanck'), f'{od}/figA1_cg_minus_original_wmapplanck.pdf', CMAP_D, -VMAX_D, VMAX_D)
    q, u, _, _ = combo('full')
    save_appendix(P(q, u) - original_P('full'), f'{od}/figA1_cg_minus_original_full.pdf', CMAP_D, -VMAX_D, VMAX_D)


def figB1():
    od = f'{FIGURES}/appendix'
    q, u, _, _ = combo('full')
    save_appendix(P(q, u), f'{od}/figB1_full_P.pdf', CMAP_P, 0, VMAX_P)
    save_appendix(ANG(q, u), f'{od}/figB1_full_polang.pdf', CMAP_ANG, -90, 90, unit='deg')
    save_appendix(q, f'{od}/figB1_full_Q.pdf', CMAP_D, -1.0, 1.0)
    save_appendix(u, f'{od}/figB1_full_U.pdf', CMAP_D, -1.0, 1.0)


def figC1():
    od = f'{FIGURES}/appendix'
    qw, uw, _, _ = combo('wmapplanck_fdec')
    for key in ('311', '313', '417', '419'):
        q, u, _, _ = single(key)
        save_appendix(P(q, u) - P(qw, uw), f'{od}/figC1_mfi{key}_subP.pdf', CMAP_D, -VMAX_D, VMAX_D)
        save_appendix(np.hypot(q - qw, u - uw), f'{od}/figC1_mfi{key}_subP2.pdf', CMAP_P, 0, VMAX_D)


def figmasks():
    od = f'{FIGURES}/masks'
    q, u, _, _ = combo('full', 512)
    under = P(q, u)
    reg = hp.read_map(f'{REGION_DIR}/mask_quijote_mfi_horn3_intensity_satband_nside512.fits') != 0

    for stem in ('loop1_only', 'fan_mask', 'loop3_mask', 'Loop3_south_mask', 'loop3_shell_correct', 'shell_south_of_cygnus',
                 'loop_2', 'loop9_only', 'loop_between_9_and_GCS_only', 'loop_GCS_only', 'regionx', 'LoopLB_mask', 'Loop11_mask'):
        mk = hp.read_map(f'{REGION_DIR}/{stem}.fits')
        if hp.get_nside(mk) != 512:
            mk = hp.ud_grade(mk, 512)
        m = under.copy()
        m[~(reg & (mk != 0))] = np.nan
        vmax = round(float(np.ceil(np.nanpercentile(m, 99)*10)/10), 1)

        fig = plt.figure(figsize=(6.5, 4.0))
        hp.mollview(m, fig=fig.number, cmap=CMAP_P, min=0, max=vmax, cbar=False, title=None, xsize=2400, notext=True)
        for ax in fig.get_axes():
            pos = ax.get_position()
            ax.set_position([pos.x0, pos.y0 + 0.01, pos.width, pos.height*0.90])
            ax.set_frame_on(False)
        cax = fig.add_axes([0.29, 0.03, 0.42, 0.035])
        sm = plt.cm.ScalarMappable(cmap=CMAP_P, norm=mcolors.Normalize(vmin=0, vmax=vmax))
        sm.set_array([])
        cb = fig.colorbar(sm, cax=cax, orientation='horizontal')
        cb.set_label(MK, fontsize=14, labelpad=-4)
        cb.set_ticks([0, vmax])
        cb.set_ticklabels(['0.0', f'{vmax:.1f}'])
        cb.ax.tick_params(labelsize=13, length=3)
        os.makedirs(od, exist_ok=True)
        fig.savefig(f'{od}/{stem}.pdf', dpi=600, bbox_inches='tight', pad_inches=0.05, transparent=True)
        plt.close(fig)
        print('  saved', f'{od}/{stem}.pdf')


if __name__ == '__main__':
    todo = sys.argv[1:] or ['fig1', 'fig2', 'fig4', 'figA1', 'figB1', 'figC1', 'figmasks']
    for f in todo:
        print(f)
        globals()[f]()
