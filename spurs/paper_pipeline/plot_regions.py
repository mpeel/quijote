''' plot_regions.py:
Top panel of Fig. 3: the full weighted map at 10 GHz (histogram-equalised colour
scale) with the best-fit circle arcs of the known loops (black) and of the arcs
and shells identified in the paper (red), the outline of the Fermi bubbles and
the labels of the main regions.

Example usage:
python plot_regions.py

Version 1.0 [Oct 2026]
Roke Cepeda-Arroita
roke.cepeda@iac.es
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
import healpy as hp
import csv
import glob
import os
from scipy.optimize import minimize
from cmcrameri import cm as cmc

from config import PRODUCTS, FIGURES, LOOP_REGION_DIR, NEW_REGION_DIR

os.makedirs(f'{FIGURES}/fig3', exist_ok=True)
OUT = f'{FIGURES}/fig3/fig6_regions'
regs_dir = LOOP_REGION_DIR
new_regs_dir = NEW_REGION_DIR

# ---------------------------------------------------------------------------
# Name -> region-file number mapping  (derived from check_loop_labels.py)
# ---------------------------------------------------------------------------
LOOP_FILE_ID = {
    'GCS':  1,
    'I':    2,
    'IX':   3,
    'X':    4,
    'III':  5,
    'IIIs': 6,
    'VIIb': 7,
    'Is':   8,
    'VIII': 9,
    'XIV':  10,
    'XIII': 11,
    'XI':   12,
    'XII':  13,
    'IV':   14,
    'II':   15,
}

# ---------------------------------------------------------------------------
# New spurs/shells: filename stem -> label text, position, colour, and
# whether to fit a circle (fit_circle=True) or just draw raw points (False).
# Edit 'label' and 'loc' (lon, lat) for each new region as needed.
# ---------------------------------------------------------------------------
NEW_REGIONS = {
    'new_loop_1':  {'file': 'new_loop_1.reg',  'label': '', 'loc': (0, 0), 'color': 'r',  'fit_circle': True,  'full_circle': False, 'stroke': False, 'zorder': None, 'alpha': 1.0, 'linestyle': '--', 'linewidth': 0.8, 'smooth': None},
    'new_loop_2':  {'file': 'new_loop_2.reg',  'label': '', 'loc': (0, 0), 'color': 'r',  'fit_circle': True,  'full_circle': False, 'stroke': False, 'zorder': None, 'alpha': 1.0, 'linestyle': '--', 'linewidth': 0.8, 'smooth': None},
    'new_loop_3':  {'file': 'new_loop_3.reg',  'label': '', 'loc': (0, 0), 'color': 'r',  'fit_circle': True,  'full_circle': False, 'stroke': False, 'zorder': None, 'alpha': 1.0, 'linestyle': '--', 'linewidth': 0.8, 'smooth': None},
    'new_loop_5':  {'file': 'new_loop_5.reg',  'label': '', 'loc': (0, 0), 'color': 'r',  'fit_circle': True,  'full_circle': False, 'stroke': False, 'zorder': None, 'alpha': 1.0, 'linestyle': '--', 'linewidth': 0.8, 'smooth': None},
    'new_loop_6':  {'file': 'new_loop_6.reg',  'label': '', 'loc': (0, 0), 'color': 'r',  'fit_circle': True,  'full_circle': False, 'stroke': False, 'zorder': None, 'alpha': 1.0, 'linestyle': '--', 'linewidth': 0.8, 'smooth': None},
    'new_shell_1': {'file': 'new_shell_1.reg', 'label': '', 'loc': (0, 0), 'color': 'r',  'fit_circle': True,  'full_circle': True,  'stroke': False, 'zorder': None, 'alpha': 1.0, 'linestyle': '--', 'linewidth': 0.8, 'smooth': None},
    'new_shell_2': {'file': 'new_shell_2.reg', 'label': '', 'loc': (0, 0), 'color': 'r',  'fit_circle': True,  'full_circle': True,  'stroke': False, 'zorder': None, 'alpha': 1.0, 'linestyle': '--', 'linewidth': 0.8, 'smooth': None},
    'new_shell_3': {'file': 'new_shell_3.reg', 'label': '', 'loc': (0, 0), 'color': 'r',  'fit_circle': True,  'full_circle': False, 'stroke': False, 'zorder': None, 'alpha': 1.0, 'linestyle': '--', 'linewidth': 0.8, 'smooth': None},
    # Fermi bubbles: two halves drawn separately, smoothed, no clipping
    'fermi_n':     {'file': f'{regs_dir}/amespur6.reg', 'label': '', 'loc': (0, 0), 'color': 'w',  'fit_circle': False, 'full_circle': False, 'stroke': False, 'zorder': 0, 'alpha': 0.35, 'linestyle': '-',  'linewidth': 0.8, 'smooth': 0.2},
    'fermi_s':     {'file': f'{regs_dir}/amespur7.reg', 'label': '', 'loc': (0, 0), 'color': 'w',  'fit_circle': False, 'full_circle': False, 'stroke': False, 'zorder': 0, 'alpha': 0.35, 'linestyle': '-',  'linewidth': 0.8, 'smooth': 0.2},
}

# Arc extension for new regions: name -> (extend_start_deg, extend_end_deg)
# Only used when full_circle=False
NEW_REGION_ARC = {k: (0, 0) for k in NEW_REGIONS}

# ---------------------------------------------------------------------------
# Arc extension: loop name -> (extend_start_deg, extend_end_deg)
# ---------------------------------------------------------------------------
LOOP_ARC = {
    'GCS':  (20, 0),
    'I':    (20, 10),
    'IX':   (0, 10),
    'X':    (0, 0),
    'III':  (10, 5),
    'IIIs': (0, 0),
    'VIIb': (0, 0),
    'Is':   (-10, 5),
    'VIII': (0, 0),
    'XIV':  (5, 0),
    'XIII': (15, 0),
    'XI':   (-40, 0),
    'XII':  (5, 0),
    'IV':   (20, -30),
    'II':   (-40, 10),
}

# Label positions (Galactic lon, lat) for each loop name
LOOP_LABEL_LOC = {
    'GCS':  (-18,  24),
    'I':    ( 43,  45),
    'IX':   ( -15,  51),
    'X':    ( 52,  14),
    'III':  (115,  30),
    'IIIs': (120, -40),
    'VIIb': ( 25, -69),
    'Is':   (-25, -34),
    'VIII': (-20, -17),
    'XIV':  (-41,  20),
    'XIII': (-33,   3),
    'XI':   (-125+7, 5),
    'XII':  (-108,-57),
    'IV':   (-41,  35),
    'II':   (140,  -65),
}

# ---------------------------------------------------------------------------
# Geometry utilities
# ---------------------------------------------------------------------------

def lb_to_vec(l_deg, b_deg):
    l = np.radians(np.asarray(l_deg, dtype=float))
    b = np.radians(np.asarray(b_deg, dtype=float))
    return np.array([np.cos(b)*np.cos(l),
                     np.cos(b)*np.sin(l),
                     np.sin(b)])


def fit_circle(lons, lats):
    """Fit a small circle to (lon, lat) points. Returns (l0, b0, r) in degrees."""
    pts = lb_to_vec(lons, lats)

    def cost(params):
        c = lb_to_vec(params[0], params[1])
        ang = np.degrees(np.arccos(np.clip(c @ pts, -1, 1)))
        return np.sum((ang - ang.mean())**2)

    mv = pts.mean(axis=1); mv /= np.linalg.norm(mv)
    res = minimize(cost,
                   [np.degrees(np.arctan2(mv[1], mv[0])) % 360,
                    np.degrees(np.arcsin(np.clip(mv[2], -1, 1)))],
                   method="Nelder-Mead",
                   options={"xatol": 1e-4, "fatol": 1e-6, "maxiter": 10000})
    l0, b0 = res.x[0] % 360, res.x[1]
    c = lb_to_vec(l0, b0)
    r = np.degrees(np.arccos(np.clip(c @ pts, -1, 1))).mean()
    return l0, b0, r


def make_basis(l0_deg, b0_deg):
    """Return centre unit vector c and orthonormal basis (u, v) perpendicular to it."""
    c = lb_to_vec(l0_deg, b0_deg)
    north = np.array([0., 0., 1.])
    if abs(np.dot(c, north)) > 0.99:
        north = np.array([1., 0., 0.])
    u = north - np.dot(north, c) * c; u /= np.linalg.norm(u)
    v = np.cross(c, u);               v /= np.linalg.norm(v)
    return c, u, v


def pt_angle(lon, lat, c, u, v):
    """Parametric angle of a sky point projected onto the circle plane."""
    vec = lb_to_vec(float(lon), float(lat))
    perp = vec - np.dot(vec, c) * c
    if np.linalg.norm(perp) < 1e-10:
        return 0.0
    perp /= np.linalg.norm(perp)
    return np.arctan2(np.dot(perp, v), np.dot(perp, u))


def make_arc(l0, b0, r, lons, lats,
             extend_start_deg=0.0, extend_end_deg=0.0,
             n=500, full_circle=False):
    """
    Return (arc_lons, arc_lats) for a circle arc.

    If full_circle=True, draws the complete 360° circle regardless of the
    data point distribution.

    Otherwise the arc spans the true angular min/max of the data points
    (sorted, not file-order), then extended by extend_start_deg before the
    first point and extend_end_deg after the last.
    """
    c, u, v = make_basis(l0, b0)
    r_rad   = np.radians(r)

    if full_circle:
        angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    else:
        # Compute parametric angles for every data point
        angles_data = np.array([pt_angle(lo, la, c, u, v)
                                 for lo, la in zip(lons, lats)])
        angles_data = np.unwrap(angles_data)

        # Diagnostic: print actual span so you can spot ordering issues
        span_deg = np.degrees(angles_data.max() - angles_data.min())
        print(f"    arc span from data: {span_deg:.1f}°  "
              f"(first={np.degrees(angles_data[0]):.1f}°, "
              f"last={np.degrees(angles_data[-1]):.1f}°, "
              f"min={np.degrees(angles_data.min()):.1f}°, "
              f"max={np.degrees(angles_data.max()):.1f}°)")

        # Use true min/max so unordered files still produce the correct arc
        a_start = angles_data.min() - np.radians(extend_start_deg)
        a_end   = angles_data.max() + np.radians(extend_end_deg)
        angles  = np.linspace(a_start, a_end, n)

    pts = (np.cos(r_rad) * c[:, None]
           + np.sin(r_rad) * np.cos(angles) * u[:, None]
           + np.sin(r_rad) * np.sin(angles) * v[:, None])
    pts /= np.linalg.norm(pts, axis=0)

    x, y, z = pts
    arc_lons = np.degrees(np.arctan2(y, x)) % 360
    arc_lats = np.degrees(np.arcsin(np.clip(z, -1, 1)))

    return arc_lons, arc_lats


def projplot_arc(arc_lons, arc_lats, *args, **kwargs):
    """
    Wrapper around hp.visufunc.projplot that splits the arc at any lon
    wraparound (jumps > 180°) to avoid a spurious line across the map.
    """
    split_idx    = np.where(np.abs(np.diff(arc_lons)) > 180)[0] + 1
    segments_lon = np.split(arc_lons, split_idx)
    segments_lat = np.split(arc_lats, split_idx)
    for slon, slat in zip(segments_lon, segments_lat):
        hp.visufunc.projplot(slon, slat, *args, lonlat=True, **kwargs)


def read_region_file(loop_name):
    """Read coordinate points for a named loop from its .reg file."""
    def parse(s):
        if ':' in s:
            p = [float(x) for x in s.split(':')]
            return p[0] + p[1]/60 + p[2]/3600
        return float(s)

    file_id = LOOP_FILE_ID[loop_name]
    matches = glob.glob(f'{regs_dir}/short_{file_id}xyrad*.reg')

    if not matches:
        return None, None

    lons, lats = [], []
    with open(matches[0]) as f:
        for line in csv.reader(f, delimiter='\t'):
            parts = [x for x in line[0].split(' ') if x]
            lons.append(parse(parts[0]))
            lats.append(parse(parts[1]))
    return lons, lats


def read_new_region_file(filename):
    """Read coordinate points from a .reg file.
    If filename is an absolute path or contains a directory separator it is
    used as-is; otherwise it is looked up under new_regs_dir."""
    def parse(s):
        if ':' in s:
            p = [float(x) for x in s.split(':')]
            return p[0] + p[1]/60 + p[2]/3600
        return float(s)

    import os
    filepath = filename if os.sep in filename or filename.startswith('.') else f'{new_regs_dir}/{filename}'
    matches = glob.glob(filepath)
    if not matches:
        return None, None

    lons, lats = [], []
    with open(matches[0]) as f:
        for line in csv.reader(f, delimiter='\t'):
            parts = [x for x in line[0].split(' ') if x]
            if len(parts) < 2:
                continue
            lons.append(parse(parts[0]))
            lats.append(parse(parts[1]))
    return lons, lats




def smooth_lonlat(lons, lats, smooth_factor=0.0, n_out=500):
    """
    Smooth an arbitrary (lon, lat) polyline using a cubic spline.

    Parameters
    ----------
    lons, lats    : input coordinate arrays (degrees)
    smooth_factor : passed to scipy.interpolate.splprep as `s`.
                    0  -> interpolating spline (passes through every point)
                    >0 -> smoothing spline
    n_out         : number of output points along the smoothed curve

    Returns
    -------
    lons_s, lats_s : smoothed coordinate arrays
    """
    from scipy.interpolate import splprep, splev

    lons = np.array(lons, dtype=float)
    lats = np.array(lats, dtype=float)

    # Work in 3-D unit-vector space to avoid longitude wraparound issues
    vecs = lb_to_vec(lons, lats)          # shape (3, N)
    x, y, z = vecs

    tck, _ = splprep([x, y, z], s=smooth_factor, k=3, per=False)
    u_new  = np.linspace(0, 1, n_out)
    xs, ys, zs = splev(u_new, tck)

    # Renormalise onto the unit sphere
    norms = np.sqrt(xs**2 + ys**2 + zs**2)
    xs, ys, zs = xs/norms, ys/norms, zs/norms

    lons_s = np.degrees(np.arctan2(ys, xs)) % 360
    lats_s = np.degrees(np.arcsin(np.clip(zs, -1, 1)))
    return lons_s, lats_s


def find_polyline_intersections(lons1, lats1, lons2, lats2):
    """
    Find all intersection points between two polylines in lon/lat space.
    Returns a list of (lon, lat, i1, t1, i2, t2) where i1/i2 are the segment
    indices and t1/t2 are fractional positions along each segment.
    """
    def seg_intersect(p1, p2, p3, p4):
        d1 = p2 - p1;  d2 = p4 - p3
        cross = d1[0]*d2[1] - d1[1]*d2[0]
        if abs(cross) < 1e-10:
            return None
        d3 = p3 - p1
        t = (d3[0]*d2[1] - d3[1]*d2[0]) / cross
        u = (d3[0]*d1[1] - d3[1]*d1[0]) / cross
        if 0 <= t <= 1 and 0 <= u <= 1:
            return t, u
        return None

    pts1 = np.stack([lons1, lats1], axis=1)
    pts2 = np.stack([lons2, lats2], axis=1)
    hits = []
    for i in range(len(pts1) - 1):
        for j in range(len(pts2) - 1):
            r = seg_intersect(pts1[i], pts1[i+1], pts2[j], pts2[j+1])
            if r is not None:
                t, u = r
                pt = pts1[i] + t * (pts1[i+1] - pts1[i])
                hits.append((pt[0], pt[1], i, t, j, u))
    return hits


def trim_fermi_curves(lons1, lats1, lons2, lats2):
    """
    Trim two overlapping curves so they meet exactly at their intersection
    points instead of crossing past each other.

    Keeps the 'outer' arc of each curve (the section between the two
    intersection points that stays away from the partner curve).
    Falls back to the original arrays if fewer than 2 intersections are found.
    """
    hits = find_polyline_intersections(lons1, lats1, lons2, lats2)
    if len(hits) < 2:
        print(f"    trim_fermi: only {len(hits)} intersection(s) found, skipping trim")
        return lons1, lats1, lons2, lats2

    # Sort hits by position along curve 1 and take the outermost two
    hits.sort(key=lambda h: h[2] + h[3])   # i1 + t1 ~ continuous parameter
    h_start, h_end = hits[0], hits[-1]

    def insert_and_trim(lons, lats, h_a, h_b, param_idx):
        """Keep the segment of `lons/lats` between intersections h_a and h_b
        (sorted by their parameter along this curve), inserting the exact
        intersection coordinate at each end."""
        # Sort the two hits by their parameter along *this* curve
        i_a, t_a = h_a[param_idx], h_a[param_idx + 1]
        i_b, t_b = h_b[param_idx], h_b[param_idx + 1]
        if i_a + t_a > i_b + t_b:
            (i_a, t_a), (i_b, t_b) = (i_b, t_b), (i_a, t_a)
            h_a, h_b = h_b, h_a

        lon_a = lons[i_a] + t_a * (lons[i_a+1] - lons[i_a])
        lat_a = lats[i_a] + t_a * (lats[i_a+1] - lats[i_a])
        lon_b = lons[i_b] + t_b * (lons[i_b+1] - lons[i_b])
        lat_b = lats[i_b] + t_b * (lats[i_b+1] - lats[i_b])

        new_lons = np.concatenate([[lon_a], lons[i_a+1:i_b+1], [lon_b]])
        new_lats = np.concatenate([[lat_a], lats[i_a+1:i_b+1], [lat_b]])
        return new_lons, new_lats

    # param_idx=2 → (i1,t1) for curve1;  param_idx=4 → (i2,t2) for curve2
    new_lons1, new_lats1 = insert_and_trim(lons1, lats1, h_start, h_end, 2)
    new_lons2, new_lats2 = insert_and_trim(lons2, lats2, h_start, h_end, 4)
    return new_lons1, new_lats1, new_lons2, new_lats2


# ---------------------------------------------------------------------------
# Load map
# ---------------------------------------------------------------------------
# Full weighted map (combine.py), P at 10 GHz, Nside 512
fits_file = f'{PRODUCTS}/combined/full_60arcmin_n512.fits'
print(f"Reading map: {fits_file}")
m = hp.read_map(fits_file, field=2)   # columns Q, U, P

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
fig = plt.figure(1, figsize=(6.5, 4))
hp.mollview(m, fig=1, unit=r'$\mu$K', cbar=False, title=None,
            cmap=cmc.roma_r, norm='hist', xsize=800*3)

# ---------------------------------------------------------------------------
# Best-fit arcs and labels — existing loops
# ---------------------------------------------------------------------------
print()
print(f"{'Loop':<6}  {'l0':>6}  {'b0':>6}  {'r':>5}")
print("-" * 35)
for name in LOOP_FILE_ID:
    lons, lats = read_region_file(name)
    loc = LOOP_LABEL_LOC[name]

    if lons is None:
        print(f"{name:<6}  (no region file found, skipping arc)")
        # Still plot the label even without a region file
        hp.visufunc.projtext(loc[0], loc[1], name, lonlat=True, fontsize=10)
        continue

    l0, b0, r  = fit_circle(lons, lats)
    print(f"{name:<6}  {l0:>6.1f}  {b0:>+6.1f}  {r:>5.1f}")

    ext_start, ext_end = LOOP_ARC[name]
    arc_lons, arc_lats = make_arc(l0, b0, r, lons, lats,
                                  extend_start_deg=ext_start,
                                  extend_end_deg=ext_end)

    projplot_arc(arc_lons, arc_lats, 'k', linewidth=0.7, linestyle='--')
    hp.visufunc.projtext(loc[0], loc[1], name, lonlat=True, fontsize=10)

# ---------------------------------------------------------------------------
# New spurs / shells
# ---------------------------------------------------------------------------
print()
print(f"{'New region':<12}  {'l0':>6}  {'b0':>6}  {'r':>5}")
print("-" * 40)

# ---------------------------------------------------------------------------
# Pre-process Fermi bubbles: smooth both curves then trim at their crossings
# ---------------------------------------------------------------------------
_fermi_cache = {}
for _fname in ('fermi_n', 'fermi_s'):
    _info = NEW_REGIONS[_fname]
    _lons, _lats = read_new_region_file(_info['file'])
    if _lons is not None and _info.get('smooth') is not None:
        _arr_l = np.array(_lons, dtype=float)
        _arr_b = np.array(_lats, dtype=float)
        _arr_l, _arr_b = smooth_lonlat(_arr_l, _arr_b,
                                       smooth_factor=_info['smooth'] * len(_lons))
        _fermi_cache[_fname] = (_arr_l, _arr_b)

if len(_fermi_cache) == 2:
    print("    trimming Fermi bubble curves at their intersections...")
    _l_n, _b_n = _fermi_cache['fermi_n']
    _l_s, _b_s = _fermi_cache['fermi_s']
    _l_n, _b_n, _l_s, _b_s = trim_fermi_curves(_l_n, _b_n, _l_s, _b_s)
    _fermi_cache['fermi_n'] = (_l_n, _b_n)
    _fermi_cache['fermi_s'] = (_l_s, _b_s)

for name, info in NEW_REGIONS.items():
    loc       = info['loc']
    label     = info['label']
    color     = info['color']
    do_fit    = info['fit_circle']
    full_circ = info.get('full_circle', False)
    do_stroke = info['stroke']
    zorder    = info['zorder']
    alpha     = info.get('alpha', 1.0)
    linestyle = info.get('linestyle', '--')
    linewidth = info.get('linewidth', 0.7)
    smooth    = info.get('smooth', None)

    lons, lats = read_new_region_file(info['file'])

    if lons is None:
        print(f"{name:<12}  (no region file found, skipping arc)")
        if label:
            hp.visufunc.projtext(loc[0], loc[1], label, lonlat=True,
                                 fontsize=10, color=color)
        continue

    if do_fit:
        # Best-fit circle arc (or full circle)
        l0, b0, r = fit_circle(lons, lats)
        circ_label = "full circle" if full_circ else "arc"
        print(f"{name:<12}  {l0:>6.1f}  {b0:>+6.1f}  {r:>5.1f}  ({circ_label})")
        ext_start, ext_end = NEW_REGION_ARC.get(name, (0, 0))
        arc_lons, arc_lats = make_arc(l0, b0, r, lons, lats,
                                      extend_start_deg=ext_start,
                                      extend_end_deg=ext_end,
                                      full_circle=full_circ)
        split_idx    = np.where(np.abs(np.diff(arc_lons)) > 180)[0] + 1
        segments_lon = np.split(arc_lons, split_idx)
        segments_lat = np.split(arc_lats, split_idx)
        for slon, slat in zip(segments_lon, segments_lat):
            kw = dict(linewidth=linewidth, linestyle=linestyle, alpha=alpha)
            if zorder is not None:
                kw['zorder'] = zorder
            line_objs = hp.visufunc.projplot(slon, slat, color, lonlat=True, **kw)
            if do_stroke:
                for ln in line_objs:
                    ln.set_path_effects([
                        PathEffects.withStroke(linewidth=1.8, foreground='k', alpha=0.6)
                    ])
    else:
        # Raw point-to-point (Fermi bubbles etc.)
        print(f"{name:<12}  (raw, no circle fit)")
        raw_lons = np.array(lons, dtype=float)
        raw_lats = np.array(lats, dtype=float)

        if name in _fermi_cache:
            # Use the pre-smoothed and intersection-trimmed arrays
            raw_lons, raw_lats = _fermi_cache[name]
        elif smooth is not None:
            raw_lons, raw_lats = smooth_lonlat(raw_lons, raw_lats,
                                               smooth_factor=smooth * len(lons))

        split_idx    = np.where(np.abs(np.diff(raw_lons)) > 180)[0] + 1
        segments_lon = np.split(raw_lons, split_idx)
        segments_lat = np.split(raw_lats, split_idx)
        for slon, slat in zip(segments_lon, segments_lat):
            kw = dict(linewidth=linewidth, linestyle=linestyle, alpha=alpha)
            if zorder is not None:
                kw['zorder'] = zorder
            line_objs = hp.visufunc.projplot(slon, slat, color, lonlat=True, **kw)
            if do_stroke:
                for ln in line_objs:
                    ln.set_path_effects([
                        PathEffects.withStroke(linewidth=1.8, foreground='k', alpha=0.6)
                    ])

    if label:
        hp.visufunc.projtext(loc[0], loc[1], label, lonlat=True,
                             fontsize=10, color=color)

# ---------------------------------------------------------------------------
# Region labels
# ---------------------------------------------------------------------------
t1 = hp.visufunc.projtext(148,   0,    'Fan Region', lonlat=True,
                           fontsize=6, color='w', horizontalalignment='center')
t2 = hp.visufunc.projtext(-4.5, -47.7+4, 'Fermi\nBubbles', lonlat=True,
                           fontsize=6, color='w', horizontalalignment='center')
t3 = hp.visufunc.projtext(76.18+10, 5.75-4, 'Cyg', lonlat=True,
                           fontsize=5, color='w', horizontalalignment='center')
t4 = hp.visufunc.projtext(255.8+8, -4.1, 'Gum', lonlat=True,
                           fontsize=5, color='w', horizontalalignment='center')
# Hα filament label — placed just below Loop II
t5 = hp.visufunc.projtext(90, -65, r'H$\alpha$', lonlat=True,
                           fontsize=7, color='w', horizontalalignment='center')
# M42
t6 = hp.visufunc.projtext(209, -19.+2, r'M42', lonlat=True,
                           fontsize=5, color='w', horizontalalignment='center')

# Tau A
t7 = hp.visufunc.projtext(185, -5+2, r'Tau A', lonlat=True,
                           fontsize=4.5, color='w', horizontalalignment='center')


HA_ARROW_ANGLE_DEG = 110
HA_ARROW_LENGTH    = 7
angle_rad = np.radians(HA_ARROW_ANGLE_DEG)
tail_lon  = 80
tail_lat  = -58
tip_lon   = tail_lon + HA_ARROW_LENGTH * np.cos(angle_rad)
tip_lat   = tail_lat + HA_ARROW_LENGTH * np.sin(angle_rad)

ax_hp = plt.gca()
tail_proj = hp.projector.MollweideProj().ang2xy(tail_lon, tail_lat, lonlat=True)
tip_proj  = hp.projector.MollweideProj().ang2xy(tip_lon,  tip_lat,  lonlat=True)

arrow = ax_hp.annotate('',
    xy=tip_proj, xytext=tail_proj,
    xycoords='data', textcoords='data',
    arrowprops=dict(arrowstyle='->', color='w', lw=0.8, mutation_scale=6),
)
arrow.arrow_patch.set_path_effects([
    PathEffects.withStroke(linewidth=1.5, foreground='k', alpha=0.5)
])


for t in [t1, t2, t3, t4, t5, t6, t7]:
    t.set_path_effects([PathEffects.withStroke(linewidth=0.7, foreground='k', alpha=0.5)])



plt.savefig(OUT + '.png', dpi=600, bbox_inches='tight', pad_inches=0, transparent=True)
plt.savefig(OUT + '.pdf',   dpi=600,    bbox_inches='tight', pad_inches=0, transparent=True)
print("\nSaved: ./regions.png  ./regions.pdf")
