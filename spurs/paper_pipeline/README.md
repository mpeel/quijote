# Polarised Loops and Spurs: Weighted Maps, Noise and Spectral Indices

Code used to produce the results of *QUIJOTE scientific results – XI. Polarized Synchrotron Loops and Spurs in Variance-Weighted Microwave Survey Maps* (Peel, Cepeda-Arroita et al. 2026). It smooths the QUIJOTE MFI, WMAP and *Planck* LFI maps to 1°, propagates their noise with Monte Carlo simulations, combines them into inverse-variance weighted polarisation maps at 10 GHz, and measures the polarised spectral indices of the loops and spurs.

---

## **Overview**

### Main Components:
1. **`config.py`**: paths to the input maps and the main settings (resolution, Nsides, number of noise realisations, spectral index).
2. **`smoothing.py`**: smooths I, Q, U to 60 arcmin with the full beam transfer function of each map (Q and U as a spin-2 field).
3. **`noise_mc.py`**: smoothed maps at Nside 512, 256 and 64, and their noise from 1000 white-noise realisations per map.
4. **`fdec.py`**: function-of-declination filtering of the WMAP and *Planck* maps, to remove the same large-scale power as in QUIJOTE.
5. **`combine.py`**: weighted Q and U maps at 10 GHz (per-pixel generalised least squares with the full QU noise covariance).
6. **`sn_table.py`**: signal-to-noise ratios of Table 1.
7. **`spectral_indices.py`**: spectral indices of Table 2 (Fuskeland et al. 2014 method), T-T plots and β(α) plots.
8. **`make_figures.py`**, **`plot_regions.py`**, **`plot_lic.py`**: the map figures of the paper.
9. **`check_combination.py`**, **`horn_correlation.py`**, **`pixel_window.py`**: checks quoted in the paper.
10. **`ttplot_fuskeland.py`**: copy of Mike Peel's [astrocode](https://github.com/mpeel/astrocode) T-T fitting code, with the changes listed at the top of the file.

---

## **Dependencies**
- `numpy`, `scipy`, `healpy`, `astropy`, `matplotlib`, `cmcrameri`
- `fastcc` for the colour corrections: clone [mpeel/fastcc](https://github.com/mpeel/fastcc) and add the folder that contains it to your `PYTHONPATH` (the code uses `from fastcc.fastcc import fastcc`)

---

## **Input Data**

Set `DATA_DIR` in `config.py`. The code expects:

| Folder | Contents |
|---|---|
| `cosmoglobe/` | Cosmoglobe DR1 maps `CG_023-WMAP_K`, `CG_030-WMAP_Ka`, `CG_040-WMAP_Q1`, `CG_040-WMAP_Q2`, `CG_030`, `CG_044` ([Watts et al. 2023](https://arxiv.org/abs/2303.08095)) |
| `quijote_mfi/` | QUIJOTE MFI DR1 maps at 11 and 13 GHz, the horn 2 and horn 4 maps at 17 and 19 GHz in `horns/`, and the satellite-band mask in `masks/` ([doi:10.26698/quijote-mfi-dr1](https://doi.org/10.26698/quijote-mfi-dr1)) |
| `beams/` | WMAP 9-year `wmap_ampl_bl_*_9yr_v5p1.txt`, NPIPE `Bl_TEB_npipe6v19_*GHzx*GHz.fits`, and the QUIJOTE RIMO `quijote_mfi/rimo_quijote_mfi_beamtf_dr1.fits` |
| `quijote_masks/` | QUIJOTE masks and the region masks of Table 2 |
| `ancillary/` | Haslam 408 MHz map (Remazeilles et al. 2015) |
| `regions/` | loop traces of Vidal et al. (2015) and of the new arcs, for Fig. 3 |
| `original/` | weighted maps of the previous version of the paper (only for Fig. A1) |

---

## **Getting Started**

1. **Set the paths** in `config.py`.

2. **Smooth the maps and propagate the noise** (about an hour per map with 1000 realisations):
   ```bash
   for band in K Ka Q1 Q2 030 044 MFI11 MFI13 MFI17H2 MFI17H4 MFI19H2 MFI19H4; do
       python noise_mc.py $band
   done
   ```

3. **Make the weighted maps**:
   ```bash
   python combine.py
   ```

4. **Tables and checks**:
   ```bash
   python sn_table.py              # Table 1
   python spectral_indices.py      # Table 2, T-T and Fuskeland plots
   python check_combination.py     # noise of the weighted maps vs the simulations
   python horn_correlation.py      # noise correlation between frequencies of the same horn
   python pixel_window.py          # effect of the pixel window
   ```

5. **Figures**:
   ```bash
   python make_figures.py
   python plot_regions.py
   python plot_lic.py
   ```

---

## **Output**

All products are written to `products/` and the figures to `figures/`:

| File | Contents |
|---|---|
| `products/<map>/<map>_60arcmin_n<nside>.fits` | smoothed I, Q, U and their errors (`I_ERR`, `Q_ERR`, `U_ERR`, `QU_COV`) |
| `products/<map>/<map>_60arcmin_n64_noise_sims.npy` | the Nside-64 noise realisations |
| `products/combined/<name>_60arcmin_n<nside>.fits` | weighted maps at 10 GHz (`Q`, `U`, `P`, errors, number of input maps) |
| `products/sn_table.json` | Table 1 |
| `products/spectral_indices.json` | Table 2 |

The weighted maps are listed at the top of `combine.py`. The main ones are `quijote` (QUIJOTE only), `wmapplanck` (WMAP and *Planck*) and `full` (all of them).

---

## **Notes**
- The simulations contain white noise only. The QUIJOTE noise includes the DR1 rescaling factors, which account for the 1/f noise at the pixel scale; the extra noise from 1/f correlations at 1° is applied in the `*_1f` maps.
- All the analyses use Nside 64. The Nside 512 and 256 maps are only used for the figures.
- The FDEC filtering is a re-implementation of the QUIJOTE one, applied to the smoothed WMAP and *Planck* maps.

---

Roke Cepeda-Arroita (roke.cepeda@iac.es)
