''' config.py:
Paths and settings for the re-analysis of the loops and spurs paper.
Change the paths below to point at your copies of the input maps.

Version 1.0 [Oct 2026]
Roke Cepeda-Arroita
roke.cepeda@iac.es
'''

import os


# Input Data

DATA_DIR = '/path/to/data'                                     # root of the input maps
CG_DIR = f'{DATA_DIR}/cosmoglobe'                              # Cosmoglobe DR1 frequency maps
MFI_DIR = f'{DATA_DIR}/quijote_mfi'                            # QUIJOTE MFI DR1 maps (horn maps in horns/)
BEAM_DIR = f'{DATA_DIR}/beams'                                 # WMAP9, NPIPE and QUIJOTE RIMO beam transfer functions
MFI_MASK = f'{MFI_DIR}/masks/mask_quijote_satband_nside512.fits'                    # DR1 satellite-band mask
QT_MASK = f'{DATA_DIR}/quijote_masks/mask_quijote_ncp_satband_nside512.fits'        # satellite band + NCP, for the T-T plots
QT_MASK_PLOT = f'{DATA_DIR}/quijote_masks/mask_quijote_satband_ncp86_nside512.fits' # same, used to blank the map figures
REGION_DIR = f'{DATA_DIR}/quijote_masks'                       # region masks of Table 2
HASLAM = f'{DATA_DIR}/ancillary/haslam408_ds_Remazeilles2014.fits'
LOOP_REGION_DIR = f'{DATA_DIR}/regions/vidal2015'              # loop traces of Vidal et al. (2015), for Fig. 3
NEW_REGION_DIR = f'{DATA_DIR}/regions/new'                     # traces of the arcs and shells identified here

# Weighted maps of the previous version of the paper (WMAP9 + NPIPE), only for Appendix A
ORIGINAL_MAPS = {
    'wmapplanck': f'{DATA_DIR}/original/512_60.0smoothed_quijotecombwei10_tqu_v1.5_noise_v1.0_-3.1_wmapplanck_combine.fits',
    'full':       f'{DATA_DIR}/original/512_60.0smoothed_quijotecombwei10_tqu_v1.5_noise_v1.0_-3.1_wmapplanck_qtmask_combine.fits',
}


# Output

HERE = os.path.dirname(os.path.abspath(__file__))
PRODUCTS = f'{HERE}/products'                                  # smoothed maps, noise, weighted maps, tables
FIGURES = f'{HERE}/figures'                                    # paper figures


# Settings

FWHM = 60.0                   # arcmin, final resolution
NSIDES = (512, 256, 64)       # output Nsides (all analyses use 64, the others are for display)
NSIM = 1000                   # noise realisations per map
BETA = -3.1                   # spectral index used to scale all maps to NU0
NU0 = 10.0                    # GHz, reference frequency of the weighted maps
