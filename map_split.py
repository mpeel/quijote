#!/usr/bin/env python
# -*- coding: utf-8  -*-
#
# Do a comparison of M31 from QUIJOTE, Planck and WMAP
# 
# Version history:
#
# 09-May-2019  M. Peel       Started
# 05-Jun-2019  M. Peel       Generalised to cope with multiple runs

import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
from astrocode.fitspectrum.astroutils import *

# quijote_map = ['/Users/mpeel/Documents/maps/quijote_201905/quijote_combinedmaps_select_test_rfi_final_may2019_final_test_rfi_final_nside512_finalmaskfdec.fits','/Users/mpeel/Documents/maps/quijote_201905/quijote_combinedmaps_select_test_rfi_final_half1_may2019_final_test_rfi_final_nside512_finalmaskfdec.fits','/Users/mpeel/Documents/maps/quijote_201905/quijote_combinedmaps_select_test_rfi_final_half2_may2019_final_test_rfi_final_nside512_finalmaskfdec.fits']
# outdir = '/Users/mpeel/Documents/maps/quijote_201905/reform/'
# prefix = ['mfi_may2019', 'mfi_may2019_half1', 'mfi_may2019_half2']
# types = ['map','nhits','weights']

# quijote_map = ['/Users/mpeel/Documents/maps/quijote_201810/quijote_combinedmaps_select_oct2018_test_rfi_alldata_highfk_nside512.fits']
# outdir = '/Users/mpeel/Documents/maps/quijote_201810/reform/'
# prefix = ['mfi_oct2018']
# types = ['map','mapsmth','nhits','weights']

# quijote_map = ['/Users/mpeel/Documents/maps/quijote_201904/NOMINAL30_period2_apr2019_sky2_nomask_perperiod_nside512_packed.fits','/Users/mpeel/Documents/maps/quijote_201904/NOMINAL35_period6_apr2019_sky2_nomask_perperiod_nside512_packed.fits','/Users/mpeel/Documents/maps/quijote_201904/NOMINAL40_period2_apr2019_sky2_nomask_perperiod_nside512_packed.fits','/Users/mpeel/Documents/maps/quijote_201904/NOMINAL40_period5_apr2019_sky2_nomask_perperiod_nside512_packed.fits','/Users/mpeel/Documents/maps/quijote_201904/NOMINAL50_period2_apr2019_sky2_nomask_perperiod_nside512_packed.fits','/Users/mpeel/Documents/maps/quijote_201904/NOMINAL50_period5_apr2019_sky2_nomask_perperiod_nside512_packed.fits','/Users/mpeel/Documents/maps/quijote_201904/NOMINAL50_period6_apr2019_sky2_nomask_perperiod_nside512_packed.fits','/Users/mpeel/Documents/maps/quijote_201904/NOMINAL60_period1_apr2019_sky2_nomask_perperiod_nside512_packed.fits','/Users/mpeel/Documents/maps/quijote_201904/NOMINAL60_period2_apr2019_sky2_nomask_perperiod_nside512_packed.fits','/Users/mpeel/Documents/maps/quijote_201904/NOMINAL60_period5_apr2019_sky2_nomask_perperiod_nside512_packed.fits','/Users/mpeel/Documents/maps/quijote_201904/NOMINAL60_period6_apr2019_sky2_nomask_perperiod_nside512_packed.fits','/Users/mpeel/Documents/maps/quijote_201904/NOMINAL65_period1_apr2019_sky2_nomask_perperiod_nside512_packed.fits','/Users/mpeel/Documents/maps/quijote_201904/NOMINAL65_period2_apr2019_sky2_nomask_perperiod_nside512_packed.fits','/Users/mpeel/Documents/maps/quijote_201904/NOMINAL65_period6_apr2019_sky2_nomask_perperiod_nside512_packed.fits','/Users/mpeel/Documents/maps/quijote_201904/NOMINAL70_period6_apr2019_sky2_nomask_perperiod_nside512_packed.fits']
# prefix = ['mfi_apr2019_30_2','mfi_apr2019_35_6','mfi_apr2019_40_2','mfi_apr2019_40_5','mfi_apr2019_50_2','mfi_apr2019_50_5','mfi_apr2019_50_6','mfi_apr2019_60_1','mfi_apr2019_60_2','mfi_apr2019_60_5','mfi_apr2019_60_6','mfi_apr2019_65_1','mfi_apr2019_65_2','mfi_apr2019_65_6','mfi_apr2019_70_6']
# types = ['pixels','map','weights','nhits']
# outdir = '/Users/mpeel/Documents/maps/quijote_201904/reform/'

# quijote_map = ['/Users/mpeel/Documents/maps/quijote_201803/NOMINAL30_period2_mar2018_test_rfi_packed.fits','/Users/mpeel/Documents/maps/quijote_201803/NOMINAL35_period6_mar2018_test_rfi_packed.fits','/Users/mpeel/Documents/maps/quijote_201803/NOMINAL40_period2_mar2018_test_rfi_packed.fits','/Users/mpeel/Documents/maps/quijote_201803/NOMINAL40_period5_mar2018_test_rfi_packed.fits','/Users/mpeel/Documents/maps/quijote_201803/NOMINAL50_period2_mar2018_test_rfi_packed.fits','/Users/mpeel/Documents/maps/quijote_201803/NOMINAL50_period5_mar2018_test_rfi_packed.fits','/Users/mpeel/Documents/maps/quijote_201803/NOMINAL50_period6_mar2018_test_rfi_packed.fits','/Users/mpeel/Documents/maps/quijote_201803/NOMINAL60_period1_mar2018_test_rfi_packed.fits','/Users/mpeel/Documents/maps/quijote_201803/NOMINAL60_period2_mar2018_test_rfi_packed.fits','/Users/mpeel/Documents/maps/quijote_201803/NOMINAL60_period5_mar2018_test_rfi_packed.fits','/Users/mpeel/Documents/maps/quijote_201803/NOMINAL60_period6_mar2018_test_rfi_packed.fits','/Users/mpeel/Documents/maps/quijote_201803/NOMINAL65_period1_mar2018_test_rfi_packed.fits','/Users/mpeel/Documents/maps/quijote_201803/NOMINAL65_period2_mar2018_test_rfi_packed.fits','/Users/mpeel/Documents/maps/quijote_201803/NOMINAL65_period6_mar2018_test_rfi_packed.fits','/Users/mpeel/Documents/maps/quijote_201803/NOMINAL70_period6_mar2018_test_rfi_packed.fits']
# prefix = ['mfi_mar2018_30_2','mfi_mar2018_35_6','mfi_mar2018_40_2','mfi_mar2018_40_5','mfi_mar2018_50_2','mfi_mar2018_50_5','mfi_mar2018_50_6','mfi_mar2018_60_1','mfi_mar2018_60_2','mfi_mar2018_60_5','mfi_mar2018_60_6','mfi_mar2018_65_1','mfi_mar2018_65_2','mfi_mar2018_65_6','mfi_mar2018_70_6']
# outdir = '/Users/mpeel/Documents/maps/quijote_201803/reform/'


# quijote_map = ['/Users/mpeel/Documents/maps/quijote_2019_sim/quijote_combinedmaps_select_sims_sky_dipole_oofnoise_may2019_sims_sky_dipole_oofnoise_nside512_finalmaskfdec.fits','/Users/mpeel/Documents/maps/quijote_2019_sim/quijote_combinedmaps_select_sims_sky_dipole_whitenoise_may2019_sims_sky_dipole_whitenoise_nside512_finalmaskfdec.fits']
# prefix = ['mfi_sim2019oof','mfi_sim2019white']
# outdir = '/Users/mpeel/Documents/maps/quijote_2019_sim/reform/'
# types = ['map','nhits','weights']
#types = ['map','weights','nhits']


# quijote_map = ['/Users/mpeel/Documents/maps/quijote_201907/quijote_combinedmaps_select_withcov_rfi_jul2019_withcov_rfi_nside512_finalmaskfdec2.fits','/Users/mpeel/Documents/maps/quijote_201907/quijote_combinedmaps_select_withcov_rfi_half1_jul2019_withcov_rfi_nside512_finalmaskfdec2.fits','/Users/mpeel/Documents/maps/quijote_201907/quijote_combinedmaps_select_withcov_rfi_half2_jul2019_withcov_rfi_nside512_finalmaskfdec2.fits']
# outdir = '/Users/mpeel/Documents/maps/quijote_201907/reform/'
# prefix = ['mfi_jul2019', 'mfi_jul2019_half1', 'mfi_jul2019_half2']
# types = ['map','mapsmth','nhits','weights']

# quijote_map = ['/Users/mpeel/Documents/maps/quijote_201911/quijote_combinedmaps_select_withcov_rfi_newflagstat_nov2019_recalib_nside512_finalmaskfdec2.fits','/Users/mpeel/Documents/maps/quijote_201911/quijote_combinedmaps_select_withcov_rfi_newflagstat_half1_nov2019_recalib_nside512_finalmaskfdec2.fits','/Users/mpeel/Documents/maps/quijote_201911/quijote_combinedmaps_select_withcov_rfi_newflagstat_half2_nov2019_recalib_nside512_finalmaskfdec2.fits']
# outdir = '/Users/mpeel/Documents/maps/quijote_201911/reform/'
# prefix = ['mfi_nov2019', 'mfi_nov2019_half1', 'mfi_nov2019_half2']
# types = ['map','mapsmth','nhits','weights']


# quijote_map = ['/Volumes/Toshiba5TB2/mfi/quijote_202003/quijote_combinedmaps_nominalselect_newang2_withcov_rfi_newflagstat_allatonce5s_mar2020_newang2_withcov_rfi_newflagstat_allatonce_nside512_finalmaskfdec2.fits']
# outdir = '/Users/mpeel/Documents/maps/quijote_202003/reform/'
# prefix = ['mfi_mar2020']#, 'mfi_nov2019_half1', 'mfi_nov2019_half2']
# types = ['map','nhits','weights']

# quijote_map = ['/Users/mpeel/Documents/maps/quijote_202004/quijote_combinedmaps_nominalselect_recalib_withcov_rfi_newflagstat_allatonce5s_apr2020_recalib_withcov_rfi_newflagstat_allatonce_nside512_finalmaskfdec2.fits']
# outdir = '/Users/mpeel/Documents/maps/quijote_202004/reform/'
# prefix = ['mfi_apr2020']#, 'mfi_nov2019_half1', 'mfi_nov2019_half2']
# types = ['map','nhits','weights']

# quijote_map = ['/Users/mpeel/Documents/maps/quijote_202004b/quijote_combinedmaps_nominalselect_recalib2_withcov_rfi_newflagstat_allatonce5s_apr2020_recalib2_withcov_rfi_newflagstat_allatonce_nside512_finalmaskfdec2.fits']
# outdir = '/Users/mpeel/Documents/maps/quijote_202004b/reform/'
# prefix = ['mfi_apr2020b']#, 'mfi_nov2019_half1', 'mfi_nov2019_half2']
# types = ['map','nhits','weights']

# quijote_map = ['/Users/mpeel/Documents/maps/quijote_202004b/quijote_combinedmaps_nominalselect_recalib6_withcov_rfi_newflagstat_raw_allatonce5s_apr2020_recalib6_withcov_rfi_newflagstat_raw_allatonce_nside512_finalmaskfdec2.fits']#'/Users/mpeel/Documents/maps/quijote_202004b/quijote_combinedmaps_nominalselect_recalib3_withcov_rfi_newflagstat_raw_allatonce5s_apr2020_recalib3_withcov_rfi_newflagstat_raw_allatonce_nside512_finalmaskfdec2.fits','/Users/mpeel/Documents/maps/quijote_202004b/quijote_combinedmaps_nominalselect_recalib4_withcov_rfi_newflagstat_raw_allatonce5s_apr2020_recalib4_withcov_rfi_newflagstat_raw_allatonce_nside512_finalmaskfdec2.fits','/Users/mpeel/Documents/maps/quijote_202004b/quijote_combinedmaps_nominalselect_recalib5_withcov_rfi_newflagstat_raw_allatonce5s_apr2020_recalib5_withcov_rfi_newflagstat_raw_allatonce_nside512_finalmaskfdec2.fits']
# outdir = '/Users/mpeel/Documents/maps/quijote_202004b/reform/'
# prefix = ['mfi_apr2020_recalib6']#,'mfi_apr2020_recalib3','mfi_apr2020_recalib4','mfi_apr2020_recalib5']#, 'mfi_nov2019_half1', 'mfi_nov2019_half2']
# types = ['map','nhits','weights']

# quijote_map = ['/Users/mpeel/Documents/maps/quijote_202004b/quijote_combinedmaps_nominalselect_recalib7_withcov_rfi_newflagstat_raw_allatonce5s_apr2020_recalib7_withcov_rfi_newflagstat_raw_allatonce_nside512_finalmaskfdec2.fits','/Users/mpeel/Documents/maps/quijote_202004b/quijote_combinedmaps_nominalselect_recalib8_withcov_rfi_newflagstat_raw_allatonce5s_apr2020_recalib8_withcov_rfi_newflagstat_raw_allatonce_nside512_finalmaskfdec2.fits']
# prefix = ['mfi_apr2020_recalib7','mfi_apr2020_recalib8']

# quijote_map = ['/Users/mpeel/Documents/maps/quijote_202004b/quijote_combinedmaps_nominalselect_recalib9_withcov_rfi_newflagstat_raw_allatonce5s_apr2020_recalib9_withcov_rfi_newflagstat_raw_allatonce_nside512_finalmaskfdec2.fits']
# prefix = ['mfi_apr2020_recalib9']

# Jack-knifes
# base = '/Volumes/Toshiba5TB2/mfi/quijote_202004b_jk/'
# quijote_map = [base+'quijote_combinedmaps_nominalselect_recalib2_withcov_rfi_newflagstat_raw_allatonce5s_apr2020_recalib2_withcov_rfi_newflagstat_raw_pwv1_nside512_finalmaskfdec2.fits',base+'quijote_combinedmaps_nominalselect_recalib2_withcov_rfi_newflagstat_raw_allatonce5s_apr2020_recalib2_withcov_rfi_newflagstat_raw_allatonce_nside512_finalmaskfdec2.fits',base+'quijote_combinedmaps_nominalselect_recalib2_withcov_rfi_newflagstat_raw_allatonce5s_apr2020_recalib2_withcov_rfi_newflagstat_raw_altsample1_nside512_finalmaskfdec2.fits',base+'quijote_combinedmaps_nominalselect_recalib2_withcov_rfi_newflagstat_raw_allatonce5s_apr2020_recalib2_withcov_rfi_newflagstat_raw_altsample2_nside512_finalmaskfdec2.fits',base+'quijote_combinedmaps_nominalselect_recalib2_withcov_rfi_newflagstat_raw_allatonce5s_apr2020_recalib2_withcov_rfi_newflagstat_raw_daynight1_nside512_finalmaskfdec2.fits',base+'quijote_combinedmaps_nominalselect_recalib2_withcov_rfi_newflagstat_raw_allatonce5s_apr2020_recalib2_withcov_rfi_newflagstat_raw_daynight2_nside512_finalmaskfdec2.fits',base+'quijote_combinedmaps_nominalselect_recalib2_withcov_rfi_newflagstat_raw_allatonce5s_apr2020_recalib2_withcov_rfi_newflagstat_raw_fivesample1_nside512_finalmaskfdec2.fits',base+'quijote_combinedmaps_nominalselect_recalib2_withcov_rfi_newflagstat_raw_allatonce5s_apr2020_recalib2_withcov_rfi_newflagstat_raw_fivesample2_nside512_finalmaskfdec2.fits',base+'quijote_combinedmaps_nominalselect_recalib2_withcov_rfi_newflagstat_raw_allatonce5s_apr2020_recalib2_withcov_rfi_newflagstat_raw_half1_nside512_finalmaskfdec2.fits',base+'quijote_combinedmaps_nominalselect_recalib2_withcov_rfi_newflagstat_raw_allatonce5s_apr2020_recalib2_withcov_rfi_newflagstat_raw_half2_nside512_finalmaskfdec2.fits',base+'quijote_combinedmaps_nominalselect_recalib2_withcov_rfi_newflagstat_raw_allatonce5s_apr2020_recalib2_withcov_rfi_newflagstat_raw_halfring1_nside512_finalmaskfdec2.fits',base+'quijote_combinedmaps_nominalselect_recalib2_withcov_rfi_newflagstat_raw_allatonce5s_apr2020_recalib2_withcov_rfi_newflagstat_raw_halfring2_nside512_finalmaskfdec2.fits',base+'quijote_combinedmaps_nominalselect_recalib2_withcov_rfi_newflagstat_raw_allatonce5s_apr2020_recalib2_withcov_rfi_newflagstat_raw_period1_nside512_finalmaskfdec2.fits',base+'quijote_combinedmaps_nominalselect_recalib2_withcov_rfi_newflagstat_raw_allatonce5s_apr2020_recalib2_withcov_rfi_newflagstat_raw_period2_nside512_finalmaskfdec2.fits',base+'quijote_combinedmaps_nominalselect_recalib2_withcov_rfi_newflagstat_raw_allatonce5s_apr2020_recalib2_withcov_rfi_newflagstat_raw_period5_nside512_finalmaskfdec2.fits',base+'quijote_combinedmaps_nominalselect_recalib2_withcov_rfi_newflagstat_raw_allatonce5s_apr2020_recalib2_withcov_rfi_newflagstat_raw_period6_nside512_finalmaskfdec2.fits',base+'quijote_combinedmaps_nominalselect_recalib2_withcov_rfi_newflagstat_raw_allatonce5s_apr2020_recalib2_withcov_rfi_newflagstat_raw_pwv2_nside512_finalmaskfdec2.fits',base+'quijote_combinedmaps_nominalselect_recalib2_withcov_rfi_newflagstat_raw_allatonce5s_apr2020_recalib2_withcov_rfi_newflagstat_raw_ring1_nside512_finalmaskfdec2.fits',base+'quijote_combinedmaps_nominalselect_recalib2_withcov_rfi_newflagstat_raw_allatonce5s_apr2020_recalib2_withcov_rfi_newflagstat_raw_ring2_nside512_finalmaskfdec2.fits',base+'quijote_combinedmaps_nominalselect_recalib2_withcov_rfi_newflagstat_raw_allatonce5s_apr2020_recalib2_withcov_rfi_newflagstat_raw_sample1_nside512_finalmaskfdec2.fits',base+'quijote_combinedmaps_nominalselect_recalib2_withcov_rfi_newflagstat_raw_allatonce5s_apr2020_recalib2_withcov_rfi_newflagstat_raw_sample2_nside512_finalmaskfdec2.fits',base+'quijote_combinedmaps_nominalselect_recalib2_withcov_rfi_newflagstat_raw_allatonce5s_apr2020_recalib2_withcov_rfi_newflagstat_raw_tbem1_nside512_finalmaskfdec2.fits',base+'quijote_combinedmaps_nominalselect_recalib2_withcov_rfi_newflagstat_raw_allatonce5s_apr2020_recalib2_withcov_rfi_newflagstat_raw_tbem2_nside512_finalmaskfdec2.fits',base+'quijote_combinedmaps_nominalselect_recalib2_withcov_rfi_newflagstat_raw_allatonce5s_apr2020_recalib2_withcov_rfi_newflagstat_raw_twosample1_nside512_finalmaskfdec2.fits',base+'quijote_combinedmaps_nominalselect_recalib2_withcov_rfi_newflagstat_raw_allatonce5s_apr2020_recalib2_withcov_rfi_newflagstat_raw_twosample2_nside512_finalmaskfdec2.fits']
# outdir = base+'reform/'
# prefix = ['mfi_apr2020b_pwv1', 'mfi_apr2020b_allatonce', 'mfi_apr2020b_altsample1', 'mfi_apr2020b_altsample2', 'mfi_apr2020b_daynight1', 'mfi_apr2020b_daynight2', 'mfi_apr2020b_fivesample1', 'mfi_apr2020b_fivesample2', 'mfi_apr2020b_half1', 'mfi_apr2020b_half2', 'mfi_apr2020b_halfring1', 'mfi_apr2020b_halfring2', 'mfi_apr2020b_period1', 'mfi_apr2020b_period2', 'mfi_apr2020b_period5', 'mfi_apr2020b_period6', 'mfi_apr2020b_pwv2', 'mfi_apr2020b_ring1', 'mfi_apr2020b_ring2', 'mfi_apr2020b_sample1', 'mfi_apr2020b_sample2', 'mfi_apr2020b_tbem1', 'mfi_apr2020b_tbem2', 'mfi_apr2020b_twosample1', 'mfi_apr2020b_twosample2']
# # prefix = ['mfi_apr2020b']#, 'mfi_nov2019_half1', 'mfi_nov2019_half2']
# types = ['map','nhits','weights']

# quijote_map = ['/Users/mpeel/Documents/maps/quijote_202007_sim/quijote_combinedmaps_nominalselect_sims_sky_dipole_whitenoise_withcov_rfi_newflagstat_allatonce5s_july2020_allatonce_nside512_finalmaskfdec2.fits']
# prefix = ['mfisim_jul2020_sim']
# outdir = '/Users/mpeel/Documents/maps/quijote_202007_sim/'

# quijote_map = ['/Users/mpeel/Documents/maps/quijote_202007_sim/quijote_combinedmaps_nominalselect_sims_sky_dipole_oofnoise_withcov_rfi_newflagstat_allatonce5s_july2020_allatonce_nside512_finalmaskfdec2.fits']
# prefix = ['mfisim_jul2020_simoof']
# outdir = '/Users/mpeel/Documents/maps/quijote_202007_sim/'


invert_weights = True
indir = '/Users/mpeel/Documents/maps/quijote_202011/'
outdir = indir+'reform/'
types = ['map','mapsmth','nhits','weights']
quijote_map = ['quijote_combinedmaps_nominalselect_recalib11_withcov_rfi_newflagstat_raw_allatonce5s_apr2020_recalib11_withcov_rfi_newflagstat_raw_allatonce_nside512_finalmaskfdec2.fits','quijote_combinedmaps_nominalselect_recalib11_withcov_rfi_newflagstat_raw_allatonce5s_apr2020_recalib11_withcov_rfi_newflagstat_raw_half1_nside512_finalmaskfdec2.fits','quijote_combinedmaps_nominalselect_recalib11_withcov_rfi_newflagstat_raw_allatonce5s_apr2020_recalib11_withcov_rfi_newflagstat_raw_half2_nside512_finalmaskfdec2.fits','quijote_combinedmaps_nominalselect_recalib11_withcov_rfi_newflagstat_raw_allatonce5s_apr2020_recalib11_withcov_rfi_newflagstat_raw_halfring1_nside512_finalmaskfdec2.fits','quijote_combinedmaps_nominalselect_recalib11_withcov_rfi_newflagstat_raw_allatonce5s_apr2020_recalib11_withcov_rfi_newflagstat_raw_halfring2_nside512_finalmaskfdec2.fits','quijote_combinedmaps_nominalselect_recalib11_withcov_rfi_newflagstat_raw_allatonce5s_apr2020_recalib11_withcov_rfi_newflagstat_raw_period1_nside512_finalmaskfdec2.fits','quijote_combinedmaps_nominalselect_recalib11_withcov_rfi_newflagstat_raw_allatonce5s_apr2020_recalib11_withcov_rfi_newflagstat_raw_period2_nside512_finalmaskfdec2.fits','quijote_combinedmaps_nominalselect_recalib11_withcov_rfi_newflagstat_raw_allatonce5s_apr2020_recalib11_withcov_rfi_newflagstat_raw_period5_nside512_finalmaskfdec2.fits','quijote_combinedmaps_nominalselect_recalib11_withcov_rfi_newflagstat_raw_allatonce5s_apr2020_recalib11_withcov_rfi_newflagstat_raw_period6_nside512_finalmaskfdec2.fits','quijote_combinedmaps_nominalselect_recalib11_withcov_rfi_newflagstat_raw_allatonce5s_apr2020_recalib11_withcov_rfi_newflagstat_raw_pwv1_nside512_finalmaskfdec2.fits','quijote_combinedmaps_nominalselect_recalib11_withcov_rfi_newflagstat_raw_allatonce5s_apr2020_recalib11_withcov_rfi_newflagstat_raw_pwv2_nside512_finalmaskfdec2.fits']
prefix = ['mfi_nov2020','mfi_nov2020_half1','mfi_nov2020_half2','mfi_nov2020_halfring1','mfi_nov2020_halfring2','mfi_nov2020_period1','mfi_nov2020_period2','mfi_nov2020_period5','mfi_nov2020_period6','mfi_nov2020_pwv1','mfi_nov2020_pwv2']

invert_weights = True
indir = '/Users/mpeel/Documents/maps/quijote_m31/'
outdir = indir+'reform/'
types = ['map','nhits','weights']
quijote_map = ['quijote_combinedmaps_nominal_m31_recalib11_withcov_rfi_newflagstat_raw_allatonce5s_jan2021_allatonce_nside512_finalmaskfdec2.fits','quijote_combinedmaps_nominal_m31_recalib11_withcov_rfi_newflagstat_raw_allatonce5s_jan2021_halfring1_nside512_finalmaskfdec2.fits','quijote_combinedmaps_nominal_m31_recalib11_withcov_rfi_newflagstat_raw_allatonce5s_jan2021_halfring2_nside512_finalmaskfdec2.fits']
prefix = ['mfi_m31','mfi_m31half1','mfi_m31half2']
scale = [[1.0, 1.0, 1.0, 1.0],[1.0, 1.0, 1.0, 1.0]]
# scale = [[1.0, 0.91, 0.88, 0.96], [1.0, 0.88, 0.88, 0.92]]

print(len(quijote_map))
print(len(prefix))
for mapnum in range(0,len(quijote_map)):
	print(quijote_map[mapnum], prefix[mapnum])
	# continue

	inputfits = fits.open(indir+quijote_map[mapnum])
	print(inputfits[1].header)

	col_names = inputfits[1].columns.names
	print(col_names)
	nmaps = len(col_names)
	maps = []
	# exit()
	# print(inputfits[1].data)
	# print(inputfits[1].data[0][0])

	nside = 256
	npix = 12*nside**2
	nummaps = 8
	pol = 3
	maps = np.zeros((nummaps,pol,npix))
	data = inputfits[1].data[0]

	# 1 = 11/13GHz, broken
	# 2 = 17/19GHz, working
	# 3 = 11/13GHz, working
	# 4 = 17/19GHz, working
	freqs = [[11.0,17.0,11.0,17.0],[13.0,19.0,13.0,19.0]]

	# data[1] = pixel numbers
	# data[2][0] = I
	# data[2][1] = Q
	# data[2][2] = U

	if types[0] == 'pixels':
		# if data[5] != nside:
		nside = data[5]
		print(nside)
		maps = np.zeros((pol,npix))
		maps[0][:] = hp.pixelfunc.UNSEEN
		maps[1][:] = hp.pixelfunc.UNSEEN
		maps[2][:] = hp.pixelfunc.UNSEEN
		# print(len(data[0]))
		# print(len(data[0][0]))
		# print(data[0][0][0][0][0])
		print(np.shape(data[1][0][0][0]))
		count_types = len(types)
		pixelnumbers = data[0]
		print(pixelnumbers)

		for i in range(1,count_types): # map/nhits/weights
			for k in range(0,2): # Low/high band
				for l in range(0,4): # Horns
					print(scale[k][l])
					maps[0][pixelnumbers[:]] = data[i][0][0][k][l][:] / scale[k][l]
					maps[1][pixelnumbers[:]] = data[i][0][1][k][l][:]
					maps[2][pixelnumbers[:]] = data[i][0][2][k][l][:]

					outname = prefix[mapnum]+'_'+types[i]+'_'+str(freqs[k][l])+'_'+str(l+1)+'.fits'
					hp.write_map(outdir+outname,maps,overwrite=True)
					hp.mollview(maps[0],norm='hist')
					plt.savefig(outdir+outname+'_I.pdf')
					plt.clf()
					hp.mollview(maps[1],norm='hist')
					plt.savefig(outdir+outname+'_Q.pdf')
					plt.clf()
					hp.mollview(maps[2],norm='hist')
					plt.savefig(outdir+outname+'_U.pdf')
					plt.clf()
	else:
		print(len(data[0]))
		print(len(data[0][0]))
		print(data[0][0][0][0][0])
		print(np.shape(data[0][0][0][0]))
		count_types = len(types)
		for i in range(0,count_types): # map/nhits/weights
			for k in range(0,2): # Low/high band
				for l in range(0,4): # Horns
					outname = prefix[mapnum]+'_'+types[i]+'_'+str(freqs[k][l])+'_'+str(l+1)+'.fits'
					print(outname)
					if types[i] in ['map','mapsmth']:
						print('Rescaling by ' + str(scale[k][l]))
						data[i][0][k][l][:] /= scale[k][l]
						data[i][1][k][l][:] /= scale[k][l]
						data[i][2][k][l][:] /= scale[k][l]
					if types[i] == 'weights':
						print('Rescaling by ' + str((1.0/scale[k][l])**2))
						data[i][0][k][l][:] *= scale[k][l]**2
						data[i][1][k][l][:] *= scale[k][l]**2
						data[i][2][k][l][:] *= scale[k][l]**2
					if types[i] == 'weights' and invert_weights:
						# Invert the weights map
						data[i][0][k][l][data[1][0][k][l]!=0]	= 1.0 / data[i][0][k][l][data[1][0][k][l]!=0]
						data[i][1][k][l][data[1][0][k][l]!=0]	= 1.0 / data[i][1][k][l][data[1][0][k][l]!=0]
						data[i][2][k][l][data[1][0][k][l]!=0]	= 1.0 / data[i][2][k][l][data[1][0][k][l]!=0]

					data[i][0][k][l][data[1][0][k][l]==0] = hp.pixelfunc.UNSEEN
					data[i][1][k][l][data[1][0][k][l]==0] = hp.pixelfunc.UNSEEN
					data[i][2][k][l][data[1][0][k][l]==0] = hp.pixelfunc.UNSEEN
					outmap = [data[i][0][k][l], data[i][1][k][l], data[i][2][k][l]]
					print(np.shape(outmap))
					hp.write_map(outdir+outname,outmap,overwrite=True)
					hp.mollview(outmap[0],norm='hist')
					plt.savefig(outdir+outname+'_I.pdf')
					plt.clf()
					hp.mollview(outmap[1],norm='hist')
					plt.savefig(outdir+outname+'_Q.pdf')
					plt.clf()
					hp.mollview(outmap[2],norm='hist')
					plt.savefig(outdir+outname+'_U.pdf')
					plt.clf()

# EOF