#!/usr/bin/env python
# -*- coding: utf-8  -*-
#
# Do some plots of the polarisation maps from the jack-knifes
#
# Version history:
#
# 27-Apr-2019  M. Peel       Started

import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
from astrocode.astroutils import *

prefixes = ['mfi_apr2020b_pwv1', 'mfi_apr2020b_allatonce', 'mfi_apr2020b_altsample1', 'mfi_apr2020b_altsample2', 'mfi_apr2020b_daynight1', 'mfi_apr2020b_daynight2', 'mfi_apr2020b_fivesample1', 'mfi_apr2020b_fivesample2', 'mfi_apr2020b_half1', 'mfi_apr2020b_half2', 'mfi_apr2020b_halfring1', 'mfi_apr2020b_halfring2', 'mfi_apr2020b_period1', 'mfi_apr2020b_period2', 'mfi_apr2020b_period5', 'mfi_apr2020b_period6', 'mfi_apr2020b_pwv2', 'mfi_apr2020b_ring1', 'mfi_apr2020b_ring2', 'mfi_apr2020b_sample1', 'mfi_apr2020b_sample2', 'mfi_apr2020b_tbem1', 'mfi_apr2020b_tbem2', 'mfi_apr2020b_twosample1', 'mfi_apr2020b_twosample2']

indirectory =  '/Volumes/Toshiba5TB2/mfi/quijote_202004b_jk/reform/'
outdirectory = '/Volumes/Toshiba5TB2/mfi/quijote_202004b_jk/plot/'


# prefixes = ['mfi','mfi_recalib3', 'mfi_recalib4','mfi_recalib5']

# indirectory =  '/Users/mpeel/Documents/maps/quijote_202004b/reform/'
# outdirectory = '/Users/mpeel/Documents/maps/quijote_202004b/plot/'

basedir = '/Users/mpeel/Documents/maps/'
mfi_mask = hp.read_map(basedir+'quijote_201907/weighted/mfi_commonmask.fits')

for prefix in prefixes:
	filename = '512_60.0smoothed_'+prefix.replace('apr2020b_','')+'2_19.0_512_202004b_jk_mKCMBunits.fits'
	# try:
	inputmap = hp.read_map(indirectory+filename,field=None)
	# except:
	# 	filename = '512_60.00smoothed_'+prefix.replace('apr2020b_','')+'4_19.0_512_202004b_mKCMBunits.fits'
	# 	inputmap = hp.read_map(indirectory+filename,field=None)

	pol = np.sqrt(inputmap[1]**2+inputmap[2]**2)
	pol[mfi_mask==0] = hp.UNSEEN
	hp.mollview(pol,min=0,max=0.5,cmap='jet',title=filename)
	plt.savefig(outdirectory+filename.replace('.fits','.pdf'))
