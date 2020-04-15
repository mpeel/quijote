#!/usr/bin/env python
# -*- coding: utf-8  -*-
#
# Do a cross-check of the polarisation angle statistics
# 
# Version history:
#
# 14-Apr-2020  M. Peel       Started

import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from scipy import odr
from astrocode.fitspectrum.astroutils import *
from compare_polang import *
from astrocode.fitspectrum.spectra import *
def dodiff(map1, map2):
	diff = map1-map2
	diff[diff < -90] += 180.0
	diff[diff > 90] -= 180.0
	return diff

indirectory = '/Volumes/Toshiba5TB2/mfi/validation/simulations_vielva/'
outdirectory = indirectory+'_analysis/'
mask_file = '/Users/mpeel/Documents/maps/quijote_202004/analyse/_quickmask.fits'
nside_out = 64

mask = hp.read_map(mask_file)

vals = []
for i in range(1,101):
	filename = 'quijote_11GHz_horn3_'+'{:04d}'.format(i)+'_sm1deg.fits'
	test = hp.read_map(indirectory+filename,field=None)
	qmap = hp.ud_grade(test[1],nside_out,order_in='RING',order_out='RING')
	umap = hp.ud_grade(test[2],nside_out,order_in='RING',order_out='RING')

	polang_map = calc_polang(qmap,umap)
	val = np.median(polang_map[mask==1])
	vals.append(val)
print(vals)
n, bins, patches = plt.hist(vals, 50, density=True, facecolor='g', alpha=0.75)
plt.savefig('_test.pdf')