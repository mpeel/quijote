#!/usr/bin/env python
# -*- coding: utf-8  -*-
#
# Do a quick comparison of the simulated standard deviations with different nsides
# 
# Version history:
#
# 03-Feb-2020  M. Peel       Started

import healpy as hp
import numpy as np
import matplotlib.pyplot as plt

# nside = [512, 256, 64]
# npix = hp.nside2npix(nside)

indirectory = '/Volumes/Toshiba5TB2/mfi/validation/Nov2019/noise_simulations/'
outdirectory = '/Volumes/Toshiba5TB2/mfi/validation/Nov2019/noise_simulations/analyse/17_2'
mask = '/Users/mpeel/Documents/maps/quijote_masks/mask_quijote_ncp_lowdec_nside512.fits'

map_512 = hp.read_map(indirectory+'sims_quijote_17GHz_horn2_new/std_H2_17_sm1deg_nside512.fits',field=None)
map_256 = hp.read_map(indirectory+'sims_quijote_17GHz_horn2_new/std_H2_17_sm1deg_nside256.fits',field=None)
map_64 = hp.read_map(indirectory+'sims_quijote_17GHz_horn2_new/std_H2_17_sm1deg_nside64.fits',field=None)
# map_512 = hp.read_map(indirectory+'sims_quijote_11GHz_horn3_new/std_H3_11_sm1deg_nside512.fits',field=None)
# map_256 = hp.read_map(indirectory+'sims_quijote_11GHz_horn3_new/std_H3_11_sm1deg_nside256.fits',field=None)
# map_64 = hp.read_map(indirectory+'sims_quijote_11GHz_horn3_new/std_H3_11_sm1deg_nside64.fits',field=None)

mask_512 = hp.read_map(mask)
mask_256 = hp.ud_grade(mask_512,nside_out=256)
mask_64 = hp.ud_grade(mask_512,nside_out=64)

power=0
map_512_64 = np.sqrt(hp.ud_grade(map_512**2, nside_out=64,power=power))
map_512_256 = np.sqrt(hp.ud_grade(map_512**2, nside_out=256,power=power))
map_256_64 = np.sqrt(hp.ud_grade(map_256**2, nside_out=64,power=power))

for i in range(0,3):
	hp.mollview(map_512[i])
	plt.savefig(outdirectory+'std_512_'+str(i)+'.pdf')
	hp.mollview(map_256[i])
	plt.savefig(outdirectory+'std_256_'+str(i)+'.pdf')
	hp.mollview(map_64[i])
	plt.savefig(outdirectory+'std_64_'+str(i)+'.pdf')

	hp.mollview(map_512_64[i])
	plt.savefig(outdirectory+'std_512_64_'+str(i)+'.pdf')

	hp.mollview(map_512_256[i])
	plt.savefig(outdirectory+'std_512_256_'+str(i)+'.pdf')

	hp.mollview(map_256_64[i])
	plt.savefig(outdirectory+'std_256_64_'+str(i)+'.pdf')

	print('std')
	print(np.std(map_512[i][mask_512==1]))
	print(np.std(map_256[i][mask_256==1]))
	print(np.std(map_64[i][mask_64==1]))

	print('median')
	print(np.median(map_512[i][mask_512==1]))
	print(np.median(map_256[i][mask_256==1]))
	print(np.median(map_64[i][mask_64==1]))

	print('512_to_64')
	print(np.std(map_512_64[i][mask_64==1]))
	print(np.median(map_512_64[i][mask_64==1]))
	print('512_to_256')
	print(np.std(map_512_256[i][mask_256==1]))
	print(np.median(map_512_256[i][mask_256==1]))

	print('256_to_64')
	print(np.std(map_256_64[i][mask_64==1]))
	print(np.median(map_256_64[i][mask_64==1]))
