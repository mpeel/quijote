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
indirectory2 = '/Users/mpeel/Documents/maps/quijote_202004/'
outdirectory = indirectory+'_analysis/'
mask_file = indirectory2+'analyse/_quickmask.fits'
nside_out = 64
dosubtract = True

mask = hp.read_map(mask_file)

data_filenames = [indirectory2+'reform/512_60.00smoothed_mfi3_11.0_512_202004_mKCMBunits.fits',indirectory2+'reform/512_60.00smoothed_mfi3_13.0_512_202004_mKCMBunits.fits',indirectory2+'reform/512_60.0smoothed_mfi2_17.0_512_202004_mKCMBunits.fits',indirectory2+'reform/512_60.0smoothed_mfi2_19.0_512_202004_mKCMBunits.fits',indirectory2+'reform/512_60.0smoothed_mfi4_17.0_512_202004_mKCMBunits.fits',indirectory2+'reform/512_60.0smoothed_mfi4_19.0_512_202004_mKCMBunits.fits']
sim_filenames_prefix = [indirectory+'quijote_11GHz_horn3_',indirectory+'quijote_13GHz_horn3_',indirectory+'quijote_17GHz_horn2_',indirectory+'quijote_19GHz_horn2_',indirectory+'quijote_17GHz_horn4_',indirectory+'quijote_19GHz_horn4_']
sim_postfix = '_sm1deg.fits'
var_filename = [indirectory2 + '../quijote_201911/noise/std_H3_11_sm1deg_nside64.fits', indirectory2 + '../quijote_201911/noise/std_H3_13_sm1deg_nside64.fits', indirectory2 + '../quijote_201911/noise/std_H2_17_sm1deg_nside64.fits', indirectory2 + '../quijote_201911/noise/std_H2_19_sm1deg_nside64.fits', indirectory2 + '../quijote_201911/noise/std_H4_17_sm1deg_nside64.fits', indirectory2 + '../quijote_201911/noise/std_H4_19_sm1deg_nside64.fits']

if dosubtract:
	mapdata = hp.read_map('/Users/mpeel/Documents/maps/wmap9/wmap_band_smth_iqumap_r9_9yr_K_v5.fits',field=None)
	qmap_wmap = hp.ud_grade(mapdata[1],nside_out,order_in='RING',order_out='RING')
	umap_wmap = hp.ud_grade(mapdata[2],nside_out,order_in='RING',order_out='RING')
	hp.mollview(qmap_wmap)
	plt.savefig(outdirectory+'_wmap_q.png')
	plt.clf()
	hp.mollview(umap_wmap)
	plt.savefig(outdirectory+'_wmap_u.png')
	plt.clf()
	polang_wmap = calc_polang(qmap_wmap,umap_wmap)
	hp.mollview(polang_wmap)
	plt.savefig(outdirectory+'_wmap_polang.png')
	plt.clf()

for j in range(0,len(data_filenames)):
	data = hp.read_map(data_filenames[j],field=None)
	data_qmap = hp.ud_grade(data[1],nside_out,order_in='RING',order_out='RING')
	data_umap = hp.ud_grade(data[2],nside_out,order_in='RING',order_out='RING')

	var = hp.read_map(var_filename[j],field=None)
	var_q = var[1].copy()#**2
	var_u = var[2].copy()#**2
	var = []

	hp.mollview(data_qmap)
	plt.savefig(outdirectory+'qmap.png')
	plt.clf()
	hp.mollview(data_umap)
	plt.savefig(outdirectory+'umap.png')
	plt.clf()

	vals = []
	vals2 = []
	vals3 = []
	for i in range(1,101):
		filename = indirectory+'quijote_11GHz_horn3_'+'{:04d}'.format(i)+'_sm1deg.fits'
		test = hp.read_map(sim_filenames_prefix[j]+'{:04d}'.format(i)+sim_postfix,field=None)
		qmap = hp.ud_grade(test[1],nside_out,order_in='RING',order_out='RING')
		umap = hp.ud_grade(test[2],nside_out,order_in='RING',order_out='RING')
		qmap[:] += data_qmap[:]
		umap[:] += data_umap[:]

		hp.mollview(qmap)
		plt.savefig(outdirectory+'qmap_'+str(j)+'_'+str(i)+'.png')
		plt.clf()
		hp.mollview(umap)
		plt.savefig(outdirectory+'umap_'+str(j)+'_'+str(i)+'.png')
		plt.clf()

		polang_map = calc_polang(qmap,umap)
		polang_unc_map = calc_polang_unc(qmap,umap,var_q,var_u)
		hp.mollview(polang_map)
		plt.savefig(outdirectory+'polang_'+str(j)+'_'+str(i)+'.png')
		plt.clf()

		val = np.median(polang_map[mask==1]-polang_wmap[mask==1])
		vals.append(val)

		val = np.average(polang_map[mask==1]-polang_wmap[mask==1],weights=1.0/(polang_unc_map[mask==1]**2.0),returned=True)
		vals2.append(val[0])
		vals3.append(np.sqrt(1.0/val[1]))

	# print(vals)
	print(np.std(vals))
	n, bins, patches = plt.hist(vals, 50, density=True, facecolor='g', alpha=0.75)
	plt.title('std = ' + str(np.std(vals)))
	plt.savefig(outdirectory+'_histogram_median_'+str(j)+'.pdf')
	plt.clf()

	print(np.std(vals2))
	n, bins, patches = plt.hist(vals2, 50, density=True, facecolor='g', alpha=0.75)
	plt.title('std = ' + str(np.std(vals2)))
	plt.savefig(outdirectory+'_histogram_wgt_'+str(j)+'.pdf')
	plt.clf()

	print(np.std(vals3))
	n, bins, patches = plt.hist(vals3, 50, density=True, facecolor='g', alpha=0.75)
	plt.title('std = ' + str(np.std(vals3)))
	plt.savefig(outdirectory+'_histogram_wgterr_'+str(j)+'.pdf')
	plt.clf()
