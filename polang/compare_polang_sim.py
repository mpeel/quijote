#!/usr/bin/env python
# -*- coding: utf-8  -*-
#
# Compare the polarisation angles from simulated maps
#
# Version history:
#
# 13-Jul-2020  M. Peel       Started

import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from scipy import odr
from astrocode.astroutils import *
from compare_polang import *
from astrocode.spectra import *

def dodiff(map1, map2):
	diff = map1-map2
	diff[diff < -90] += 180.0
	diff[diff > 90] -= 180.0
	return diff

indirectory = '/Users/mpeel/Documents/maps/quijote_202007_sim/'
outdirectory = indirectory+'_analysis/'
mask_file = indirectory+'_quickmask.fits'
nside_out = 64
dosubtract = True

mask = hp.read_map(mask_file)

data_filenames = [indirectory+'mfisim_jul2020_simoof_map_11.0_3.fits',indirectory+'mfisim_jul2020_simoof_map_13.0_3.fits',indirectory+'mfisim_jul2020_simoof_map_17.0_2.fits',indirectory+'mfisim_jul2020_simoof_map_17.0_4.fits',indirectory+'mfisim_jul2020_simoof_map_19.0_2.fits',indirectory+'mfisim_jul2020_simoof_map_19.0_4.fits']
sim_filenames = [indirectory+'total_11GHz_0512.fits',indirectory+'total_13GHz_0512.fits',indirectory+'total_17GHz_0512.fits',indirectory+'total_17GHz_0512.fits',indirectory+'total_19GHz_0512.fits',indirectory+'total_19GHz_0512.fits']
var_filename = [indirectory+'mfisim_jul2020_simoof_weights_11.0_3.fits',indirectory+'mfisim_jul2020_simoof_weights_13.0_3.fits',indirectory+'mfisim_jul2020_simoof_weights_17.0_2.fits',indirectory+'mfisim_jul2020_simoof_weights_17.0_4.fits',indirectory+'mfisim_jul2020_simoof_weights_19.0_2.fits',indirectory+'mfisim_jul2020_simoof_weights_19.0_4.fits']

for j in range(0,len(data_filenames)):
	print(data_filenames[j])
	data = hp.read_map(data_filenames[j],field=None)
	data_qmap = hp.ud_grade(data[1],nside_out,order_in='RING',order_out='RING')
	data_umap = hp.ud_grade(data[2],nside_out,order_in='RING',order_out='RING')

	sim = hp.read_map(sim_filenames[j],field=None)
	sim_qmap = hp.ud_grade(sim[1],nside_out,order_in='RING',order_out='RING')
	sim_umap = hp.ud_grade(sim[2],nside_out,order_in='RING',order_out='RING')

	var = hp.read_map(var_filename[j],field=None)
	var_q = hp.ud_grade(1.0/np.sqrt(var[1]),nside_out,order_in='RING',order_out='RING',power=2)
	var_u = hp.ud_grade(1.0/np.sqrt(var[2]),nside_out,order_in='RING',order_out='RING',power=2)

	hp.mollview(data_qmap)
	plt.savefig(outdirectory+'qmap_'+str(j)+'.png')
	plt.clf()
	hp.mollview(data_umap)
	plt.savefig(outdirectory+'umap_'+str(j)+'.png')
	plt.clf()

	polang_map = calc_polang(data_qmap,data_umap)
	polang_unc_map = calc_polang_unc(data_qmap,data_umap,var_q,var_u)
	hp.mollview(polang_map)
	plt.savefig(outdirectory+'polang_'+str(j)+'.png')
	plt.clf()

	polang_sim = calc_polang(sim_qmap,sim_umap)
	polang_unc_sim = calc_polang_unc(sim_qmap,sim_umap,var_q,var_u)
	hp.mollview(polang_sim)
	plt.savefig(outdirectory+'polang_sim_'+str(j)+'.png')
	plt.clf()

	val = np.median(polang_map[mask==1]-polang_sim[mask==1])
	std = calc_std_over_n(polang_map[mask==1]-polang_sim[mask==1])
	print([val, std])
	val = np.average(polang_map[mask==1]-polang_sim[mask==1],weights=1.0/(polang_unc_sim[mask==1]**2.0),returned=True)
	print(val)
