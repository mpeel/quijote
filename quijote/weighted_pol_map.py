#!/usr/bin/env python
# -*- coding: utf-8  -*-
#
# Do a quick analysis of the MFI maps
# 
# Version history:
#
# 31-May-2019  M. Peel       Started
# 05-Jun-2019  M. Peel       Generalised to cope with multiple runs

import healpy as hp
import numpy as np
import matplotlib.pyplot as plt

def noiserealisation(inputmap, numpixels):
    newmap = np.zeros(numpixels)
    newmap = np.random.normal(scale=1.0, size=numpixels) * inputmap
    return newmap

nside = 512
npix = hp.nside2npix(nside)

indirectory = '/Users/mpeel/Documents/maps/quijote_201905/smooth/'
outdirectory = '/Users/mpeel/Documents/maps/quijote_201905/analyse/'
date='201905'

# indirectory = '/Users/mpeel/Documents/maps/quijote_201810/smooth/'
# outdirectory = '/Users/mpeel/Documents/maps/quijote_201810/analyse/'
# prefix='mfi'
# date='201810'

prefix='half1mfi'
maps_half1 = [str(nside)+'_60.00smoothed_'+prefix+'1_11.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'1_13.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'2_17.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'2_19.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'3_11.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'3_13.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'4_17.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'4_19.0_512_'+date+'_mKCMBunits.fits']

prefix='half2mfi'
maps_half2 = [str(nside)+'_60.00smoothed_'+prefix+'1_11.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'1_13.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'2_17.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'2_19.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'3_11.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'3_13.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'4_17.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'4_19.0_512_'+date+'_mKCMBunits.fits']

prefix='mfi'
maps = [str(nside)+'_60.00smoothed_'+prefix+'1_11.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'1_13.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'2_17.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'2_19.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'3_11.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'3_13.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'4_17.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'4_19.0_512_'+date+'_mKCMBunits.fits']

use_halfrings = False
use_weights = False
use_reweight_by_rms = True

nummaps = len(maps)
freqs = [11,13,17,19,11,13,17,19]
normfreq = 10.0
index = -3.0
commonmask = np.ones(npix)
combine_q = np.zeros(npix)
combine_u = np.zeros(npix)
weight_q = np.zeros(npix)
weight_u = np.zeros(npix)
rescale_vals = np.zeros((3,nummaps))
for i in range(2,nummaps):
	print(maps[i])
	mapdata = hp.read_map(indirectory+maps[i],field=None)
	commonmask[mapdata[0][:] == hp.UNSEEN] = 0
	mapdata[0][mapdata[0][:] == hp.UNSEEN] = 0.0
	mapdata[1][mapdata[1][:] == hp.UNSEEN] = 0.0
	mapdata[2][mapdata[2][:] == hp.UNSEEN] = 0.0
	hp.mollview(mapdata[0],norm='hist')
	plt.savefig(outdirectory+maps[i]+'_0.pdf')
	hp.mollview(mapdata[1],norm='hist')
	plt.savefig(outdirectory+maps[i]+'_1.pdf')
	hp.mollview(mapdata[2],norm='hist')
	plt.savefig(outdirectory+maps[i]+'_2.pdf')
	hp.mollview(np.sqrt(mapdata[1]**2+mapdata[2]**2),min=0,max=1)
	plt.savefig(outdirectory+maps[i]+'_P.pdf')

	if use_halfrings:
		map_half1 = hp.read_map(indirectory+maps_half1[i],field=None)
		map_half2 = hp.read_map(indirectory+maps_half2[i],field=None)
		var_i = (np.abs(map_half1[0] - map_half2[0])/2.0)**2
		var_q = (np.abs(map_half1[1] - map_half2[1])/2.0)**2
		var_u = (np.abs(map_half1[2] - map_half2[2])/2.0)**2
		var_i[var_i == 0.0] = 1e4
		var_q[var_q == 0.0] = 1e4
		var_u[var_u == 0.0] = 1e4
		var_i[var_i < np.median(var_i)] = np.median(var_i)
		var_q[var_q < np.median(var_q)] = np.median(var_q)
		var_u[var_u < np.median(var_u)] = np.median(var_u)
		var_i[var_i > 1e4] = 1e4
		var_q[var_q > 1e4] = 1e4
		var_u[var_u > 1e4] = 1e4
	elif use_weights:
		# Get the variance maps
		var_i = hp.read_map(indirectory+maps[i].replace('60.0s','60.00s').replace('_mKCMBunits','_weight_0_variance'),field=None)
		var_q = hp.read_map(indirectory+maps[i].replace('60.0s','60.00s').replace('_mKCMBunits','_weight_1_variance'),field=None)
		var_u = hp.read_map(indirectory+maps[i].replace('60.0s','60.00s').replace('_mKCMBunits','_weight_2_variance'),field=None)
	elif use_reweight_by_rms:
		map_half1 = hp.read_map(indirectory+maps_half1[i],field=None)
		map_half2 = hp.read_map(indirectory+maps_half2[i],field=None)
		diff_i = (np.abs(map_half1[0] - map_half2[0])/2.0)**2
		diff_q = (np.abs(map_half1[1] - map_half2[1])/2.0)**2
		diff_u = (np.abs(map_half1[2] - map_half2[2])/2.0)**2
		diff_i[diff_i > 1e4] = 0.0
		diff_q[diff_q > 1e4] = 0.0
		diff_u[diff_u > 1e4] = 0.0
		# print(np.std(diff_i[diff_i != 0.0]))
		# print(np.std(diff_q[diff_q != 0.0]))
		# print(np.std(diff_u[diff_u != 0.0]))
		var_i = hp.read_map(indirectory+maps[i].replace('60.0s','60.00s').replace('_mKCMBunits','_weight_0_variance'),field=None)
		var_q = hp.read_map(indirectory+maps[i].replace('60.0s','60.00s').replace('_mKCMBunits','_weight_1_variance'),field=None)
		var_u = hp.read_map(indirectory+maps[i].replace('60.0s','60.00s').replace('_mKCMBunits','_weight_2_variance'),field=None)
		var_i[var_i < -1e4] = 0.0
		var_q[var_q < -1e4] = 0.0
		var_u[var_u < -1e4] = 0.0
		# print(np.max(var_i))
		# print(np.min(var_i))
		noise_i = noiserealisation(np.sqrt(var_i[var_i != 0.0]),len(var_i[var_i != 0.0]))
		noise_q = noiserealisation(np.sqrt(var_q[var_q != 0.0]),len(var_q[var_q != 0.0]))
		noise_u = noiserealisation(np.sqrt(var_u[var_u != 0.0]),len(var_u[var_u != 0.0]))
		# print(np.std(noise_i[noise_i != 0.0]))
		# print(np.std(noise_q[noise_q != 0.0]))
		# print(np.std(noise_u[noise_u != 0.0]))

		rescale_vals[0,i] = np.std(diff_i[diff_i != 0.0])/np.std(noise_i[noise_i != 0.0])
		rescale_vals[1,i] = np.std(diff_q[diff_q != 0.0])/np.std(noise_q[noise_q != 0.0])
		rescale_vals[2,i] = np.std(diff_u[diff_u != 0.0])/np.std(noise_u[noise_u != 0.0])


		var_i[:] = var_i[:] * np.std(diff_i[diff_i != 0.0])/np.std(noise_i[noise_i != 0.0])
		var_q[:] = var_q[:] * np.std(diff_q[diff_q != 0.0])/np.std(noise_q[noise_q != 0.0])
		var_u[:] = var_u[:] * np.std(diff_u[diff_u != 0.0])/np.std(noise_u[noise_u != 0.0])
		var_i[var_i == 0.0] = 1e4
		var_q[var_q == 0.0] = 1e4
		var_u[var_u == 0.0] = 1e4


	hp.mollview(var_i,norm='hist')
	plt.savefig(outdirectory+maps[i]+'_0_var.pdf')
	hp.mollview(var_q,norm='hist')
	plt.savefig(outdirectory+maps[i]+'_1_var.pdf')
	hp.mollview(var_u,norm='hist')
	plt.savefig(outdirectory+maps[i]+'_2_var.pdf')
	# print(maps[i])
	# print(np.median(np.sqrt(var_i[var_i[:] >=0])))
	# print(np.median(np.sqrt(var_q[var_q[:] >=0])))
	# print(np.median(np.sqrt(var_u[var_u[:] >=0])))

	var_i = var_i * ((normfreq/freqs[i])**index)**2
	var_q = var_q * ((normfreq/freqs[i])**index)**2
	var_u = var_u * ((normfreq/freqs[i])**index)**2

	# print(np.median(np.sqrt(var_i[var_i[:] >=0])))
	# print(np.median(np.sqrt(var_q[var_q[:] >=0])))
	# print(np.median(np.sqrt(var_u[var_u[:] >=0])))


	if i != 0 and i != 1:
		if i == 2:
			combine_q = (mapdata[1].copy()*(normfreq/freqs[i])**index)/var_q
			combine_u = (mapdata[2].copy()*(normfreq/freqs[i])**index)/var_u
			weight_q = 1.0/(var_q.copy())
			weight_u = 1.0/(var_u.copy())
		else:
			combine_q = combine_q+(mapdata[1]*(normfreq/freqs[i])**index)/var_q
			combine_u = combine_u+(mapdata[2]*(normfreq/freqs[i])**index)/var_u
			weight_q = weight_q+1.0/var_q
			weight_u = weight_u+1.0/var_u

combine_q /= weight_q
combine_u /= weight_u

hp.write_map(outdirectory+prefix+'_commonmask.fits',commonmask,overwrite=True)
hp.mollview(commonmask)
plt.savefig(outdirectory+prefix+'_commonmask.pdf')

hp.write_map(outdirectory+prefix+'_combine_q.fits',combine_q*commonmask,overwrite=True)
hp.mollview(combine_q*commonmask,min=-1,max=1)
plt.savefig(outdirectory+prefix+'_combine_q.pdf')
hp.write_map(outdirectory+prefix+'_combine_u.fits',combine_u*commonmask,overwrite=True)
hp.mollview(combine_u*commonmask,min=-1,max=1)
plt.savefig(outdirectory+prefix+'_combine_u.pdf')
hp.write_map(outdirectory+prefix+'_combine_P.fits',np.sqrt(combine_q**2+combine_u**2)*commonmask,overwrite=True)
hp.mollview(np.sqrt(combine_q**2+combine_u**2)*commonmask,min=0,max=1.2,cmap=plt.get_cmap('jet'))
plt.savefig(outdirectory+prefix+'_combine_P.pdf')

hp.write_map(outdirectory+prefix+'_combine_P_nomask.fits',np.sqrt(combine_q**2+combine_u**2),overwrite=True)
hp.mollview(np.sqrt(combine_q**2+combine_u**2),min=0,max=1.2,cmap=plt.get_cmap('jet'))
plt.savefig(outdirectory+prefix+'_combine_P_nomask.pdf')

mapdata = hp.read_map('/Users/mpeel/Documents/maps/wmap_planck_pol/weighted_both_debwk_feb2015_tqu.fits',field=None)
mapdata = hp.pixelfunc.reorder(mapdata[0],n2r=True)

commonmask2 = hp.ud_grade(commonmask,256,order_in='RING',order_out='RING')

hp.mollview(mapdata*1000.0*commonmask2,min=0,max=0.03,cmap=plt.get_cmap('jet'))
plt.savefig(outdirectory+'combine_P_planck.pdf')
hp.mollview(mapdata*1000.0,min=0,max=0.06,cmap=plt.get_cmap('jet'))
plt.savefig(outdirectory+'combine_P_planck_nomask.pdf')

np.set_printoptions(formatter={'float': '{: 0.5f}'.format})
print(rescale_vals)