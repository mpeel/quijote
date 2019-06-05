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

nside = 512
npix = hp.nside2npix(nside)

indirectory = '/Users/mpeel/Documents/maps/quijote_201905/smooth/'
outdirectory = '/Users/mpeel/Documents/maps/quijote_201905/analyse/'
prefix='mfi'
date='201905'
# prefix='half1mfi'
# prefix='half2mfi'

# indirectory = '/Users/mpeel/Documents/maps/quijote_201810/smooth/'
# outdirectory = '/Users/mpeel/Documents/maps/quijote_201810/analyse/'
# prefix='mfi'
# date='201810'

maps = [str(nside)+'_60.00smoothed_'+prefix+'1_11.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'1_13.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'2_17.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'2_19.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'3_11.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'3_13.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'4_17.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'4_19.0_512_'+date+'_mKCMBunits.fits']


nummaps = len(maps)
freqs = [11,13,17,19,11,13,17,19]
normfreq = 28.4
index = 3.0
commonmask = np.ones(npix)
combine_q = np.zeros(npix)
combine_u = np.zeros(npix)
for i in range(2,nummaps):
	mapdata = hp.read_map(indirectory+maps[i],field=None)
	commonmask[mapdata[0][:] == hp.UNSEEN] = 0
	hp.mollview(mapdata[0],norm='hist')
	plt.savefig(outdirectory+maps[i]+'_0.pdf')
	hp.mollview(mapdata[1],norm='hist')
	plt.savefig(outdirectory+maps[i]+'_1.pdf')
	hp.mollview(mapdata[2],norm='hist')
	plt.savefig(outdirectory+maps[i]+'_2.pdf')
	hp.mollview(np.sqrt(mapdata[1]**2+mapdata[2]**2),min=0,max=1)
	plt.savefig(outdirectory+maps[i]+'_P.pdf')

	if i != 0 and i != 1:
		if i == 2:
			combine_q = mapdata[1].copy()*(normfreq/freqs[i])**index
			combine_u = mapdata[2].copy()*(normfreq/freqs[i])**index
		else:
			combine_q = combine_q+mapdata[1]*(normfreq/freqs[i])**index
			combine_u = combine_u+mapdata[2]*(normfreq/freqs[i])**index

combine_q /= 6.0
combine_u /= 6.0

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
hp.mollview(np.sqrt(combine_q**2+combine_u**2)*commonmask,min=0,max=2.0,cmap=plt.get_cmap('jet'))
plt.savefig(outdirectory+prefix+'_combine_P.pdf')

mapdata = hp.read_map('/Users/mpeel/Documents/maps/wmap_planck_pol/weighted_both_debwk_feb2015_tqu.fits',field=None)
mapdata = hp.pixelfunc.reorder(mapdata[0],n2r=True)

commonmask2 = hp.ud_grade(commonmask,256,order_in='RING',order_out='RING')

hp.mollview(mapdata*1000.0*commonmask2,min=0,max=0.03,cmap=plt.get_cmap('jet'))
plt.savefig(outdirectory+'combine_P_planck.pdf')
