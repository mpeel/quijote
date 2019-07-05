#!/usr/bin/env python
# -*- coding: utf-8  -*-
#
# Do a quick comparison of the weights and nobs maps
# 
# Version history:
#
# 05-Jul-2019  M. Peel       Started

import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy.coordinates import Galactic, FK5
from astropy import units as u

nside = 512
npix = hp.nside2npix(nside)

indirectory = '/Users/mpeel/Documents/maps/quijote_201904/reform/'
outdirectory = '/Users/mpeel/Documents/maps/quijote_201904/analyse/'
date='201905'

# Create a map that is the declination for given Galactic pixels
declinationmap = np.zeros(npix)
pos = hp.pixelfunc.pix2ang(nside,range(0,npix),lonlat=True)
galcoord = SkyCoord(l=pos[0]*u.degree, b=pos[1]*u.degree, frame='galactic')
newpos = galcoord.transform_to(FK5(equinox='J2000'))
declinationmap = newpos.dec.degree+90.0
hp.mollview(declinationmap)
plt.savefig(outdirectory+'declinationmap.pdf')
numbins = 180
declinationmap_int = declinationmap.astype(int)

mapends = ['_map_13.0_1.fits','_map_13.0_3.fits','_map_17.0_2.fits','_map_17.0_4.fits','_map_19.0_2.fits','_map_19.0_4.fits','_map_11.0_1.fits','_map_11.0_3.fits']
# prefixes = ['mfi_may2019']
prefixes = ['mfi_apr2019_30_2','mfi_apr2019_35_6','mfi_apr2019_40_2','mfi_apr2019_40_5','mfi_apr2019_50_2','mfi_apr2019_50_5','mfi_apr2019_50_6','mfi_apr2019_60_1','mfi_apr2019_60_2','mfi_apr2019_60_5','mfi_apr2019_60_6','mfi_apr2019_65_1','mfi_apr2019_65_2','mfi_apr2019_65_6','mfi_apr2019_70_6']
maps = []
for prefix in prefixes:
	for end in mapends:
		maps.append(prefix+end)

for mapname in maps:
	ratio = np.zeros((3,npix))
	decaverage = np.zeros((3,numbins))
	decaverage_n = np.zeros((3,numbins))
	weights = hp.read_map(indirectory+mapname.replace('map','weights'),field=None)
	nobs = hp.read_map(indirectory+mapname.replace('map','nhits'),field=None)
	print(weights[0])
	ratio = (1.0/(weights**0.5))/(nobs**0.5)
	for i in range(0,3):
		ratio[i][~np.isfinite(ratio[i])] = hp.UNSEEN
		print(np.median(ratio[i][ratio[i] != hp.UNSEEN]))
		for j in range(0,npix):
			if ratio[i][j] != hp.UNSEEN:
				decaverage[i][declinationmap_int[j]] += ratio[i][j]
				decaverage_n[i][declinationmap_int[j]] += 1
		hp.mollview(ratio[i])
		plt.savefig(outdirectory+mapname+'_ratio_'+str(i)+'.pdf')
		plt.close()
		plt.clf()
	for i in range(0,3):
		plt.plot(range(-90,90),decaverage[i]/decaverage_n[i],label=str(i))
	plt.title(mapname)
	plt.savefig(outdirectory+mapname+'_decaverage.pdf')
