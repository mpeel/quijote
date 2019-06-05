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
types = ['map','nhits','weights']
weights_id = 2
invert_weights = True

quijote_map = ['/Users/mpeel/Documents/maps/quijote_201810/quijote_combinedmaps_select_oct2018_test_rfi_alldata_highfk_nside512.fits']
outdir = '/Users/mpeel/Documents/maps/quijote_201810/reform/'
prefix = ['mfi_oct2018']
types = ['map','mapsmth','nhits','weights']
weights_id = 3

for mapnum in range(0,len(quijote_map)):
	inputfits = fits.open(quijote_map[mapnum])
	print(inputfits[1].header)
	col_names = inputfits[1].columns.names
	print(col_names)
	nmaps = len(col_names)
	maps = []
	# print(inputfits[1].data)
	print(inputfits[1].data[0][0])

	nside = 512
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

	print(len(data[0]))
	print(len(data[0][0]))
	print(data[0][0][0][0][0])
	print(np.shape(data[0][0][0][0]))
	count_types = len(types)
	for i in range(0,count_types): # map/nhits/weights
		for k in range(0,2): # Low/high band
			for l in range(0,4): # Horns
				if i == weights_id and invert_weights:
					# Invert the weights map
					data[i][0][k][l][data[1][0][k][l]!=0]	= 1.0 / data[i][0][k][l][data[1][0][k][l]!=0]
					data[i][1][k][l][data[1][0][k][l]!=0]	= 1.0 / data[i][1][k][l][data[1][0][k][l]!=0]
					data[i][2][k][l][data[1][0][k][l]!=0]	= 1.0 / data[i][2][k][l][data[1][0][k][l]!=0]

				data[i][0][k][l][data[1][0][k][l]==0] = hp.pixelfunc.UNSEEN
				data[i][1][k][l][data[1][0][k][l]==0] = hp.pixelfunc.UNSEEN
				data[i][2][k][l][data[1][0][k][l]==0] = hp.pixelfunc.UNSEEN
				outmap = [data[i][0][k][l], data[i][1][k][l], data[i][2][k][l]]
				print(np.shape(outmap))
				outname = prefix[mapnum]+'_'+types[i]+'_'+str(freqs[k][l])+'_'+str(l+1)+'.fits'
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

	# for i in range(0,len(data[1])):
	# 	for k in range(0,2): # Low and high
	# 		for l in range(0,4): # Horns
	# 			# print(i)
	# 			# print(k)
	# 			# print(l)
	# 			maps[k*4+l][0][i] = data[0][i][k][l]
	# 			maps[k*4+l][1][i] = data[1][i][k][l]
	# 			maps[k*4+l][2][i] = data[2][i][k][l]

	# for i in range(0,nummaps):
	# 	for j in range(0,2): # Pol

	# 		newmap = hp.ud_grade(maps[i][j],256,order_in='NEST',order_out='RING')
	# 		hp.write_map(outdir+'mfi_'+str(i+1)+'_'+str(j+1)+'.fits',newmap,overwrite=True)
	# 		hp.mollview(newmap)
	# 		plt.savefig(outdir+'mfi_'+str(i+1)+'_'+str(j+1)+'.pdf')
	# 		plt.close()
	# 		plt.clf()


