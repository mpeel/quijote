#!/usr/bin/env python
# -*- coding: utf-8  -*-
#
# Split the maps
# 
# Version history:
#
# 09-May-2019  M. Peel       Started
# 05-Jun-2019  M. Peel       Generalised to cope with multiple runs
# 22-Mar-2021  M. Peel       Reformat with multiple maps in the same file

import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
from astrocode.fitspectrum.astroutils import *
import os

invert_weights = True
indir = '/Users/mpeel/Documents/maps/quijote_202103/'
outdir = indir+'reform/'
# types = ['map','mapsmth','nhits','weights','covqu']
quijote_map = ['quijote_combinedmaps_nominalselect_recalib11_withcov2_rfi_newflagstat_raw_allatonce5s_mar2021_allatonce_nside512_finalmaskfdec2_newwei.fits','quijote_combinedmaps_nominalselect_recalib11_withcov2_rfi_newflagstat_raw_allatonce5s_mar2021_daynight1_nside512_finalmaskfdec2_newwei.fits','quijote_combinedmaps_nominalselect_recalib11_withcov2_rfi_newflagstat_raw_allatonce5s_mar2021_daynight2_nside512_finalmaskfdec2_newwei.fits','quijote_combinedmaps_nominalselect_recalib11_withcov2_rfi_newflagstat_raw_allatonce5s_mar2021_half1_nside512_finalmaskfdec2_newwei.fits','quijote_combinedmaps_nominalselect_recalib11_withcov2_rfi_newflagstat_raw_allatonce5s_mar2021_half2_nside512_finalmaskfdec2_newwei.fits','quijote_combinedmaps_nominalselect_recalib11_withcov2_rfi_newflagstat_raw_allatonce5s_mar2021_halfring1_nside512_finalmaskfdec2_newwei.fits','quijote_combinedmaps_nominalselect_recalib11_withcov2_rfi_newflagstat_raw_allatonce5s_mar2021_halfring2_nside512_finalmaskfdec2_newwei.fits','quijote_combinedmaps_nominalselect_recalib11_withcov2_rfi_newflagstat_raw_allatonce5s_mar2021_period1_nside512_finalmaskfdec2_newwei.fits','quijote_combinedmaps_nominalselect_recalib11_withcov2_rfi_newflagstat_raw_allatonce5s_mar2021_period2_nside512_finalmaskfdec2_newwei.fits','quijote_combinedmaps_nominalselect_recalib11_withcov2_rfi_newflagstat_raw_allatonce5s_mar2021_period5_nside512_finalmaskfdec2_newwei.fits','quijote_combinedmaps_nominalselect_recalib11_withcov2_rfi_newflagstat_raw_allatonce5s_mar2021_period6_nside512_finalmaskfdec2_newwei.fits','quijote_combinedmaps_nominalselect_recalib11_withcov2_rfi_newflagstat_raw_allatonce5s_mar2021_pwv1_nside512_finalmaskfdec2_newwei.fits','quijote_combinedmaps_nominalselect_recalib11_withcov2_rfi_newflagstat_raw_allatonce5s_mar2021_pwv2_nside512_finalmaskfdec2_newwei.fits']
prefix = ['mfi_mar2021','mfi_mar2021_daynight1','mfi_mar2021_daynight2','mfi_mar2021_half1','mfi_mar2021_half2','mfi_mar2021_halfring1','mfi_mar2021_halfring2','mfi_mar2021_period1','mfi_mar2021_period2','mfi_mar2021_period5','mfi_mar2021_period6','mfi_mar2021_pwv1','mfi_mar2021_pwv2']
unit = 'mKCMB'
scale = [[1.0, 1.0, 1.0, 1.0],[1.0, 1.0, 1.0, 1.0]]

print(len(quijote_map))
print(len(prefix))
os.makedirs(outdir, exist_ok=True)
for mapnum in range(0,len(quijote_map)):
	print(quijote_map[mapnum], prefix[mapnum])

	inputfits = fits.open(indir+quijote_map[mapnum])
	print(inputfits[1].header)

	col_names = inputfits[1].columns.names
	print(col_names)
	nmaps = len(col_names)
	maps = []

	nside = 512
	npix = 12*nside**2
	nummaps = 8
	pol = 3
	maps = np.zeros((nummaps,pol,npix))
	data = inputfits[1].data[0]
	print(np.shape(data))

	# 1 = 11/13GHz, broken
	# 2 = 17/19GHz, working
	# 3 = 11/13GHz, working
	# 4 = 17/19GHz, working
	freqs = [[11.0,17.0,11.0,17.0],[13.0,19.0,13.0,19.0]]

	# data[1] = pixel numbers
	# data[2][0] = I
	# data[2][1] = Q
	# data[2][2] = U

	# print(len(data[0]))
	# print(len(data[0][0]))
	# print(data[0][0][0][0][0])
	# print(np.shape(data[0][0][0][0]))
	for k in range(0,2): # Low/high band
		for l in range(0,4): # Horns
			outname = prefix[mapnum]+'_'+str(freqs[k][l])+'_'+str(l+1)+'.fits'
			indexmap = 0
			indexwgt = 0
			indexcov = 0
			indexsmth = -1
			# We want to alter the maps a bit
			for i in range(0,len(col_names)): # map/nhits/weights
				if col_names[i] == 'MAP':
					indexmap = i
				elif col_names[i] == 'WEI':
					indexwgt = i
				elif col_names[i] == 'COV_QU':
					indexcov = i
				elif col_names[i] == 'MAP_SM1D':
					indexsmth = i

			for i in range(0,len(col_names)): # map/nhits/weights
				print(outname)
				if col_names[i] in ['MAP','MAP_SM1D']:
					print('Rescaling by ' + str(scale[k][l]))
					data[i][0][k][l][:] /= scale[k][l]
					data[i][1][k][l][:] /= scale[k][l]
					data[i][2][k][l][:] /= scale[k][l]
				if col_names[i] == 'WEI':
					print('Rescaling by ' + str((1.0/scale[k][l])**2))
					data[i][0][k][l][:] *= scale[k][l]**2
					data[i][1][k][l][:] *= scale[k][l]**2
					data[i][2][k][l][:] *= scale[k][l]**2
				if col_names[i] == 'WEI' and invert_weights:
					# Invert the weights map
					data[i][0][k][l][data[1][0][k][l]!=0] = 1.0 / data[i][0][k][l][data[1][0][k][l]!=0]
					data[i][1][k][l][data[1][0][k][l]!=0] = 1.0 / data[i][1][k][l][data[1][0][k][l]!=0]
					data[i][2][k][l][data[1][0][k][l]!=0] = 1.0 / data[i][2][k][l][data[1][0][k][l]!=0]

				if col_names[i] != 'COV_QU':
					data[i][0][k][l][data[indexwgt][0][k][l]==0] = hp.pixelfunc.UNSEEN
					data[i][1][k][l][data[indexwgt][1][k][l]==0] = hp.pixelfunc.UNSEEN
					data[i][2][k][l][data[indexwgt][2][k][l]==0] = hp.pixelfunc.UNSEEN
					data[i][0][k][l][data[indexwgt][0][k][l]==hp.pixelfunc.UNSEEN] = hp.pixelfunc.UNSEEN
					data[i][1][k][l][data[indexwgt][1][k][l]==hp.pixelfunc.UNSEEN] = hp.pixelfunc.UNSEEN
					data[i][2][k][l][data[indexwgt][2][k][l]==hp.pixelfunc.UNSEEN] = hp.pixelfunc.UNSEEN
					data[i][0][k][l][~np.isfinite(data[i][0][k][l])] = hp.pixelfunc.UNSEEN
					data[i][1][k][l][~np.isfinite(data[i][1][k][l])] = hp.pixelfunc.UNSEEN
					data[i][2][k][l][~np.isfinite(data[i][2][k][l])] = hp.pixelfunc.UNSEEN
				else:
					data[i][k][l][data[indexwgt][1][k][l]==0] = hp.pixelfunc.UNSEEN
					data[i][k][l][data[indexwgt][1][k][l]==hp.pixelfunc.UNSEEN] = hp.pixelfunc.UNSEEN
					data[i][k][l][~np.isfinite(data[i][k][l])] = hp.pixelfunc.UNSEEN

			# Assemble the output map
			cols = []
			cols.append(fits.Column(name='I', format='E', array=np.asarray(data[indexmap][0][k][l])))
			cols.append(fits.Column(name='Q', format='E', array=np.asarray(data[indexmap][1][k][l])))
			cols.append(fits.Column(name='U', format='E', array=np.asarray(data[indexmap][2][k][l])))
			cols.append(fits.Column(name='II_cov', format='E', array=np.asarray(data[indexwgt][0][k][l])))
			cols.append(fits.Column(name='QQ_cov', format='E', array=np.asarray(data[indexwgt][1][k][l])))
			cols.append(fits.Column(name='QU_cov', format='E', array=np.asarray(data[indexcov][k][l])))
			cols.append(fits.Column(name='UU_cov', format='E', array=np.asarray(data[indexwgt][2][k][l])))

			cols = fits.ColDefs(cols)
			bin_hdu = fits.BinTableHDU.from_columns(cols)
			bin_hdu.header['ORDERING']='RING'
			bin_hdu.header['POLCONV']='COSMO'
			bin_hdu.header['PIXTYPE']='HEALPIX'
			bin_hdu.header['COMMENT']=quijote_map[mapnum]
			bin_hdu.header['NSIDE'] = nside
			bin_hdu.header['TUNIT1'] = unit
			bin_hdu.header['TUNIT2'] = unit
			bin_hdu.header['TUNIT3'] = unit
			bin_hdu.header['TUNIT4'] = unit
			bin_hdu.header['TUNIT5'] = unit
			bin_hdu.header['TUNIT6'] = unit
			bin_hdu.header['TUNIT7'] = unit
			bin_hdu.header.comments['TTYPE1'] = 'Intensity brightness temperature'
			bin_hdu.header.comments['TTYPE2'] = 'Q polarisation brightness temperature'
			bin_hdu.header.comments['TTYPE3'] = 'U polarisation brightness temperature'
			bin_hdu.header.comments['TTYPE4'] = 'Intensity covariance'
			bin_hdu.header.comments['TTYPE5'] = 'QQ covariance'
			bin_hdu.header.comments['TTYPE6'] = 'QU covariance'
			bin_hdu.header.comments['TTYPE7'] = 'UU covariance'
			bin_hdu.header['POLAR'] = 'T'
			bin_hdu.header['POLCCONV'] = 'COSMO'
			bin_hdu.header['OBJECT'] = 'FULLSKY'
			bin_hdu.header['FIRSTPIX'] = 0
			bin_hdu.header['LASTPIX'] = len(data[indexmap][0][k][l])
			bin_hdu.header['INDXSCHM'] = 'IMPLICIT'
			bin_hdu.header['BAD_DATA'] = (hp.UNSEEN, 'standard healpix value')
			bin_hdu.header['TELESCOP'] = 'QTMFI1'

			bin_hdu.writeto(outdir+outname,overwrite=True)
			# exit()
			if indexsmth > 0:
				outmap = [data[indexsmth][0][k][l], data[indexsmth][1][k][l], data[indexsmth][2][k][l]]
				hp.write_map(outdir+outname.replace('.fits','_smth.fits'),outmap,overwrite=True)


# EOF