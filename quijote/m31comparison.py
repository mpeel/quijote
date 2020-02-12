#!/usr/bin/env python
# -*- coding: utf-8  -*-
#
# Do a comparison of M31 from QUIJOTE, Planck and WMAP
# 
# Version history:
#
# 09-May-2019  M. Peel       Started

import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
from astrocode.fitspectrum.astroutils import *

qfolder = 'quijote_201911'
qext = 'nov2019'
qfolder = 'quijote_201905'
qext = 'may2019'
map_prefix = '/Users/mpeel/Documents/maps/wmap_planck_smoothed_maps_2015/'
filenames = ['256_60.00smoothed_wmap9decbeamCMBSmicasub_22.8_512_2013_mKCMBunits.fits','256_60.00smoothed_wmap9decbeamCMBSmicasub_33.0_512_2013_mKCMBunits.fits','256_60.00smoothed_wmap9decbeamCMBSmicasub_40.7_512_2013_mKCMBunits.fits','256_60.00smoothed_wmap9decbeamCMBSmicasub_60.7_512_2013_mKCMBunits.fits','256_60.00smoothed_wmap9decbeamCMBSmicasub_93.5_512_2013_mKCMBunits.fits','256_60.00smoothed_PlanckR2fullbeamCMBSmicasub_28.4_1024_2015_mKCMBunits.fits','256_60.00smoothed_PlanckR2fullbeamCMBSmicasub_44.1_1024_2015_mKCMBunits.fits','256_60.00smoothed_PlanckR2fullbeamCMBSmicasub_70.4_1024_2015_mKCMBunits.fits','256_60.00smoothed_PlanckR2fullbeamCMBSmicasub_100_2048_2015_mKCMBunits.fits','256_60.00smoothed_PlanckR2fullbeamCMBSmicasub_143_2048_2015_mKCMBunits.fits','256_60.00smoothed_PlanckR2fullbeamCMBSmicasub_217_2048_2015_mKCMBunits.fits','256_60.00smoothed_PlanckR2fullbeamCMBSmicasub_353_2048_2015_mKCMBunits.fits','../'+qfolder+'/reform/mfi_'+qext+'_map_11.0_1.fits','../'+qfolder+'/reform/mfi_'+qext+'_map_13.0_1.fits','../'+qfolder+'/reform/mfi_'+qext+'_map_17.0_2.fits','../'+qfolder+'/reform/mfi_'+qext+'_map_19.0_2.fits','../'+qfolder+'/reform/mfi_'+qext+'_map_11.0_3.fits','../'+qfolder+'/reform/mfi_'+qext+'_map_11.0_3.fits','../'+qfolder+'/reform/mfi_'+qext+'_map_17.0_4.fits','../'+qfolder+'/reform/mfi_'+qext+'_map_19.0_4.fits']
freqs = [22.8,33.0,40.7,60.7,93.5,28.4,44.1,70.4,100.0,143.0,217.0,353.0,11.0,13.0,17.0,19.0,11.0,13.0,17.0,19.0]
# quijote_map = '/Users/mpeel/Documents/maps/quijote_201907/quijote_combinedmaps_select_withcov_rfi_jul2019_withcov_rfi_nside512_finalmaskfdec2.fits'
# quijote_map = '/Users/mpeel/Documents/QUIJOTE/m31/mapsds_iqu_m31_comb_rasters_931td_rgs_nominal_oct2018_ctod_ns512.fits'
quijote_map = '/Users/mpeel/Documents/QUIJOTE/m31/quijote_combinedmaps_nominal_and_m31_test_rfi_newflags2_may2019_withper5_test_rfi_newflags2_nside512.fits'

outdir = '/Users/mpeel/Documents/QUIJOTE/m31/comparison_zoom_202001/'

cbass_map = '/Users/mpeel/Documents/maps/cbass2019/cbass_global8p8deg_swapQU_NIGHT_v28allelsNs_37_noiseCut_masked5pc_G_1024_ol500_lessTol_g_map_g_1deg_0256.fits'

res_arcmin = 60.0
lon = 121.2
lat = -21.6
aper_inner_radius = 100.0
aper_outer_radius1 = 100.0
aper_outer_radius2 = 140.0
noise_model = 1
units='mK_CMB'

freqs = [0]
# for filename in filenames:
for i in range(0,len(freqs)):
	print(freqs[i])
	print(filenames[i])
	mapdata = hp.read_map(map_prefix + filenames[i])
	print(haperflux(mapdata, freqs[i], res_arcmin, lon, lat, aper_inner_radius, aper_outer_radius1, aper_outer_radius2, units, column=0, dopol=False, nested=False, noise_model=noise_model, abratio=0.7, angle=45.0, silent=False))
	hp.gnomview(mapdata,rot=[10.68, 41.27],coord=['G','C'],title='M31',reso=2)
	plt.savefig(outdir+filenames[i].replace('../'+qfolder+'/reform/','')+'_m31.pdf')
	plt.close()
	plt.clf()

# pixels = query_ellipse(256, 121.2, -21.6, 100.0/60.0, 0.7, 45.0)
# pixels2 = query_ellipse(256, 121.2, -21.6, 140.0/60.0, 0.7, 45.0)

mapdata = hp.read_map(cbass_map)
# mapdata2=mapdata.copy()
# mapdata2[pixels2] = 0.0
# mapdata2[pixels] = mapdata[pixels]
# mapdata=mapdata2
units='K_CMB'
print(haperflux(mapdata, 4.76, res_arcmin, lon, lat, aper_inner_radius, aper_outer_radius1, aper_outer_radius2, units, column=0, dopol=False, nested=False, noise_model=noise_model, abratio=0.7, angle=45.0, silent=False))

hp.mollview(mapdata)
plt.savefig(outdir+'cbass.pdf')
hp.gnomview(mapdata,rot=[10.68, 41.27],coord=['G','C'],title='M31',reso=2)
plt.savefig(outdir+'cbass_m31.pdf')
plt.close()
plt.clf()

inputfits = fits.open(quijote_map)
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
freqs = [11.0,17.0,11.0,17.0,13.0,19.0,13.0,19.0]

# data[1] = pixel numbers
# data[2][0] = I
# data[2][1] = Q
# data[2][2] = U

print(np.shape(data))
print(len(data[1]))
print(len(data[2][0]))

# for i in range(0,len(data[1])):
# 	# print(data[1][i])
# 	# print(data[2][0][i])
# 	# print(len(data[2][0][i]))
# 	# maps[0][data[1][i]] = data[2][0][i][0][0]
# 	# maps[1][data[1][i]] = data[2][0][i][0][1]
# 	# maps[2][data[1][i]] = data[2][0][i][0][2]
# 	# maps[3][data[1][i]] = data[2][0][i][0][3]
# 	# maps[4][data[1][i]] = data[2][0][i][1][0]
# 	# maps[5][data[1][i]] = data[2][0][i][1][1]
# 	# maps[6][data[1][i]] = data[2][0][i][1][2]
# 	# maps[7][data[1][i]] = data[2][0][i][1][3]
# 	for k in range(0,2): # Low and high
# 		for l in range(0,4): # Horns
# 			# maps[k*4+l][0][data[1][i]] = data[0][0][k][l][:]
# 			# maps[k*4+l][1][data[1][i]] = data[0][1][k][l][:]
# 			# maps[k*4+l][2][data[1][i]] = data[0][2][k][l][:]
# 			# These were for the smoothed maps
# 			maps[k*4+l][0][data[1][i]] = data[2][0][i][k][l]
# 			maps[k*4+l][1][data[1][i]] = data[2][1][i][k][l]
# 			maps[k*4+l][2][data[1][i]] = data[2][2][i][k][l]

for k in range(0,2): # Low and high
	for l in range(0,4): # Horns
		maps[k*4+l][0] = data[0][0][k][l]
		maps[k*4+l][1] = data[0][1][k][l]
		maps[k*4+l][2] = data[0][2][k][l]

print(np.shape(maps))
for i in range(0,nummaps):
	for j in range(0,2): # Pol

		# newmap = hp.ud_grade(maps[i][j],256,order_in='NEST',order_out='RING')
		newmap = hp.ud_grade(maps[i][j],256,order_in='RING',order_out='RING')
		hp.write_map(outdir+'mfi_'+str(i+1)+'_'+str(j+1)+'.fits',newmap,overwrite=True)
		hp.mollview(newmap)
		plt.savefig(outdir+'mfi_'+str(i+1)+'_'+str(j+1)+'.pdf')
		plt.close()
		plt.clf()
		hp.gnomview(newmap,rot=[10.68, 41.27],coord=['G','C'],title='M31',reso=2)
		plt.savefig(outdir+'mfi_'+str(i+1)+'_'+str(j+1)+'_m31.pdf')
		plt.close()
		plt.clf()
		units='mK_CMB'
		print(haperflux(newmap/0.9, freqs[i], res_arcmin, lon, lat, aper_inner_radius, aper_outer_radius1, aper_outer_radius2, units, column=0, dopol=False, nested=False, noise_model=noise_model, abratio=0.7, angle=45.0, silent=False))


