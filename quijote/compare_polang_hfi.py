#!/usr/bin/env python
# -*- coding: utf-8  -*-
#
# Do a quick comparison of the Planck HFI polarisation angles
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

indirectory = '/Volumes/Toshiba5TB2/maps/planck2018_tqu/'
outdirectory = '/Users/mpeel/Documents/maps/planck2018_tqu_weight/'

beta_d = 1.53 # From Planck 2018 IV
T_d = 19.6 # Ditto
const = get_spectrum_constants()
val_353 = thermaldust_comm(const, 353.0, 1.0, beta_d, T_d)*planckcorr(const, 353.0)
val_217 = thermaldust_comm(const, 217.0, 1.0, beta_d, T_d)*planckcorr(const, 217.0)
val_143 = thermaldust_comm(const, 143.0, 1.0, beta_d, T_d)*planckcorr(const, 143.0)

# wgt_filename = outdirectory+'planck2018_hfi_tqu_60arcmin_v0_combine.fits'
# wgt = hp.read_map(wgt_filename,field=None)
# polang_map_wgt = calc_polang(wgt[1],wgt[2])
p353_filename = indirectory+'2048_60.0smoothed_PlanckR3fullbeam_353_2048_2018_mKCMBunits.fits'
p353 = hp.read_map(p353_filename,field=None)
# polang_map_p353 = calc_polang(p353[1],p353[2])
p217_filename = indirectory+'2048_60.0smoothed_PlanckR3fullbeam_217_2048_2018_mKCMBunits.fits'
p217 = hp.read_map(p217_filename,field=None)
# polang_map_p217 = calc_polang(p217[1],p217[2])
p143_filename = indirectory+'2048_60.0smoothed_PlanckR3fullbeam_143_2048_2018_mKCMBunits.fits'
p143 = hp.read_map(p143_filename,field=None)
# polang_map_p143 = calc_polang(p143[1],p143[2])

# hp.mollview(polang_map_wgt,cmap=plt.get_cmap('hsv'))
# plt.savefig(outdirectory+'_polang_wgt.pdf')
# hp.mollview(polang_map_p353,cmap=plt.get_cmap('hsv'))
# plt.savefig(outdirectory+'_polang_p353.pdf')
# hp.mollview(polang_map_p217,cmap=plt.get_cmap('hsv'))
# plt.savefig(outdirectory+'_polang_p217.pdf')
# hp.mollview(polang_map_p143,cmap=plt.get_cmap('hsv'))
# plt.savefig(outdirectory+'_polang_p143.pdf')
# plt.clf()

# diff_353_217 = dodiff(polang_map_p353,polang_map_p217)
# diff_353_143 = dodiff(polang_map_p353,polang_map_p143)
# diff_217_143 = dodiff(polang_map_p217,polang_map_p143)

# diff_353_217_val = np.median(diff_353_217)
# diff_353_143_val = np.median(diff_353_143)
# diff_217_143_val = np.median(diff_217_143)
# print(diff_353_217_val)
# print(diff_353_143_val)
# print(diff_217_143_val)

# hp.mollview(diff_353_217,cmap=plt.get_cmap('hsv'))
# plt.savefig(outdirectory+'_polang_p353_p217.pdf')
# hp.mollview(diff_353_143,cmap=plt.get_cmap('hsv'))
# plt.savefig(outdirectory+'_polang_p353_p143.pdf')
# hp.mollview(diff_217_143,cmap=plt.get_cmap('hsv'))
# plt.savefig(outdirectory+'_polang_p217_p143.pdf')

# threshold = 30
# hp.mollview(diff_353_217,cmap=plt.get_cmap('hsv'),min=-threshold,max=threshold)
# plt.savefig(outdirectory+'_polang_p353_p217_cut.pdf')
# hp.mollview(diff_353_143,cmap=plt.get_cmap('hsv'),min=-threshold,max=threshold)
# plt.savefig(outdirectory+'_polang_p353_p143_cut.pdf')
# hp.mollview(diff_217_143,cmap=plt.get_cmap('hsv'),min=-threshold,max=threshold)
# plt.savefig(outdirectory+'_polang_p217_p143_cut.pdf')

# hp.gnomview(diff_353_143,reso=20,min=-15,max=15)
# plt.savefig(outdirectory+'_polang_p353_p143_zoom.pdf')



p217 *= val_353/val_217
p143 *= val_353/val_143

hp.mollview(np.sqrt(((p353[1]-p217[1])**2)+((p353[2]-p217[2])**2)),min=0,max=0.1)
plt.savefig(outdirectory+'_diff_p353_p217.pdf')
plt.clf()
hp.mollview(np.sqrt(((p353[1]-p143[1])**2)+((p353[2]-p143[2])**2)),min=0,max=0.1)
plt.savefig(outdirectory+'_diff_p353_p143.pdf')
plt.clf()
hp.mollview(np.sqrt(((p217[1]-p143[1])**2)+((p217[2]-p143[2])**2)),min=0,max=0.1)
plt.savefig(outdirectory+'_diff_p217_p143.pdf')
plt.clf()

hp.write_map(outdirectory+'_diff_p353_p217.fits',np.sqrt(((p353[1]-p217[1])**2)+((p353[2]-p217[2])**2)),overwrite=True)
hp.write_map(outdirectory+'_diff_p353_p143.fits',np.sqrt(((p353[1]-p143[1])**2)+((p353[2]-p143[2])**2)),overwrite=True)
hp.write_map(outdirectory+'_diff_p217_p143.fits',np.sqrt(((p217[1]-p143[1])**2)+((p217[2]-p143[2])**2)),overwrite=True)

Pp353 = np.sqrt(p353[1]**2+p353[2]**2)
Pp217 = np.sqrt(p217[1]**2+p217[2]**2)
Pp143 = np.sqrt(p143[1]**2+p143[2]**2)
hp.mollview(Pp353-Pp217,min=-0.1,max=0.1)
plt.savefig(outdirectory+'_diff2_p353_p217.pdf')
plt.clf()
hp.mollview(Pp353-Pp143,min=-0.1,max=0.1)
plt.savefig(outdirectory+'_diff2_p353_p143.pdf')
plt.clf()
hp.mollview(Pp217-Pp143,min=-0.1,max=0.1)
plt.savefig(outdirectory+'_diff2_p217_p143.pdf')
plt.clf()

hp.write_map(outdirectory+'_diff2_p353_p217.fits',Pp353-Pp217,overwrite=True)
hp.write_map(outdirectory+'_diff2_p353_p143.fits',Pp353-Pp143,overwrite=True)
hp.write_map(outdirectory+'_diff2_p217_p143.fits',Pp217-Pp143,overwrite=True)
