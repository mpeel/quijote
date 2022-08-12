#!/usr/bin/env python
# -*- coding: utf-8  -*-
#
# Create various plots for the spurs paper
#
# Version history:
#
# 27-Apr-2019  M. Peel       Started
# 20-Sep-2021  M. Peel		 Tidying
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
from astrocode.astroutils import *
import matplotlib
from configure import *
from astrocode.polfunc import *

# Location of data
basedir = '/Users/mpeel/Documents/maps/'
outdir=basedir+'quijote_202103_tqu_v1.5_noise_v1.0_newwf/spurs/'

## Planck/WMAP
plancknpipe30 = 'planck2020_tqu_v1.5_noise_v1.0_10k/256_60.0smoothed_PlanckR4fullbeamnodpNoise_28.4_1024_2020_mKCMBunits.fits'
wmapk9 = 'wmap9_tqu_v1.5_noise_v1.0_10k/256_60.0smoothed_wmap9beamNoise_22.8_512_2013_mKCMBunits.fits'
planck_wmap_weighted = 'weighted_wmap_planck_qu_tqu_v1.5_noise_v1.0/256_60.0smoothed_wmap9_planck2020_tqu_v1.5_noise_v1.0_-3.1_combine.fits'
planck_wmap_weighted_2018 = 'weighted_wmap_planck_qu_tqu_v1.5_noise_v1.0/256_60.0smoothed_wmap9_planck2018_tqu_v1.5_noise_v1.0_-3.1_combine.fits'
planck_wmap_weighted_2015 = 'weighted_wmap_planck_qu_tqu_v1.5_noise_v1.0/256_60.0smoothed_wmap9_planck2015_tqu_v1.5_noise_v1.0_-3.1_combine.fits'
planck_wmap_weighted_fdec = 'quijote_202103_tqu_v1.5_noise_v1.0_weighted_fdec/wmapplanck_combine.fits'

## MFI
mfi_weighted = 'quijote_202103_tqu_v1.5_noise_v1.0_weighted/256_60.0smoothed_quijotecombwei10_tqu_v1.5_noise_v1.0_-3.1_combine.fits'
# mfi_planck_weighted = 'quijote_202103_tqu_v1.5_noise_v1.0_weighted/256_60.0smoothed_quijotecombwei10_tqu_v1.5_noise_v1.0_-3.1_wmapplanck_combine.fits'
mfi_planck_weighted = 'quijote_202103_tqu_v1.5_noise_v1.0_weighted_fdec/512_60.0smoothed_quijotecombwei10_tqu_v1.5_noise_v1.0_-3.1_wmapplanck_qtmask_combine.fits'
mfi_maps = ['quijote_202103_tqu_v1.5_noise_v1.0_newwf/256_60.0smoothed_QUIJOTEMFI2_17.0_2021_mKCMBunits.fits',\
'quijote_202103_tqu_v1.5_noise_v1.0_newwf/256_60.0smoothed_QUIJOTEMFI2_19.0_2021_mKCMBunits.fits',\
'quijote_202103_tqu_v1.5_noise_v1.0_newwf/256_60.0smoothed_QUIJOTEMFI3_11.0_2021_mKCMBunits.fits',\
'quijote_202103_tqu_v1.5_noise_v1.0_newwf/256_60.0smoothed_QUIJOTEMFI3_13.0_2021_mKCMBunits.fits',\
'quijote_202103_tqu_v1.5_noise_v1.0_newwf/256_60.0smoothed_QUIJOTEMFI4_17.0_2021_mKCMBunits.fits',\
'quijote_202103_tqu_v1.5_noise_v1.0_newwf/256_60.0smoothed_QUIJOTEMFI4_19.0_2021_mKCMBunits.fits']
mfi_mask_filename = basedir+'quijote_masks/mask_quijote_ncp_satband_nside512.fits'
# mfi_mask_filename = basedir+'quijote_201907/weighted/mfi_commonmask.fits'

# Configuration
mfi_freqs = [16.7, 18.7, 11.1, 12.9, 17, 19]
plotmax_p = 1.5
plotmax_qu = 1.0
plotmax = 0.5
plotmax_sub = 0.5
spectrum = -3.0
fontsize = 18
reffreq = 10.0
planck_wmap_weighted_reffreq = 28.4
matplotlib.rcParams.update({'font.size':fontsize})

dofigs = [2]
dofigs = [1,2,3,4,5]

regions = get_regions()

# Read in the MFI mask
mfi_mask_512 = hp.read_map(mfi_mask_filename)
mfi_mask = hp.ud_grade(mfi_mask_512,256)
mfi_mask[mfi_mask != 1.0] = 0.0
mfi_mask[:] = 1.0

# Figure 1 WMAP+Planck
if 1 in dofigs:
	# WMAP K
	wmapk9_map = hp.read_map(basedir+wmapk9,field=None)
	wmapk9_pol_10ghz = np.sqrt(wmapk9_map[1]**2+wmapk9_map[2]**2)*(reffreq/22.8)**spectrum
	hp.mollview(wmapk9_pol_10ghz,min=0,max=plotmax_p,cmap='jet',unit='mK CMB',title='')#,title='WMAP 22.8GHz polarised intensity (rescaled to 10GHz)')
	plt.savefig(outdir+'fig1_wmapk_pol.pdf')
	plt.clf()
	plt.close()
	# Planck 30
	plancknpipe30_map = hp.read_map(basedir+plancknpipe30,field=None)
	plancknpipe30_pol_10ghz = np.sqrt(plancknpipe30_map[1]**2+plancknpipe30_map[2]**2)*(reffreq/28.4)**spectrum
	hp.mollview(plancknpipe30_pol_10ghz,min=0,max=plotmax_p,cmap='jet',unit='mK CMB',title='')#,title='Planck 28.4GHz polarised intensity (rescaled to 10GHz)')
	plt.savefig(outdir+'fig1_planck30_pol.pdf')
	plt.clf()
	plt.close()

# Figure 2
if 2 in dofigs:
	# Weighted MFI map
	mfiw = hp.read_map(basedir+mfi_weighted,field=None)
	mfiw_pol = np.sqrt(mfiw[1]**2+mfiw[2]**2)
	mfiw_pol[mfiw[1]==0] = hp.UNSEEN
	hp.mollview(mfiw_pol,min=0,max=plotmax_p,cmap='jet',unit='mK CMB',title='')#,title='MFI weighted polarised intensity'
	plt.savefig(outdir+'fig2_mfi_combine.pdf')
	plt.clf()
	plt.close()
	mfiw_polang = calc_polang(mfiw[1],mfiw[2])
	mfiw_polang[mfiw[1]==0] = hp.UNSEEN
	hp.mollview(mfiw_polang,min=-90,max=90,cmap='hsv',unit='deg',title='')
	plt.savefig(outdir+'fig2_mfi_combine_polang.pdf')
	plt.clf()
	plt.close()

# Figure 3
if 3 in dofigs:
	# Weighted Planck+WMAP map
	planck_wmap = hp.read_map(basedir+planck_wmap_weighted,field=None)
	planck_wmap_pol = np.sqrt(planck_wmap[1]**2+planck_wmap[2]**2)
	planck_wmap_pol_rescale = np.sqrt((planck_wmap[1] * (reffreq/planck_wmap_weighted_reffreq)**spectrum)**2+(planck_wmap[2] * (reffreq/planck_wmap_weighted_reffreq)**spectrum)**2)

	planck_wmap_2018 = hp.read_map(basedir+planck_wmap_weighted_2018,field=None)
	planck_wmap_pol_2018_rescale = np.sqrt((planck_wmap_2018[1] * (reffreq/planck_wmap_weighted_reffreq)**spectrum)**2+(planck_wmap_2018[2] * (reffreq/planck_wmap_weighted_reffreq)**spectrum)**2)

	planck_wmap_2015 = hp.read_map(basedir+planck_wmap_weighted_2015,field=None)
	planck_wmap_pol_2015_rescale = np.sqrt((planck_wmap_2015[1] * (reffreq/planck_wmap_weighted_reffreq)**spectrum)**2+(planck_wmap_2015[2] * (reffreq/planck_wmap_weighted_reffreq)**spectrum)**2)

	# NB: fdec map is at Nside 512, and is already at 10GHz.
	planck_wmap_fdec = hp.read_map(basedir+planck_wmap_weighted_fdec,field=None)
	planck_wmap_fdec = hp.ud_grade(planck_wmap_fdec,nside_out=256)
	planck_wmap_pol_fdec_rescale = np.sqrt((planck_wmap_fdec[1])**2+(planck_wmap_fdec[2])**2)

	hp.mollview(planck_wmap_pol_rescale,min=0,max=plotmax_p,cmap='jet',title='',unit='mK CMB')#,title='WMAP+Planck weighted polarised intensity'
	plt.savefig(outdir+'fig3_wmap_planck.pdf')
	plt.clf()
	hp.mollview(planck_wmap_pol_rescale-planck_wmap_pol_2018_rescale,min=-plotmax_sub,max=plotmax_sub,cmap='jet',title='',unit='mK CMB')#,title='WMAP+Planck weighted polarised intensity'
	plt.savefig(outdir+'fig3_wmap_planck_2020_2018.pdf')
	plt.clf()
	hp.mollview(planck_wmap_pol_2018_rescale-planck_wmap_pol_2015_rescale,min=-plotmax_sub,max=plotmax_sub,cmap='jet',title='',unit='mK CMB')#,title='WMAP+Planck weighted polarised intensity'
	plt.savefig(outdir+'fig3_wmap_planck_2018_2015.pdf')
	plt.clf()
	hp.mollview(planck_wmap_pol_rescale-planck_wmap_pol_2015_rescale,min=-plotmax_sub,max=plotmax_sub,cmap='jet',title='',unit='mK CMB')#,title='WMAP+Planck weighted polarised intensity'
	plt.savefig(outdir+'fig3_wmap_planck_2020_2015.pdf')
	plt.clf()

# for region in regions:
# 	reso = 5.0
# 	hp.gnomview(mfipw_pol,rot=[region[1],region[2]],title=region[0],reso=reso)
# 	hp.graticule()
# 	plt.savefig(outdir+region[0]+'_mfipw.pdf')
# exit()

# Figures 1 and 4, MFI parts
if 1 in dofigs or 4 in dofigs:
	# Read in MFI
	mfi311 = hp.read_map(basedir+mfi_maps[2],field=None)
	mfi311_pol = np.sqrt(mfi311[1]**2+mfi311[2]**2)*(reffreq/mfi_freqs[2])**spectrum
	mfi311_pol[mfi311[1]==hp.UNSEEN] = hp.UNSEEN
	if 1 in dofigs:
		hp.mollview(mfi311_pol,min=0,max=plotmax_p,cmap='jet',unit='mK CMB',title='')#,title='MFI horn 3 11GHz polarised intensity'
		plt.savefig(outdir+'fig1_mfi311_pol.pdf')
		plt.clf()
		plt.close()
	# polfrac = mfi311_pol/mfi311[0]
	# polfrac[mfi_mask==0] = hp.UNSEEN
	# polfrac[polfrac < 0.0] = hp.UNSEEN
	# polfrac[polfrac > 20.0] = hp.UNSEEN
	# # hp.mollview(polfrac,min=0,max=0.7,cmap='jet',title='MFI horn 3 11GHz polarised fraction',unit='mK CMB')
	# plt.savefig(outdir+'polfracmfi311.pdf')
	# exit()
	if 4 in dofigs:
		mfi311_sub = mfi311_pol - planck_wmap_pol_fdec_rescale#planck_wmap_pol_rescale
		mfi311_sub[mfi311[1]==hp.UNSEEN] = hp.UNSEEN
		hp.mollview(mfi311_sub,min=-plotmax_sub,max=plotmax_sub,cmap='jet',unit='mK CMB',title='')#,title='MFI horn 3 11GHz polarised intensity @ 10GHz - Planck+WMAP extrapolated'
		plt.savefig(outdir+'fig4_mfi311_subP.pdf')
		plt.clf()

		mfi311_sub = np.sqrt(((mfi311[1] * (reffreq/mfi_freqs[2])**spectrum) - planck_wmap_fdec[1])**2+((mfi311[2] * (reffreq/mfi_freqs[2])**spectrum) - planck_wmap_fdec[2])**2)
		# mfi311_sub = np.sqrt(((mfi311[1] * (reffreq/mfi_freqs[2])**spectrum) - planck_wmap[1] * (reffreq/planck_wmap_weighted_reffreq)**spectrum)**2+((mfi311[2] * (reffreq/mfi_freqs[2])**spectrum) - planck_wmap[2] * (reffreq/planck_wmap_weighted_reffreq)**spectrum)**2)
		mfi311_sub[mfi311[1]==hp.UNSEEN] = hp.UNSEEN
		hp.mollview(mfi311_sub,min=0,max=plotmax,cmap='jet',unit='mK CMB',title='')#,title='MFI horn 3 11GHz polarised intensity @ 10GHz - Planck+WMAP extrapolated (Q-Q, U-U)'
		plt.savefig(outdir+'fig4_mfi311_subP2.pdf')
		plt.clf()

	# mfi311_sub = hp.sphtfunc.smoothing(mfi311_sub,fwhm=np.sqrt(2.0**2-1.0**2)*np.pi/180.0)
	# # mfi311_sub[mfi_mask==0] = hp.UNSEEN
	# mfi311_sub[mfi311[1]==hp.UNSEEN] = hp.UNSEEN
	# hp.mollview(mfi311_sub,min=0,max=plotmax,cmap='jet',title='MFI horn 3 11GHz polarised intensity @ 10GHz - Planck+WMAP extrapolated (Q-Q, U-U)',unit='mK CMB')
	# plt.savefig(outdir+'mfi311_subP2_smth.pdf')
	# plt.clf()
	# exit()

	mfi313 = hp.read_map(basedir+mfi_maps[3],field=None)
	mfi313_pol = np.sqrt(mfi313[1]**2+mfi313[2]**2)*(reffreq/mfi_freqs[3])**spectrum
	mfi313_pol[mfi313[1]==hp.UNSEEN] = hp.UNSEEN
	if 1 in dofigs:
		hp.mollview(mfi313_pol,min=0,max=plotmax_p,cmap='jet',unit='mK CMB',title='')#,title='MFI horn 3 13GHz polarised intensity'
		plt.savefig(outdir+'fig1_mfi313_pol.pdf')
		plt.clf()
	if 4 in dofigs:
		mfi313_sub = mfi313_pol - planck_wmap_pol_fdec_rescale#planck_wmap_pol_rescale
		mfi313_sub[mfi313[1]==hp.UNSEEN] = hp.UNSEEN
		hp.mollview(mfi313_sub,min=-plotmax_sub,max=plotmax_sub,cmap='jet',unit='mK CMB',title='')#,title='MFI horn 3 13GHz polarised intensity @ 10GHz - Planck+WMAP extrapolated'
		plt.savefig(outdir+'fig4_mfi313_subP.pdf')
		plt.clf()
		mfi313_sub = np.sqrt(((mfi313[1] * (reffreq/mfi_freqs[3])**spectrum) - planck_wmap_fdec[1])**2+((mfi313[2] * (reffreq/mfi_freqs[3])**spectrum) - planck_wmap_fdec[2])**2)
		# mfi313_sub = np.sqrt(((mfi313[1] * (reffreq/mfi_freqs[3])**spectrum) - planck_wmap[1] * (reffreq/planck_wmap_weighted_reffreq)**spectrum)**2+((mfi313[2] * (reffreq/mfi_freqs[3])**spectrum) - planck_wmap[2] * (reffreq/planck_wmap_weighted_reffreq)**spectrum)**2)
		mfi313_sub[mfi313[1]==hp.UNSEEN] = hp.UNSEEN
		hp.mollview(mfi313_sub,min=0,max=plotmax,cmap='jet',unit='mK CMB',title='')#,title='MFI horn 3 13GHz polarised intensity @ 10GHz - Planck+WMAP extrapolated (Q-Q, U-U)'
		plt.savefig(outdir+'fig4_mfi313_subP2.pdf')
		plt.clf()

	mfi217 = hp.read_map(basedir+mfi_maps[0],field=None)
	mfi217_pol = np.sqrt(mfi217[1]**2+mfi217[2]**2)*(reffreq/mfi_freqs[0])**spectrum
	mfi217_pol[mfi217[1]==hp.UNSEEN] = hp.UNSEEN
	if 1 in dofigs:
		hp.mollview(mfi217_pol,min=0,max=plotmax_p,cmap='jet',unit='mK CMB',title='')#,title='MFI horn 2 17GHz polarised intensity'
		plt.savefig(outdir+'fig1_mfi217_pol.pdf')
		plt.clf()

	mfi219 = hp.read_map(basedir+mfi_maps[1],field=None)
	mfi219_pol = np.sqrt(mfi219[1]**2+mfi219[2]**2)*(reffreq/mfi_freqs[1])**spectrum
	mfi219_pol[mfi219[1]==hp.UNSEEN] = hp.UNSEEN
	if 1 in dofigs:
		hp.mollview(mfi219_pol,min=0,max=plotmax_p,cmap='jet',unit='mK CMB',title='')#,title='MFI horn 2 19GHz polarised intensity'
		plt.savefig(outdir+'fig1_mfi219_pol.pdf')
		plt.clf()

	mfi417 = hp.read_map(basedir+mfi_maps[4],field=None)
	mfi417_pol = np.sqrt(mfi417[1]**2+mfi417[2]**2)*(reffreq/mfi_freqs[4])**spectrum
	mfi417_pol[mfi417[1]==hp.UNSEEN] = hp.UNSEEN
	if 1 in dofigs:
		hp.mollview(mfi417_pol,min=0,max=plotmax_p,cmap='jet',unit='mK CMB',title='')#,title='MFI horn 4 17GHz polarised intensity'
		plt.savefig(outdir+'fig1_mfi417_pol.pdf')
		plt.clf()
	if 4 in dofigs:
		mfi417_sub = mfi417_pol - planck_wmap_pol_fdec_rescale#planck_wmap_pol_rescale
		mfi417_sub[mfi417[1]==hp.UNSEEN] = hp.UNSEEN
		hp.mollview(mfi417_sub,min=-plotmax_sub,max=plotmax_sub,cmap='jet',unit='mK CMB',title='')#,title='MFI horn 4 17GHz polarised intensity @ 10GHz - Planck+WMAP extrapolated'
		plt.savefig(outdir+'fig4_mfi417_subP.pdf')
		plt.clf()
		mfi417_sub = np.sqrt(((mfi417[1] * (reffreq/mfi_freqs[4])**spectrum) - planck_wmap_fdec[1])**2+((mfi417[2] * (reffreq/mfi_freqs[4])**spectrum) - planck_wmap_fdec[2])**2)
		# mfi417_sub = np.sqrt(((mfi417[1] * (reffreq/mfi_freqs[4])**spectrum) - planck_wmap[1] * (reffreq/planck_wmap_weighted_reffreq)**spectrum)**2+((mfi417[2] * (reffreq/mfi_freqs[4])**spectrum) - planck_wmap[2] * (reffreq/planck_wmap_weighted_reffreq)**spectrum)**2)
		mfi417_sub[mfi417[1]==hp.UNSEEN] = hp.UNSEEN
		hp.mollview(mfi417_sub,min=0,max=plotmax,cmap='jet',unit='mK CMB',title='')#,title='MFI horn 4 17GHz polarised intensity @ 10GHz - Planck+WMAP extrapolated (Q-Q, U-U)'
		plt.savefig(outdir+'fig4_mfi417_subP2.pdf')
		plt.clf()

	mfi419 = hp.read_map(basedir+mfi_maps[5],field=None)
	mfi419_pol = np.sqrt(mfi419[1]**2+mfi419[2]**2)*(reffreq/mfi_freqs[5])**spectrum
	mfi419_pol[mfi419[1]==hp.UNSEEN] = hp.UNSEEN
	if 1 in dofigs:
		hp.mollview(mfi419_pol,min=0,max=plotmax_p,cmap='jet',unit='mK CMB',title='')#,title='MFI horn 4 19GHz polarised intensity'
		plt.savefig(outdir+'fig1_mfi419_pol.pdf')
		plt.clf()
		plt.close()
	if 4 in dofigs:
		mfi419_sub = mfi419_pol - planck_wmap_pol_fdec_rescale#planck_wmap_pol_rescale
		mfi419_sub[mfi419[1]==hp.UNSEEN] = hp.UNSEEN
		hp.mollview(mfi419_sub,min=-plotmax_sub,max=plotmax_sub,cmap='jet',unit='mK CMB',title='')#,title='MFI horn 4 19GHz polarised intensity @ 10GHz - Planck+WMAP extrapolated'
		plt.savefig(outdir+'fig4_mfi419_subP.pdf')
		plt.clf()
		plt.close()
		mfi419_sub = np.sqrt(((mfi419[1] * (reffreq/mfi_freqs[5])**spectrum) - planck_wmap_fdec[1])**2+((mfi419[2] * (reffreq/mfi_freqs[5])**spectrum) - planck_wmap_fdec[2])**2)
		# mfi419_sub = np.sqrt(((mfi419[1] * (reffreq/mfi_freqs[5])**spectrum) - planck_wmap[1] * (reffreq/planck_wmap_weighted_reffreq)**spectrum)**2+((mfi419[2] * (reffreq/mfi_freqs[5])**spectrum) - planck_wmap[2] * (reffreq/planck_wmap_weighted_reffreq)**spectrum)**2)
		mfi419_sub[mfi419[1]==hp.UNSEEN] = hp.UNSEEN
		hp.mollview(mfi419_sub,min=0,max=plotmax,cmap='jet',unit='mK CMB',title='')#,title='MFI horn 4 19GHz polarised intensity @ 10GHz - Planck+WMAP extrapolated (Q-Q, U-U)'
		plt.savefig(outdir+'fig4_mfi419_subP2.pdf')
		plt.clf()
		plt.close()

# Figure 5
if 5 in dofigs:
	# Weighted Planck+WMAP+MFI map
	mfipw = hp.read_map(basedir+mfi_planck_weighted,field=None)
	mfipw_pol = np.sqrt(mfipw[1]**2+mfipw[2]**2)
	# mfipw_pol[mfi_mask_512==0] = hp.UNSEEN
	hp.mollview(mfipw_pol,min=0,max=plotmax_p,cmap='jet',unit='mK CMB',title='')#,title='MFI+WMAP+Planck weighted polarised intensity'
	plt.savefig(outdir+'fig5_wmap_planck_mfi.pdf')
	plt.clf()
	temp = mfipw[1]
	# temp[mfi_mask_512==0]=hp.UNSEEN
	hp.mollview(temp,min=-plotmax_qu,max=plotmax_qu,cmap='RdYlBu_r',title='',unit='mK CMB')#,title='MFI+WMAP+Planck weighted polarised intensity'
	plt.savefig(outdir+'fig5_wmap_planck_mfi_Q.pdf')
	plt.clf()
	temp = mfipw[2]
	# temp[mfi_mask_512==0]=hp.UNSEEN
	hp.mollview(temp,min=-plotmax_qu,max=plotmax_qu,cmap='RdYlBu_r',title='',unit='mK CMB')#,title='MFI+WMAP+Planck weighted polarised intensity'
	plt.savefig(outdir+'fig5_wmap_planck_mfi_U.pdf')
	plt.clf()

	hp.mollview(mfipw_pol,min=0,max=plotmax_p,cmap='jet',unit='mK CMB',title='',xsize=1600)#,title='MFI+WMAP+Planck weighted polarised intensity'
	plt.savefig(outdir+'fig5_wmap_planck_mfi_large.pdf')
	plt.clf()

	mfiw_polang = calc_polang(mfipw[1],mfipw[2])
	mfiw_polang[mfipw[1]==0] = hp.UNSEEN
	hp.mollview(mfiw_polang,min=-90,max=90,cmap='hsv',unit='deg',title='')
	plt.savefig(outdir+'fig5_wmap_planck_mfi_combine_polang.pdf')
	plt.clf()
	plt.close()


# EOF
