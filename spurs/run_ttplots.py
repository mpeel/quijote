#!/usr/bin/env python
# -*- coding: utf-8  -*-
#
# Make tt plots for the Spurs paper
#
# Version history:
#
# 09-Jun-2022  M. Peel       Started
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
from astrocode.astroutils import *
import matplotlib
from configure import *
from astrocode.ttplot import *
from astrocode.ttplot_fuskeland import *

nside = 64
npix = 12*nside*nside

# Location of data
basedir = '/Users/mpeel/Documents/maps/'
# outdir=basedir+'quijote_202103_tqu_v1.5_noise_v1.0_newwf/spurs/'
outdir=basedir+'2026_loops_spurs/'

## Planck/WMAP
plancknpipe30 = 'planck2020_tqu_v1.5_noise_v1.0_10k/'+str(nside)+'_60.0smoothed_PlanckR4fullbeamnodpNoise_28.4_1024_2020_mKCMBunits.fits'
wmapk9 = 'wmap9_tqu_v1.5_noise_v1.0_10k/'+str(nside)+'_60.0smoothed_wmap9beamNoise_22.8_512_2013_mKCMBunits.fits'
planck_wmap_weighted = 'quijote_202103_tqu_v1.5_noise_v1.0_weighted/wmapplanck_combine.fits'

## MFI
mfi_weighted = 'quijote_202103_tqu_v1.5_noise_v1.0_weighted/'+str(nside)+'_60.0smoothed_quijotecombwei10_tqu_v1.5_noise_v1.0_-3.0_combine.fits'
mfi_planck_weighted = 'quijote_202103_tqu_v1.5_noise_v1.0_weighted/'+str(nside)+'_60.0smoothed_quijotecombwei10_tqu_v1.5_noise_v1.0_-3.0_wmapplanck_combine.fits'
mfi_maps = ['quijote_202103_tqu_v1.5_noise_v1.0_newwf/'+str(nside)+'_60.0smoothed_QUIJOTEMFI2_17.0_2021_mKCMBunits.fits',\
'quijote_202103_tqu_v1.5_noise_v1.0_newwf/'+str(nside)+'_60.0smoothed_QUIJOTEMFI2_19.0_2021_mKCMBunits.fits',\
'quijote_202103_tqu_v1.5_noise_v1.0_newwf/'+str(nside)+'_60.0smoothed_QUIJOTEMFI3_11.0_2021_mKCMBunits.fits',\
'quijote_202103_tqu_v1.5_noise_v1.0_newwf/'+str(nside)+'_60.0smoothed_QUIJOTEMFI3_13.0_2021_mKCMBunits.fits',\
'quijote_202103_tqu_v1.5_noise_v1.0_newwf/'+str(nside)+'_60.0smoothed_QUIJOTEMFI4_17.0_2021_mKCMBunits.fits',\
'quijote_202103_tqu_v1.5_noise_v1.0_newwf/'+str(nside)+'_60.0smoothed_QUIJOTEMFI4_19.0_2021_mKCMBunits.fits']
mfi_mask_filename = 'quijote_masks/mask_quijote_ncp_satband_nside512.fits'
# mfi_mask_filename = basedir+'quijote_201907/weighted/mfi_commonmask.fits'

regions = get_regions()

mfi_mask = hp.read_map(basedir+mfi_mask_filename)
mfi_mask = hp.ud_grade(mfi_mask,nside_out=nside)
mfi_mask[np.where(mfi_mask>0.5)] = 1.0
mfi_mask[np.where(mfi_mask<=0.5)] = 0.0
hp.mollview(mfi_mask)
plt.savefig(outdir+'mfi_mask.png')
plt.clf()
plt.close()

mfi311 = hp.read_map(basedir+mfi_maps[2],field=None)
p30 = hp.read_map(basedir+plancknpipe30,field=None)
wk9 = hp.read_map(basedir+wmapk9,field=None)

fontsize = 18
matplotlib.rcParams.update({'font.size':fontsize})


def run_plots(mask,mfi311,p30,wk9,outdir,name=''):

	plt.clf()
	plt.close()

	hp.mollview(mask)
	plt.savefig(outdir+name+'_mask.png')
	plt.clf()
	plt.close()

	tempmap = np.sqrt(mfi311[1]**2+mfi311[2]**2)
	tempmap[np.where(mask==0)]=0.0
	hp.mollview(tempmap,max=0.5,unit='mK CMB',title='')
	plt.savefig(outdir+name+'_q311pol.png')
	tempmap = np.sqrt(p30[1]**2+p30[2]**2)
	tempmap[np.where(mask==0)]=0.0
	hp.mollview(tempmap,max=0.05,unit='mK CMB',title='')
	plt.savefig(outdir+name+'_p30pol.png')
	tempmap = p30[1]
	tempmap[np.where(mask==0)]=0.0
	hp.mollview(tempmap,max=0.05,unit='mK CMB',title='')
	plt.savefig(outdir+name+'_p30pol_q.png')
	tempmap = p30[2]
	tempmap[np.where(mask==0)]=0.0
	hp.mollview(tempmap,max=0.05,unit='mK CMB',title='')
	plt.savefig(outdir+name+'_p30pol_u.png')
	plt.clf()
	plt.close()

	plot_tt(mfi311[1][np.where(mask==1)],p30[1][np.where(mask==1)],outdir+name+'_tt_311_30_q.png',freq1=11.0,freq2=28.4,xlabel='MFI 311 Q',ylabel='Planck N30 Q',sigma_x=np.sqrt(mfi311[4][np.where(mask==1)]),sigma=np.sqrt(p30[4][np.where(mask==1)]))
	plot_tt(mfi311[2][np.where(mask==1)],p30[2][np.where(mask==1)],outdir+name+'_tt_311_30_u.png',freq1=11.0,freq2=28.4,xlabel='MFI 311 U',ylabel='Planck N30 U',sigma_x=np.sqrt(mfi311[6][np.where(mask==1)]),sigma=np.sqrt(p30[6][np.where(mask==1)]))
	plot_tt(mfi311[1][np.where(mask==1)],wk9[1][np.where(mask==1)],outdir+name+'_tt_311_K9_q.png',freq1=11.0,freq2=22.8,xlabel='MFI 311 Q',ylabel='WMAP K9 Q',sigma_x=np.sqrt(mfi311[4][np.where(mask==1)]),sigma=np.sqrt(wk9[4][np.where(mask==1)]))
	plot_tt(mfi311[2][np.where(mask==1)],wk9[2][np.where(mask==1)],outdir+name+'_tt_311_K9_u.png',freq1=11.0,freq2=22.8,xlabel='MFI 311 U',ylabel='WMAP K9 U',sigma_x=np.sqrt(mfi311[6][np.where(mask==1)]),sigma=np.sqrt(wk9[6][np.where(mask==1)]))
	plt.clf()
	plt.close()

	beta, sigma_beta, sigma_beta2, sigma_beta3, q, sigmaq, chi = p_fusk14_cc(mask, mfi311, mfi311[4:], 11.1, 'Q11', False, p30, p30[4:], 28.4, 'P30', False, 0.0, outdir, name)


# # # NPS
# mask = hp.read_map('/Users/mpeel/Documents/maps/quijote_masks/bob_R13_5reg_mfi_ns64.fits')
# mask[mask > 1] = 1
# run_plots(mask,mfi311,p30.copy(),wk9,outdir,name='nps')

# # NPS diffuse nearby
# mask = hp.read_map('/Users/mpeel/Documents/maps/quijote_masks/bob_R13_2reg_mfi_ns64.fits')
# mask[mask==1] = 0
# mask[mask==2] = 1
# run_plots(mask,mfi311,p30.copy(),wk9,outdir,name='nps_diffuse')

# # NPS diffuse nearby minus CGS
# mask = hp.read_map('/Users/mpeel/Documents/maps/quijote_masks/bob_R13_2reg_mfi_ns64.fits')
# mask[mask==1] = 0
# mask[mask==2] = 1
# mask2 = hp.read_map('/Users/mpeel/Documents/maps/quijote_masks/bob_CGS_msk.fits')
# mask[mask2==1] = 0
# run_plots(mask,mfi311,p30.copy(),wk9,outdir,name='nps_diffuse_minus_cgs')

# # CGS diffuse nearby
# mask = hp.read_map('/Users/mpeel/Documents/maps/quijote_masks/bob_CGS_msk.fits')
# run_plots(mask,mfi311,p30.copy(),wk9,outdir,name='cgs')

# # Fan
# mask = healpixmask(nside, 100.0, 180.0, -15, 15)
# ps = query_ellipse(nside, 111.8, -2.4, 2.0, 1.0, 0.0)
# mask[ps] = 0.0
# ps = query_ellipse(nside, 133, 1.9, 2.0, 1.0, 0.0)
# mask[ps] = 0.0
# run_plots(mask,mfi311,p30.copy(),wk9,outdir,name='fan')

# Loop 3
mask = np.zeros(npix)
# Make an annulus mask
angle=30.0*np.pi/180.0
abratio=0.6
outer = query_ellipse(nside, 115, 27, 40, abratio, angle)
inner = query_ellipse(nside, 115, 27, 20, abratio, angle)
mask[outer]=1.0
mask[inner]=0.0
# Mask where we have no data
mask[np.where(mfi311[1] == hp.UNSEEN)] = 0.0
# Mask the plane
for i in range(0,npix):
	pos = hp.pixelfunc.pix2ang(nside, i)
	if np.abs(((90.0-(pos[0]*180.0/pi)-0.3*(pos[1]*180.0/pi))) <= -15.0):
		mask[i] = 0
run_plots(mask,mfi311,p30.copy(),wk9,outdir,name='loop3')

# mask = hp.read_map('/Users/mpeel/Documents/maps/quijote_masks/regionx.fits')
# mask = hp.ud_grade(mask,nside_out=nside)
# run_plots(mask,mfi311,p30.copy(),wk9,outdir,name='regionx')

# mask = hp.read_map('/Users/mpeel/Documents/maps/quijote_masks/SPS_mask.fits')
# mask = hp.ud_grade(mask,nside_out=nside)
# mask = mask * mfi_mask
# run_plots(mask,mfi311,p30.copy(),wk9,outdir,name='sps')

# mask = hp.read_map('/Users/mpeel/Documents/maps/quijote_masks/Loop3_south_mask.fits')
# mask = hp.ud_grade(mask,nside_out=nside)
# mask = mask * mfi_mask
# run_plots(mask,mfi311,p30.copy(),wk9,outdir,name='loop3south')

# mask = hp.read_map('/Users/mpeel/Documents/maps/quijote_masks/Loop3_shell_mask.fits')
# mask = hp.ud_grade(mask,nside_out=nside)
# mask = mask * mfi_mask
# run_plots(mask,mfi311,p30.copy(),wk9,outdir,name='loop3shell')

# mask = hp.read_map('/Users/mpeel/Documents/maps/quijote_masks/loop3_shell_correct.fits')
# mask = hp.ud_grade(mask,nside_out=nside)
# mask = mask * mfi_mask
# run_plots(mask,mfi311,p30.copy(),wk9,outdir,name='loop3_shell')

# mask = hp.read_map('/Users/mpeel/Documents/maps/quijote_masks/shell_south_of_cygnus.fits')
# mask = hp.ud_grade(mask,nside_out=nside)
# mask = mask * mfi_mask
# run_plots(mask,mfi311,p30.copy(),wk9,outdir,name='shell_south_of_cygnus')

# mask = hp.read_map('/Users/mpeel/Documents/maps/quijote_masks/loop_near_Halpha_filament.fits')
# mask = hp.ud_grade(mask,nside_out=nside)
# mask = mask * mfi_mask
# run_plots(mask,mfi311,p30.copy(),wk9,outdir,name='loop_near_Halpha_filament')

# mask = hp.read_map('/Users/mpeel/Documents/maps/quijote_masks/Loop11_mask.fits')
# mask = hp.ud_grade(mask,nside_out=nside)
# mask = mask * mfi_mask
# run_plots(mask,mfi311,p30.copy(),wk9,outdir,name='loop11')

# mask = hp.read_map('/Users/mpeel/Documents/maps/quijote_masks/LoopLB_mask.fits')
# mask = hp.ud_grade(mask,nside_out=nside)
# mask = mask * mfi_mask
# run_plots(mask,mfi311,p30.copy(),wk9,outdir,name='looplb')

# mask = hp.read_map('/Users/mpeel/Documents/maps/quijote_masks/loop_IV_mask.fits')
# mask = hp.ud_grade(mask,nside_out=nside)
# mask = mask * mfi_mask
# run_plots(mask,mfi311,p30.copy(),wk9,outdir,name='loopiv')

# mask = hp.read_map('/Users/mpeel/Documents/maps/quijote_masks/loop9_only.fits')
# mask = hp.ud_grade(mask,nside_out=nside)
# mask = mask * mfi_mask
# run_plots(mask,mfi311,p30.copy(),wk9,outdir,name='loop9')

# mask = hp.read_map('/Users/mpeel/Documents/maps/quijote_masks/loop_between_9_and_GCS_only.fits')
# mask = hp.ud_grade(mask,nside_out=nside)
# mask = mask * mfi_mask
# run_plots(mask,mfi311,p30.copy(),wk9,outdir,name='loop_between_9_and_cgs')

# mask = hp.read_map('/Users/mpeel/Documents/maps/quijote_masks/loop_GCS_only.fits')
# mask = hp.ud_grade(mask,nside_out=nside)
# mask = mask * mfi_mask
# run_plots(mask,mfi311,p30.copy(),wk9,outdir,name='cgs')

# mask = hp.read_map('/Users/mpeel/Documents/maps/quijote_masks/loop_near_halpha_filament (1).fits')
# mask = hp.ud_grade(mask,nside_out=nside)
# mask = mask * mfi_mask
# run_plots(mask,mfi311,p30.copy(),wk9,outdir,name='loop_near_halpha_2')

# mask = hp.read_map('/Users/mpeel/Documents/maps/quijote_masks/loop1_only.fits')
# mask = hp.ud_grade(mask,nside_out=nside)
# mask = mask * mfi_mask
# run_plots(mask,mfi311,p30.copy(),wk9,outdir,name='loop1')

# mask = hp.read_map('/Users/mpeel/Documents/maps/quijote_masks/loop_viib.fits')
# mask = hp.ud_grade(mask,nside_out=nside)
# mask = mask * mfi_mask
# run_plots(mask,mfi311,p30.copy(),wk9,outdir,name='loop_viib')


mask = hp.read_map('/Users/mpeel/Documents/maps/quijote_masks/loop_2.fits')
mask = hp.ud_grade(mask,nside_out=nside)
mask = mask * mfi_mask
run_plots(mask,mfi311,p30.copy(),wk9,outdir,name='loop_2')


