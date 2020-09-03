#!/usr/bin/env python
# -*- coding: utf-8  -*-
#
# Create various plots for the spurs paper
# 
# Version history:
#
# 27-Apr-2019  M. Peel       Started

import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
from astrocode.fitspectrum.astroutils import *

# Location of data
basedir = '/Users/mpeel/Documents/maps/'
planck2018 = 'wmap9_planck2018_weight/wmap9_planck2018_tqu_10ghz_combine'
planck30 = 'wmap9_planck2018_tqu/512_60.0smoothed_PlanckR3fullbeam_28.4_1024_2018_mKCMBunits.fits'
planck30old = 'wmap9_planck2015_tqu/512_60.0smoothed_PlanckR2fullbeambpcorr_28.4_256_2015_mKCMBunits.fits'
wmapk = 'wmap9_planck2018_tqu/512_60.0smoothed_wmap9beam_22.8_512_20132018_mKCMBunits.fits'

mfi_weighted = 'quijote_201907/weighted/mfi_combine.fits'
mfi_planck_weighted = 'quijote_201907/weighted/mfi_planck_10ghz_combine.fits'
mfi_maps = ['quijote_202004b/reform/512_60.0smoothed_mfi2_17.0_512_202004b_mKCMBunits.fits','quijote_202004b/reform/512_60.0smoothed_mfi2_19.0_512_202004b_mKCMBunits.fits','quijote_202004b/reform/512_60.00smoothed_mfi3_11.0_512_202004b_mKCMBunits.fits','quijote_202004b/reform/512_60.00smoothed_mfi3_13.0_512_202004b_mKCMBunits.fits','quijote_202004b/reform/512_60.0smoothed_mfi4_17.0_512_202004b_mKCMBunits.fits','quijote_202004b/reform/512_60.0smoothed_mfi4_19.0_512_202004b_mKCMBunits.fits']
mfi_freqs = [17, 19, 11, 13, 17, 19]
outdir=basedir+'quijote_202004b/spurs/'

planck_wmap_weighted = 'quijote_201907/weighted/wmap9_planck2018_10ghz_combine.fits'
cbass_map = 'cbass2019/512_60.00smoothed_cbass_4.76_512_mKCMBunits.fits'

plotmax = 1.0
plotmax_sub = 0.3
spectrum = -3.0

regions = [['California',160.60,-12.05],['Perseus',160.26,-18.62],['Rhooph',353.05,16.90]]

# Read in the MFI mask
# mfi_mask = hp.read_map(basedir+'quijote_masks/mask_quijote_ncp_satband_nside512.fits')
mfi_mask = hp.read_map(basedir+'quijote_201907/weighted/mfi_commonmask.fits')

# Weighted MFI map
mfiw = hp.read_map(basedir+mfi_weighted,field=None)
mfiw_pol = np.sqrt(mfiw[1]**2+mfiw[2]**2)
mfiw_pol[mfi_mask==0] = hp.UNSEEN
hp.mollview(mfiw_pol,min=0,max=plotmax*(10.0/mfi_freqs[2])**spectrum,cmap='jet',title='MFI weighted polarised intensity',unit='mK CMB')
plt.savefig(outdir+'mfiw.pdf')
plt.clf()

# Weighted Planck+WMAP+MFI map
mfipw = hp.read_map(basedir+mfi_planck_weighted,field=None)
mfipw_pol = np.sqrt(mfipw[1]**2+mfipw[2]**2)
mfipw_pol[mfi_mask==0] = hp.UNSEEN
hp.mollview(mfipw_pol,min=0,max=plotmax*(10.0/mfi_freqs[2])**spectrum,cmap='jet',title='MFI+WMAP+Planck weighted polarised intensity',unit='mK CMB')
plt.savefig(outdir+'wmap_planck_mfi.pdf')
plt.clf()
temp = mfipw[1]
temp[mfi_mask==0]=hp.UNSEEN
hp.mollview(temp,min=-plotmax*(10.0/mfi_freqs[2])**spectrum,max=plotmax*(10.0/mfi_freqs[2])**spectrum,cmap='RdYlBu_r',title='MFI+WMAP+Planck weighted polarised intensity',unit='mK CMB')
plt.savefig(outdir+'wmap_planck_mfi_Q.pdf')
plt.clf()
temp = mfipw[2]
temp[mfi_mask==0]=hp.UNSEEN
hp.mollview(temp,min=-plotmax*(10.0/mfi_freqs[2])**spectrum,max=plotmax*(10.0/mfi_freqs[2])**spectrum,cmap='RdYlBu_r',title='MFI+WMAP+Planck weighted polarised intensity',unit='mK CMB')
plt.savefig(outdir+'wmap_planck_mfi_U.pdf')
plt.clf()

# Weighted Planck+WMAP map
planck_wmap = hp.read_map(basedir+planck_wmap_weighted,field=None)
planck_wmap_pol = np.sqrt(planck_wmap[1]**2+planck_wmap[2]**2)
planck_wmap_pol[mfi_mask==0] = hp.UNSEEN
hp.mollview(planck_wmap_pol,min=0,max=plotmax*(10.0/mfi_freqs[2])**spectrum,cmap='jet',title='WMAP+Planck weighted polarised intensity',unit='mK CMB')
plt.savefig(outdir+'wmap_planck.pdf')
plt.clf()


for region in regions:
	reso = 5.0
	hp.gnomview(mfipw_pol,rot=[region[1],region[2]],title=region[0],reso=reso)
	hp.graticule()
	plt.savefig(outdir+region[0]+'_mfipw.pdf')
exit()

# Read in MFI
mfi311 = hp.read_map(basedir+mfi_maps[2],field=None)
mfi311_pol = np.sqrt(mfi311[1]**2+mfi311[2]**2)
mfi311_pol[mfi_mask==0] = hp.UNSEEN
hp.mollview(mfi311_pol,min=0,max=plotmax*(mfi_freqs[2]/mfi_freqs[2])**spectrum,cmap='jet',title='MFI horn 3 11GHz polarised intensity',unit='mK CMB')
plt.savefig(outdir+'mfi311.pdf')
# polfrac = mfi311_pol/mfi311[0]
# polfrac[mfi_mask==0] = hp.UNSEEN
# polfrac[polfrac < 0.0] = hp.UNSEEN
# polfrac[polfrac > 20.0] = hp.UNSEEN
# # hp.mollview(polfrac,min=0,max=0.7,cmap='jet',title='MFI horn 3 11GHz polarised fraction',unit='mK CMB')
# plt.savefig(outdir+'polfracmfi311.pdf')
# exit()
plt.clf()
mfi311_sub = (mfi311_pol * (10.0/11.0)**spectrum) - planck_wmap_pol
mfi311_sub[mfi_mask==0] = hp.UNSEEN
hp.mollview(mfi311_sub,min=-plotmax_sub,max=plotmax_sub,cmap='jet',title='MFI horn 3 11GHz polarised intensity @ 10GHz - Planck+WMAP extrapolated',unit='mK CMB')
plt.savefig(outdir+'mfi311_subP.pdf')
plt.clf()
mfi311_sub = np.sqrt(((mfi311[1] * (10.0/11.0)**spectrum) - planck_wmap[1])**2+((mfi311[2] * (10.0/11.0)**spectrum) - planck_wmap[2])**2)
mfi311_sub[mfi_mask==0] = hp.UNSEEN
hp.mollview(mfi311_sub,min=0,max=plotmax,cmap='jet',title='MFI horn 3 11GHz polarised intensity @ 10GHz - Planck+WMAP extrapolated (Q-Q, U-U)',unit='mK CMB')
plt.savefig(outdir+'mfi311_subP2.pdf')
plt.clf()

mfi311_sub = hp.sphtfunc.smoothing(mfi311_sub,fwhm=np.sqrt(2.0**2-1.0**2)*np.pi/180.0)
mfi311_sub[mfi_mask==0] = hp.UNSEEN
hp.mollview(mfi311_sub,min=0,max=plotmax,cmap='jet',title='MFI horn 3 11GHz polarised intensity @ 10GHz - Planck+WMAP extrapolated (Q-Q, U-U)',unit='mK CMB')
plt.savefig(outdir+'mfi311_subP2_smth.pdf')
plt.clf()
# exit()

mfi313 = hp.read_map(basedir+mfi_maps[3],field=None)
mfi313_pol = np.sqrt(mfi313[1]**2+mfi313[2]**2)
mfi313_pol[mfi_mask==0] = hp.UNSEEN
hp.mollview(mfi313_pol,min=0,max=plotmax*(mfi_freqs[3]/mfi_freqs[2])**spectrum,cmap='jet',title='MFI horn 3 13GHz polarised intensity',unit='mK CMB')
plt.savefig(outdir+'mfi313.pdf')
plt.clf()
mfi313_sub = (mfi313_pol * (10.0/13.0)**spectrum) - planck_wmap_pol
mfi313_sub[mfi_mask==0] = hp.UNSEEN
hp.mollview(mfi313_sub,min=-plotmax_sub,max=plotmax_sub,cmap='jet',title='MFI horn 3 13GHz polarised intensity @ 10GHz - Planck+WMAP extrapolated',unit='mK CMB')
plt.savefig(outdir+'mfi313_subP.pdf')
plt.clf()
mfi313_sub = np.sqrt(((mfi313[1] * (10.0/13.0)**spectrum) - planck_wmap[1])**2+((mfi313[2] * (10.0/13.0)**spectrum) - planck_wmap[2])**2)
mfi313_sub[mfi_mask==0] = hp.UNSEEN
hp.mollview(mfi313_sub,min=0,max=plotmax,cmap='jet',title='MFI horn 3 13GHz polarised intensity @ 10GHz - Planck+WMAP extrapolated (Q-Q, U-U)',unit='mK CMB')
plt.savefig(outdir+'mfi313_subP2.pdf')
plt.clf()

mfi217 = hp.read_map(basedir+mfi_maps[0],field=None)
mfi217_pol = np.sqrt(mfi217[1]**2+mfi217[2]**2)
mfi217_pol[mfi_mask==0] = hp.UNSEEN
hp.mollview(mfi217_pol,min=0,max=plotmax*(mfi_freqs[0]/mfi_freqs[2])**spectrum,cmap='jet',title='MFI horn 2 17GHz polarised intensity',unit='mK CMB')
plt.savefig(outdir+'mfi217.pdf')
plt.clf()
mfi219 = hp.read_map(basedir+mfi_maps[1],field=None)
mfi219_pol = np.sqrt(mfi219[1]**2+mfi219[2]**2)
mfi219_pol[mfi_mask==0] = hp.UNSEEN
hp.mollview(mfi219_pol,min=0,max=plotmax*(mfi_freqs[1]/mfi_freqs[2])**spectrum,cmap='jet',title='MFI horn 2 19GHz polarised intensity',unit='mK CMB')
plt.savefig(outdir+'mfi219.pdf')
plt.clf()
mfi417 = hp.read_map(basedir+mfi_maps[4],field=None)
mfi417_pol = np.sqrt(mfi417[1]**2+mfi417[2]**2)
mfi417_pol[mfi_mask==0] = hp.UNSEEN
hp.mollview(mfi417_pol,min=0,max=plotmax*(mfi_freqs[4]/mfi_freqs[2])**spectrum,cmap='jet',title='MFI horn 4 17GHz polarised intensity',unit='mK CMB')
plt.savefig(outdir+'mfi417.pdf')
plt.clf()

mfi417_sub = (mfi417_pol * (10.0/17.0)**spectrum) - planck_wmap_pol
mfi417_sub[mfi_mask==0] = hp.UNSEEN
hp.mollview(mfi417_sub,min=-plotmax_sub,max=plotmax_sub,cmap='jet',title='MFI horn 4 17GHz polarised intensity @ 10GHz - Planck+WMAP extrapolated',unit='mK CMB')
plt.savefig(outdir+'mfi417_subP.pdf')
plt.clf()
mfi417_sub = np.sqrt(((mfi417[1] * (10.0/17.0)**spectrum) - planck_wmap[1])**2+((mfi417[2] * (10.0/17.0)**spectrum) - planck_wmap[2])**2)
mfi417_sub[mfi_mask==0] = hp.UNSEEN
hp.mollview(mfi417_sub,min=0,max=plotmax,cmap='jet',title='MFI horn 4 17GHz polarised intensity @ 10GHz - Planck+WMAP extrapolated (Q-Q, U-U)',unit='mK CMB')
plt.savefig(outdir+'mfi417_subP2.pdf')
plt.clf()
mfi417_sub = hp.sphtfunc.smoothing(mfi417_sub,fwhm=np.sqrt(2.0**2-1.0**2)*np.pi/180.0)
mfi417_sub[mfi_mask==0] = hp.UNSEEN
hp.mollview(mfi417_sub,min=0,max=plotmax,cmap='jet',title='MFI horn 4 17GHz polarised intensity @ 10GHz - Planck+WMAP extrapolated (Q-Q, U-U)',unit='mK CMB')
plt.savefig(outdir+'mfi417_subP2_smth.pdf')
plt.clf()


mfi419 = hp.read_map(basedir+mfi_maps[5],field=None)
mfi419_pol = np.sqrt(mfi419[1]**2+mfi419[2]**2)
mfi419_pol[mfi_mask==0] = hp.UNSEEN
hp.mollview(mfi419_pol,min=0,max=plotmax*(mfi_freqs[5]/mfi_freqs[2])**spectrum,cmap='jet',title='MFI horn 4 19GHz polarised intensity',unit='mK CMB')
plt.savefig(outdir+'mfi419.pdf')
plt.clf()

mfi419_sub = (mfi419_pol * (10.0/19.0)**spectrum) - planck_wmap_pol
mfi419_sub[mfi_mask==0] = hp.UNSEEN
hp.mollview(mfi419_sub,min=-plotmax_sub,max=plotmax_sub,cmap='jet',title='MFI horn 4 19GHz polarised intensity @ 10GHz - Planck+WMAP extrapolated',unit='mK CMB')
plt.savefig(outdir+'mfi419_subP.pdf')
plt.clf()
mfi419_sub = np.sqrt(((mfi419[1] * (10.0/19.0)**spectrum) - planck_wmap[1])**2+((mfi419[2] * (10.0/19.0)**spectrum) - planck_wmap[2])**2)
mfi419_sub[mfi_mask==0] = hp.UNSEEN
hp.mollview(mfi419_sub,min=0,max=plotmax,cmap='jet',title='MFI horn 4 19GHz polarised intensity @ 10GHz - Planck+WMAP extrapolated (Q-Q, U-U)',unit='mK CMB')
plt.savefig(outdir+'mfi419_subP2.pdf')
plt.clf()



exit()

# C-BASS
cbass = hp.read_map(basedir+cbass_map,field=None)
cbass_pol = np.sqrt(cbass[1]**2+cbass[2]**2)*1e3
# cbass_pol[mfi_mask==0] = hp.UNSEEN
cbass_pol[cbass_pol>1e10] = hp.UNSEEN
hp.mollview(cbass_pol,min=0,max=plotmax*(4.76/mfi_freqs[2])**spectrum,cmap='jet',title='C-BASS polarised intensity',unit='mK CMB')
plt.savefig(outdir+'zz_cbass.pdf')
plt.clf()



# # Read in the planck maps
# planckq = hp.read_map(basedir+planck2018+'_q.fits')
# planckqerr = hp.read_map(basedir+planck2018+'_q_unc.fits')
# hp.mollview(planckqerr)
# plt.savefig(outdir+'Q_err_planck.pdf')
# plancku = hp.read_map(basedir+planck2018+'_u.fits')
# planckuerr = hp.read_map(basedir+planck2018+'_u_unc.fits')
# hp.mollview(np.sqrt(planckq**2+plancku**2),min=0,max=1.2)
# plt.savefig(outdir+'P_planck.pdf')

# # Read in the mfi maps
# mfiq = hp.read_map(basedir+mfi+'_q.fits')
# mfiq[mfiq == hp.UNSEEN] = 0.0
# mfiqerr = hp.read_map(basedir+mfi+'_q_unc.fits')
# mfiqerr[mfiqerr==0.0] = hp.UNSEEN
# hp.mollview(mfiqerr,max=0.05)
# plt.savefig(outdir+'Q_err_mfi.pdf')
# mfiu = hp.read_map(basedir+mfi+'_u.fits')
# mfiu[mfiu == hp.UNSEEN] = 0.0
# mfiuerr = hp.read_map(basedir+mfi+'_u_unc.fits')
# mfiuerr[mfiuerr==0.0] = hp.UNSEEN
# hp.mollview(np.sqrt(mfiq**2+mfiu**2),min=0,max=1.2)
# plt.savefig(outdir+'P_mfi.pdf')

# # Do a quick combination of both
# combq = ((planckq/planckqerr)+(mfiq/mfiqerr))/((1.0/planckqerr)+(1.0/mfiqerr))
# combu = ((plancku/planckuerr)+(mfiu/mfiuerr))/((1.0/planckuerr)+(1.0/mfiuerr))
# hp.write_map(outdir+'Q_comb.fits',combq,overwrite=True)
# hp.write_map(outdir+'U_comb.fits',combu,overwrite=True)
# hp.write_map(outdir+'P_comb.fits',np.sqrt(combq**2+combu**2),overwrite=True)
# hp.mollview(np.sqrt(combq**2+combu**2),min=0,max=1.2)
# plt.savefig(outdir+'P_comb.pdf')

# # Also plot just WMAPK, Planck30 on their own
# planck30map = hp.read_map(basedir+planck30, field=None)
# planck30oldmap = hp.read_map(basedir+planck30old, field=None)
# wmapkmap = hp.read_map(basedir+wmapk, field=None)

# polmapplanck30 = np.sqrt(planck30map[1]**2+planck30map[2]**2)*(10.0/28.4)**(-3.0)
# hp.mollview(polmapplanck30,min=0,max=1.2,title='Planck 30 2018 (rescaled to 10GHz)')
# plt.savefig(outdir+'P_planck30.pdf')

# polmapplanck30old = np.sqrt(planck30oldmap[1]**2+planck30oldmap[2]**2)*(10.0/28.4)**(-3.0)
# hp.mollview(polmapplanck30old,min=0,max=1.2,title='Planck 30 2015 (rescaled to 10GHz)')
# plt.savefig(outdir+'P_planck30old.pdf')

# polmapwmapk = np.sqrt(wmapkmap[1]**2+wmapkmap[2]**2)*(10.0/22.8)**(-3.0)
# hp.mollview(polmapwmapk,min=0,max=1.2,title='WMAPK (rescaled to 10GHz)')
# plt.savefig(outdir+'P_wmapk.pdf')

# hp.mollview(polmapplanck30-polmapwmapk,min=-0.5,max=0.5,title='Planck 30 2018 - WMAP K (rescaled to 10GHz)')
# plt.savefig(outdir+'P_planck30_minus_wmapk.pdf')

# hp.mollview(polmapplanck30old-polmapwmapk,min=-0.5,max=0.5,title='Planck 30 2015 - WMAP K (rescaled to 10GHz)')
# plt.savefig(outdir+'P_planck30old_minus_wmapk.pdf')

# hp.mollview(polmapplanck30-polmapplanck30old,min=-0.5,max=0.5,title='Planck 30 2018 - Planck 30 2015')
# plt.savefig(outdir+'P_planck30_minus_planck30.pdf')

