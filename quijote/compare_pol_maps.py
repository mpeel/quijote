#!/usr/bin/env python
# -*- coding: utf-8  -*-
#
# Compare the weighted polarisation maps from Planck+WMAP+QUIJOTE
# 
# Version history:
#
# 18-Feb-2020  M. Peel       Started
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt

basedir = '/Users/mpeel/Documents/maps/'
planck2018 = 'wmap9_planck2018_weight/wmap9_planck2018_tqu_10ghz_combine'
mfi = 'quijote_201907/analysetemp/mfi_combine'
planck30 = 'wmap9_planck2018_tqu/512_60.0smoothed_PlanckR3fullbeam_28.4_1024_2018_mKCMBunits.fits'
planck30old = 'wmap9_planck2015_tqu/512_60.0smoothed_PlanckR2fullbeambpcorr_28.4_256_2015_mKCMBunits.fits'
wmapk = 'wmap9_planck2018_tqu/512_60.0smoothed_wmap9beam_22.8_512_20132018_mKCMBunits.fits'

outdir=basedir+'quijote_201907/weightedanalysis/'

# Read in the planck maps
planckq = hp.read_map(basedir+planck2018+'_q.fits')
planckqerr = hp.read_map(basedir+planck2018+'_q_unc.fits')
hp.mollview(planckqerr)
plt.savefig(outdir+'Q_err_planck.pdf')
plancku = hp.read_map(basedir+planck2018+'_u.fits')
planckuerr = hp.read_map(basedir+planck2018+'_u_unc.fits')
hp.mollview(np.sqrt(planckq**2+plancku**2),min=0,max=1.2)
plt.savefig(outdir+'P_planck.pdf')

# Read in the mfi maps
mfiq = hp.read_map(basedir+mfi+'_q.fits')
mfiq[mfiq == hp.UNSEEN] = 0.0
mfiqerr = hp.read_map(basedir+mfi+'_q_unc.fits')
mfiqerr[mfiqerr==0.0] = hp.UNSEEN
hp.mollview(mfiqerr,max=0.05)
plt.savefig(outdir+'Q_err_mfi.pdf')
mfiu = hp.read_map(basedir+mfi+'_u.fits')
mfiu[mfiu == hp.UNSEEN] = 0.0
mfiuerr = hp.read_map(basedir+mfi+'_u_unc.fits')
mfiuerr[mfiuerr==0.0] = hp.UNSEEN
hp.mollview(np.sqrt(mfiq**2+mfiu**2),min=0,max=1.2)
plt.savefig(outdir+'P_mfi.pdf')

# Do a quick combination of both
combq = ((planckq/planckqerr)+(mfiq/mfiqerr))/((1.0/planckqerr)+(1.0/mfiqerr))
combu = ((plancku/planckuerr)+(mfiu/mfiuerr))/((1.0/planckuerr)+(1.0/mfiuerr))
hp.write_map(outdir+'Q_comb.fits',combq,overwrite=True)
hp.write_map(outdir+'U_comb.fits',combu,overwrite=True)
hp.write_map(outdir+'P_comb.fits',np.sqrt(combq**2+combu**2),overwrite=True)
hp.mollview(np.sqrt(combq**2+combu**2),min=0,max=1.2)
plt.savefig(outdir+'P_comb.pdf')

# Also plot just WMAPK, Planck30 on their own
planck30map = hp.read_map(basedir+planck30, field=None)
planck30oldmap = hp.read_map(basedir+planck30old, field=None)
wmapkmap = hp.read_map(basedir+wmapk, field=None)

polmapplanck30 = np.sqrt(planck30map[1]**2+planck30map[2]**2)*(10.0/28.4)**(-3.0)
hp.mollview(polmapplanck30,min=0,max=1.2,title='Planck 30 2018 (rescaled to 10GHz)')
plt.savefig(outdir+'P_planck30.pdf')

polmapplanck30old = np.sqrt(planck30oldmap[1]**2+planck30oldmap[2]**2)*(10.0/28.4)**(-3.0)
hp.mollview(polmapplanck30old,min=0,max=1.2,title='Planck 30 2015 (rescaled to 10GHz)')
plt.savefig(outdir+'P_planck30old.pdf')

polmapwmapk = np.sqrt(wmapkmap[1]**2+wmapkmap[2]**2)*(10.0/22.8)**(-3.0)
hp.mollview(polmapwmapk,min=0,max=1.2,title='WMAPK (rescaled to 10GHz)')
plt.savefig(outdir+'P_wmapk.pdf')

hp.mollview(polmapplanck30-polmapwmapk,min=-0.5,max=0.5,title='Planck 30 2018 - WMAP K (rescaled to 10GHz)')
plt.savefig(outdir+'P_planck30_minus_wmapk.pdf')

hp.mollview(polmapplanck30old-polmapwmapk,min=-0.5,max=0.5,title='Planck 30 2015 - WMAP K (rescaled to 10GHz)')
plt.savefig(outdir+'P_planck30old_minus_wmapk.pdf')

hp.mollview(polmapplanck30-polmapplanck30old,min=-0.5,max=0.5,title='Planck 30 2018 - Planck 30 2015')
plt.savefig(outdir+'P_planck30_minus_planck30.pdf')

