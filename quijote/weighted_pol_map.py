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
import matplotlib.colors as colors

def noiserealisation(inputmap, numpixels):
    newmap = np.zeros(numpixels)
    newmap = np.random.normal(scale=1.0, size=numpixels) * inputmap
    return newmap

nside = 512
npix = hp.nside2npix(nside)

indirectory = '/Users/mpeel/Documents/maps/quijote_201907/smooth/'
outdirectory = '/Users/mpeel/Documents/maps/quijote_201907/analyse/'
date='201907'

# indirectory = '/Users/mpeel/Documents/maps/quijote_201905/smooth/'
# outdirectory = '/Users/mpeel/Documents/maps/quijote_201905/analyse/'
# date='201905'

prefix='half1mfi'
maps_half1 = [str(nside)+'_60.0smoothed_'+prefix+'2_17.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'2_19.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'3_11.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'3_13.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'4_17.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'4_19.0_512_'+date+'_mKCMBunits.fits']#str(nside)+'_60.00smoothed_'+prefix+'1_11.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'1_13.0_512_'+date+'_mKCMBunits.fits',

prefix='half2mfi'
maps_half2 = [str(nside)+'_60.0smoothed_'+prefix+'2_17.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'2_19.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'3_11.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'3_13.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'4_17.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'4_19.0_512_'+date+'_mKCMBunits.fits']#str(nside)+'_60.00smoothed_'+prefix+'1_11.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'1_13.0_512_'+date+'_mKCMBunits.fits',

prefix='mfi'
maps = [str(nside)+'_60.0smoothed_'+prefix+'2_17.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'2_19.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'3_11.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'3_13.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'4_17.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'4_19.0_512_'+date+'_mKCMBunits.fits']#str(nside)+'_60.00smoothed_'+prefix+'1_11.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'1_13.0_512_'+date+'_mKCMBunits.fits',

doing_quijote = True
use_halfrings = False
use_weights = False
use_reweight_by_rms = True
use_reweight_by_rms_method = 2 # 1 = ricardo, 2 = alberto
# index = -2.15
index = -3.0
use_planck = True
use_cbass = False
freqs = [17,19,11,13,17,19]#11,13,
normfreq = 10.0

# These are the ones for combining Planck+WMAP
doing_quijote = True
if doing_quijote != True:
	prefix='wmap9_planck2018'
	freqs = [28.4, 44.1, 22.8, 33.0, 40.7]
	maps = ['512_60.0smoothed_PlanckR3fullbeam_28.4_1024_2018_mKCMBunits.fits','512_60.0smoothed_PlanckR3fullbeam_44.1_1024_2018_mKCMBunits.fits','512_60.0smoothed_wmap9beam_22.8_512_20132018_mKCMBunits.fits','512_60.0smoothed_wmap9beam_33.0_512_20132018_mKCMBunits.fits','512_60.0smoothed_wmap9beam_40.7_512_20132018_mKCMBunits.fits']

	prefix='wmap9'
	freqs = [22.8, 33.0, 40.7]
	maps = ['512_60.0smoothed_wmap9beam_22.8_512_20132018_mKCMBunits.fits','512_60.0smoothed_wmap9beam_33.0_512_20132018_mKCMBunits.fits','512_60.0smoothed_wmap9beam_40.7_512_20132018_mKCMBunits.fits']

	prefix='planck2018'
	freqs = [28.4, 44.1]
	maps = ['512_60.0smoothed_PlanckR3fullbeam_28.4_1024_2018_mKCMBunits.fits','512_60.0smoothed_PlanckR3fullbeam_44.1_1024_2018_mKCMBunits.fits']

	normfreq = 28.4
	indirectory = '/Users/mpeel/Documents/maps/wmap9_planck2018_tqu/'
	outdirectory = '/Users/mpeel/Documents/maps/wmap9_planck2018_weight/'

if use_planck:
	planckmap = hp.read_map('/Users/mpeel/Documents/maps/wmap_planck_pol/weighted_both_debwk_feb2015_tqu.fits',field=None)
	# planckmap = hp.pixelfunc.reorder(planckmap[0],n2r=True)
	planckmap = hp.ud_grade(planckmap[0],nside,order_in='NEST',order_out='RING')

	planck_iqu = hp.read_map('/Users/mpeel/Documents/maps/wmap9_planck2018_tqu/512_60.0smoothed_PlanckR3fullbeam_28.4_1024_2018_mKCMBunits.fits',field=None)
	# planckmap = hp.pixelfunc.reorder(planckmap[0],n2r=True)
	planck_iqu = hp.ud_grade(planck_iqu,nside,order_in='RING',order_out='RING')

if use_cbass:
	# cbass = hp.read_map('/Users/mpeel/Documents/maps/cbass2019/512_60.00smoothed_cbass_4.76_512_mKCMBunits.fits',field=None)
	cbass = hp.read_map('/Users/mpeel/Documents/maps/cbass2019/cbass_global8p8deg_swapQU_NIGHT_v28allelsNs_37_noiseCut_masked5pc_G_1024_ol500_lessTol_g_map_g_1deg_0256.fits',field=None)
	# planckmap = hp.pixelfunc.reorder(planckmap[0],n2r=True)
	cbass = hp.ud_grade(cbass,nside,order_in='RING',order_out='RING')
	newmask = np.ones(npix)
	newmask[cbass[0] < -100.0] = 0
	hp.mollview(cbass[1]*1e3*((28.4/4.76)**index)*newmask,min=-0.05,max=0.05, title='CBASS Q')
	plt.savefig(outdirectory+'cbassq.png')
	hp.mollview(cbass[2]*1e3*((28.4/4.76)**index)*newmask,min=-0.05,max=0.05,title='CBASS U')
	plt.savefig(outdirectory+'cbassu.png')
	hp.mollview(planck_iqu[1]*newmask,min=-0.05,max=0.05,title='Planck Q')
	plt.savefig(outdirectory+'planckq.png')
	hp.mollview(planck_iqu[2]*newmask,min=-0.05,max=0.05, title='Planck U')
	plt.savefig(outdirectory+'plancku.png')
	print(max(planck_iqu[1]))
	print(max(cbass[1]))
	print(max(cbass[1])*(28.4/4.76)**index)
	cbass = cbass * 1e3
	if use_planck:
		# Output a comparison between C-BASS and the Planck data
		hp.mollview(((np.sqrt(cbass[1]**2+cbass[2]**2)*(28.4/4.76)**index)-1000.0*planckmap)*newmask,min=-0.05,max=0.05,title='CBASS - (Planck+WMAP)')
		plt.savefig(outdirectory+'cbass_diff_to_planckwmap.png')
		hp.mollview((1000.0*planckmap-(np.sqrt(cbass[1]**2+cbass[2]**2)*(28.4/4.76)**index))*newmask,min=-0.05,max=0.05)
		plt.savefig(outdirectory+'cbass_diff_to_planckwmap_inverse.pdf')

		hp.mollview(((np.sqrt(cbass[1]**2+cbass[2]**2)*(28.4/4.76)**index))*newmask,min=0,max=0.05,title='CBASS P')
		plt.savefig(outdirectory+'cbass_P.png')
		hp.mollview((np.sqrt(planck_iqu[1]**2+planck_iqu[2]**2))*newmask,min=0,max=0.05,title='Planck P')
		plt.savefig(outdirectory+'planck_P.png')

		hp.mollview(((np.sqrt(cbass[1]**2+cbass[2]**2)*(28.4/4.76)**index)-np.sqrt(planck_iqu[1]**2+planck_iqu[2]**2))*newmask,min=-0.03,max=0.03,title='CBASS - Planck')
		plt.savefig(outdirectory+'cbass_diff_to_planck.png')

		hp.mollview((np.sqrt((cbass[1]*(28.4/4.76)**index - planck_iqu[1])**2+(cbass[2]*(28.4/4.76)**index - planck_iqu[2])**2))*newmask,max=0.05,title='CBASS - Planck via sqrt((QCB-QPlanck)**2 + (UCB-UPlanck)**2)')#,norm=colors.PowerNorm(gamma=0.2))
		plt.savefig(outdirectory+'cbass_diff_to_planckQU.png')
		hp.mollview((cbass[1]*(28.4/4.76)**index - planck_iqu[1])*newmask,min=-0.05,max=0.05,title='CBASS Q - Planck Q')
		plt.savefig(outdirectory+'cbass_diff_to_planckQ.png')

		hp.mollview((cbass[2]*(28.4/4.76)**index - planck_iqu[2])*newmask,min=-0.05,max=0.05,title='CBASS U - Planck U')
		plt.savefig(outdirectory+'cbass_diff_to_planckU.png')

nummaps = len(maps)
commonmask = np.ones(npix)
combine_q = np.zeros(npix)
combine_u = np.zeros(npix)
weight_q = np.zeros(npix)
weight_u = np.zeros(npix)
rescale_vals = np.zeros((3,nummaps))
for i in range(0,nummaps):
	print(maps[i])
	mapdata = hp.read_map(indirectory+maps[i],field=None)
	mapdata = hp.ud_grade(mapdata,nside,order_in='RING',order_out='RING')
	commonmask[mapdata[0][:] == hp.UNSEEN] = 0
	mapdata[0][mapdata[0][:] == hp.UNSEEN] = 0.0
	mapdata[1][mapdata[1][:] == hp.UNSEEN] = 0.0
	mapdata[2][mapdata[2][:] == hp.UNSEEN] = 0.0
	hp.mollview(mapdata[0],norm='hist')
	plt.savefig(outdirectory+maps[i]+'_0.pdf')
	hp.mollview(mapdata[1],norm='hist')
	plt.savefig(outdirectory+maps[i]+'_1.pdf')
	hp.mollview(mapdata[2],norm='hist')
	plt.savefig(outdirectory+maps[i]+'_2.pdf')
	hp.mollview(np.sqrt(mapdata[1]**2+mapdata[2]**2),min=0,max=1)
	plt.savefig(outdirectory+maps[i]+'_P.pdf')

	if doing_quijote == False:
		if 'Planck' in maps[i]:
			var_i = mapdata[4].copy()
			var_q = mapdata[7].copy()
			var_u = mapdata[9].copy()
		else:
			var_i = mapdata[3].copy()
			var_q = mapdata[3].copy()
			var_u = mapdata[3].copy()

	if use_halfrings:
		map_half1 = hp.read_map(indirectory+maps_half1[i],field=None)
		map_half2 = hp.read_map(indirectory+maps_half2[i],field=None)
		var_i = (np.abs(map_half1[0] - map_half2[0])/2.0)**2
		var_q = (np.abs(map_half1[1] - map_half2[1])/2.0)**2
		var_u = (np.abs(map_half1[2] - map_half2[2])/2.0)**2
		var_i[var_i == 0.0] = 1e4
		var_q[var_q == 0.0] = 1e4
		var_u[var_u == 0.0] = 1e4
		var_i[var_i < np.median(var_i)] = np.median(var_i)
		var_q[var_q < np.median(var_q)] = np.median(var_q)
		var_u[var_u < np.median(var_u)] = np.median(var_u)
		var_i[var_i > 1e4] = 1e4
		var_q[var_q > 1e4] = 1e4
		var_u[var_u > 1e4] = 1e4
	elif use_weights:
		# Get the variance maps
		var_i = hp.read_map(indirectory+maps[i].replace('60.0s','60.00s').replace('_mKCMBunits','_weight_0_variance'),field=None)
		var_q = hp.read_map(indirectory+maps[i].replace('60.0s','60.00s').replace('_mKCMBunits','_weight_1_variance'),field=None)
		var_u = hp.read_map(indirectory+maps[i].replace('60.0s','60.00s').replace('_mKCMBunits','_weight_2_variance'),field=None)
	elif use_reweight_by_rms:
		map_half1 = hp.read_map(indirectory+maps_half1[i],field=None)
		map_half2 = hp.read_map(indirectory+maps_half2[i],field=None)
		diff_i = (np.abs(map_half1[0] - map_half2[0])/2.0)#**2
		diff_q = (np.abs(map_half1[1] - map_half2[1])/2.0)#**2
		diff_u = (np.abs(map_half1[2] - map_half2[2])/2.0)#**2
		diff_i[diff_i > 1e4] = 0.0
		diff_q[diff_q > 1e4] = 0.0
		diff_u[diff_u > 1e4] = 0.0
		# print(np.std(diff_i[diff_i != 0.0]))
		# print(np.std(diff_q[diff_q != 0.0]))
		# print(np.std(diff_u[diff_u != 0.0]))
		var_i = hp.read_map(indirectory+maps[i].replace('60.0s','60.00s').replace('_mKCMBunits','_weight_0_variance'),field=None)
		var_q = hp.read_map(indirectory+maps[i].replace('60.0s','60.00s').replace('_mKCMBunits','_weight_1_variance'),field=None)
		var_u = hp.read_map(indirectory+maps[i].replace('60.0s','60.00s').replace('_mKCMBunits','_weight_2_variance'),field=None)
		var_i[var_i < -1e4] = 0.0
		var_q[var_q < -1e4] = 0.0
		var_u[var_u < -1e4] = 0.0
		# print(np.max(var_i))
		# print(np.min(var_i))

		if use_reweight_by_rms_method == 1:
			noise_i = noiserealisation(np.sqrt(var_i[var_i != 0.0]),len(var_i[var_i != 0.0]))
			noise_q = noiserealisation(np.sqrt(var_q[var_q != 0.0]),len(var_q[var_q != 0.0]))
			noise_u = noiserealisation(np.sqrt(var_u[var_u != 0.0]),len(var_u[var_u != 0.0]))
			rescale_vals[0,i] = np.std(diff_i[diff_i != 0.0])/np.std(noise_i[noise_i != 0.0])
			rescale_vals[1,i] = np.std(diff_q[diff_q != 0.0])/np.std(noise_q[noise_q != 0.0])
			rescale_vals[2,i] = np.std(diff_u[diff_u != 0.0])/np.std(noise_u[noise_u != 0.0])
		else:
			rescale_vals[0,i] = np.std(diff_i[var_i != 0.0] / np.sqrt(var_i[var_i != 0.0]))
			rescale_vals[1,i] = np.std(diff_q[var_q != 0.0] / np.sqrt(var_q[var_q != 0.0]))
			rescale_vals[2,i] = np.std(diff_u[var_u != 0.0] / np.sqrt(var_u[var_u != 0.0]))

		var_i[:] = var_i[:] * (rescale_vals[0,i])**2.0
		var_q[:] = var_q[:] * (rescale_vals[1,i])**2.0
		var_u[:] = var_u[:] * (rescale_vals[2,i])**2.0
		var_i[var_i == 0.0] = 1e4
		var_q[var_q == 0.0] = 1e4
		var_u[var_u == 0.0] = 1e4


	hp.mollview(var_i,norm='hist')
	plt.savefig(outdirectory+maps[i]+'_0_var.pdf')
	hp.mollview(var_q,norm='hist')
	plt.savefig(outdirectory+maps[i]+'_1_var.pdf')
	hp.mollview(var_u,norm='hist')
	plt.savefig(outdirectory+maps[i]+'_2_var.pdf')
	# print(maps[i])
	# print(np.median(np.sqrt(var_i[var_i[:] >=0])))
	# print(np.median(np.sqrt(var_q[var_q[:] >=0])))
	# print(np.median(np.sqrt(var_u[var_u[:] >=0])))

	var_i = var_i * ((normfreq/freqs[i])**index)**2
	var_q = var_q * ((normfreq/freqs[i])**index)**2
	var_u = var_u * ((normfreq/freqs[i])**index)**2

	# print(np.median(np.sqrt(var_i[var_i[:] >=0])))
	# print(np.median(np.sqrt(var_q[var_q[:] >=0])))
	# print(np.median(np.sqrt(var_u[var_u[:] >=0])))

	# Output a comparison with the Planck data
	if use_planck:
		hp.mollview(((np.sqrt(mapdata[1]**2+mapdata[2]**2)*(28.4/freqs[i])**index)-1000.0*planckmap)*commonmask,min=-0.1,max=0.1,title=maps[i] + ' - (Planck+WMAP)')
		plt.savefig(outdirectory+maps[i]+'diff_to_planckwmap.png')
		hp.mollview((1000.0*planckmap-(np.sqrt(mapdata[1]**2+mapdata[2]**2)*(28.4/freqs[i])**index))*commonmask,min=-0.1,max=0.1)
		plt.savefig(outdirectory+maps[i]+'diff_to_planckwmap_inverse.pdf')

		hp.mollview((np.sqrt((mapdata[1]*(28.4/freqs[i])**index)**2+(mapdata[2]*(28.4/freqs[i])**index)**2)-np.sqrt(planck_iqu[1]**2+planck_iqu[2]**2))*commonmask,min=-0.05,max=0.05,title=maps[i] + ' P - Planck P')
		plt.savefig(outdirectory+maps[i]+'diff_to_planckP.png')
		hp.mollview(np.sqrt((mapdata[1]*(28.4/freqs[i])**index - planck_iqu[1])**2+(mapdata[2]*(28.4/freqs[i])**index - planck_iqu[2])**2)*commonmask,min=0,max=0.1,title=maps[i] + ' - Planck via sqrt((QMFI-QPlanck)**2 + (UMFI-UPlanck)**2)')#,norm=colors.PowerNorm(gamma=0.2))
		plt.savefig(outdirectory+maps[i]+'diff_to_planckQU.png')
		hp.mollview((mapdata[1]*(28.4/freqs[i])**index - planck_iqu[1])*commonmask,min=-0.05,max=0.05,title=maps[i] + ' Q - Planck Q')
		plt.savefig(outdirectory+maps[i]+'diff_to_planckQ.png')

		hp.mollview((mapdata[2]*(28.4/freqs[i])**index - planck_iqu[2])*commonmask,min=-0.05,max=0.05,title=maps[i] + ' U - Planck U')
		plt.savefig(outdirectory+maps[i]+'diff_to_planckU.png')

		hp.gnomview(((mapdata[0]*(28.4/freqs[i])**index))*commonmask,title=maps[i] + ' I (rescaled to 28.4)',rot=[80,0],max=10.0,reso=10.0)
		plt.savefig(outdirectory+maps[i]+'_I_cyg.png')
		hp.gnomview(((mapdata[0]*(28.4/freqs[i])**index)-planck_iqu[0])*commonmask,min=-2.0,max=2.0,title=maps[i] + ' I - Planck I',rot=[80,0],reso=10.0)
		plt.savefig(outdirectory+maps[i]+'diff_to_planckI_cyg.png')
		hp.gnomview((np.sqrt((mapdata[1]*(28.4/freqs[i])**index)**2+(mapdata[2]*(28.4/freqs[i])**index)**2)-np.sqrt(planck_iqu[1]**2+planck_iqu[2]**2))*commonmask,min=-0.05,max=0.05,rot=[80,0],reso=10.0,title=maps[i] + ' P - Planck P')
		plt.savefig(outdirectory+maps[i]+'diff_to_planckP_cyg.png')
		hp.gnomview(np.sqrt((mapdata[1]*(28.4/freqs[i])**index - planck_iqu[1])**2+(mapdata[2]*(28.4/freqs[i])**index - planck_iqu[2])**2)*commonmask,min=0,max=0.1,rot=[80,0],reso=10.0,title=maps[i] + ' - Planck via sqrt((QMFI-QPlanck)**2 + (UMFI-UPlanck)**2)')#,norm=colors.PowerNorm(gamma=0.2))
		plt.savefig(outdirectory+maps[i]+'diff_to_planckQU_cyg.png')
		hp.gnomview((mapdata[1]*(28.4/freqs[i])**index - planck_iqu[1])*commonmask,min=-0.1,max=0.1,rot=[80,0],reso=10.0,title=maps[i] + ' Q - Planck Q')
		plt.savefig(outdirectory+maps[i]+'diff_to_planckQ_cyg.png')

		hp.gnomview((mapdata[2]*(28.4/freqs[i])**index - planck_iqu[2])*commonmask,min=-0.1,max=0.1,rot=[80,0],reso=10.0,title=maps[i] + ' U - Planck U')
		plt.savefig(outdirectory+maps[i]+'diff_to_planckU_cyg.png')
		# exit()



	if i == 0:
		combine_q = (mapdata[1].copy()*(normfreq/freqs[i])**index)/var_q
		combine_u = (mapdata[2].copy()*(normfreq/freqs[i])**index)/var_u
		weight_q = 1.0/(var_q.copy())
		weight_u = 1.0/(var_u.copy())
	else:
		combine_q = combine_q+(mapdata[1]*(normfreq/freqs[i])**index)/var_q
		combine_u = combine_u+(mapdata[2]*(normfreq/freqs[i])**index)/var_u
		weight_q = weight_q+1.0/var_q
		weight_u = weight_u+1.0/var_u

combine_q /= weight_q
combine_u /= weight_u

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
hp.mollview(np.sqrt(combine_q**2+combine_u**2)*commonmask,min=0,max=1.2,cmap=plt.get_cmap('jet'))
plt.savefig(outdirectory+prefix+'_combine_P.pdf')

hp.write_map(outdirectory+prefix+'_combine_P_nomask.fits',np.sqrt(combine_q**2+combine_u**2),overwrite=True)
hp.mollview(np.sqrt(combine_q**2+combine_u**2),min=0,max=1.2,cmap=plt.get_cmap('jet'))
plt.savefig(outdirectory+prefix+'_combine_P_nomask.pdf')

hp.write_map(outdirectory+prefix+'_combine.fits',[np.sqrt(combine_q**2+combine_u**2),combine_q,combine_u,1.0/weight_q,1.0/weight_u],overwrite=True)

# commonmask2 = hp.ud_grade(commonmask,256,order_in='RING',order_out='RING')
if use_planck:
	hp.write_map(outdirectory+'wmapplanck2015.fits',[np.sqrt(planck_iqu[1]**2+planck_iqu[2]**2),planck_iqu[1],planck_iqu[2]],overwrite=True)
	hp.mollview(planckmap*1000.0*commonmask,min=0,max=0.03,cmap=plt.get_cmap('jet'))
	plt.savefig(outdirectory+'combine_P_planck.pdf')
	hp.mollview(planckmap*1000.0,min=0,max=0.06,cmap=plt.get_cmap('jet'))
	plt.savefig(outdirectory+'combine_P_planck_nomask.pdf')

np.set_printoptions(formatter={'float': '{: 0.5f}'.format})
print(rescale_vals)

