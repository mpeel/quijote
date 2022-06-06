#!/usr/bin/env python
# -*- coding: utf-8  -*-
#
# Do a quick analysis of the polarisation of the MFI maps
#
# Version history:
#
# 03-Jun-2019  M. Peel       Started

import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from scipy import odr
from astrocode.astroutils import *
from astrocode.ttplot import *
from astrocode.polfunc import *

def run_offset_map(map1,map2,sigmamap,sigmamap_x=[],nsidemask=[],nside=8,outputtt='temp',outputmap='temp',sigmacut=3.0):
	offsetmap = np.zeros(hp.nside2npix(nside))
	uncmap = np.zeros(hp.nside2npix(nside))
	for i in range(0,int(np.max(nsidemask))):
		if len(map1[nsidemask==i]) > 3:
			if sigmamap_x != []:
				print(sigmamap_x[nside==i])
				try:
					fit,fiterr=plot_tt(map1[nsidemask==i],map2[nsidemask==i],outputtt+str(i)+'.png',sigma=sigmamap[nsidemask==i],sigma_x=sigmamap_x[nside==i])
				except:
					fit = [0.0, hp.UNSEEN]
					fiterr = [0.0, 0.0]
			else:
				try:
					fit,fiterr=plot_tt(map1[nsidemask==i],map2[nsidemask==i],outputtt+str(i)+'.png',sigma=sigmamap[nsidemask==i])
				except:
					fit = [0.0, hp.UNSEEN]
					fiterr = [0.0, 0.0]
			offsetmap[i] = fit[1]
			uncmap[i] = fiterr[1]
	avgerr = np.nanmedian(uncmap[uncmap != 0.0])
	# print(avgerr)
	# input('Continue?')
	offsetmap[np.abs(uncmap) > sigmacut * avgerr] = hp.UNSEEN
	uncmap[np.abs(uncmap) > sigmacut * avgerr] = hp.UNSEEN
	offsetmap[uncmap==0] = hp.UNSEEN
	uncmap[uncmap==0] = hp.UNSEEN
	if outputmap != 'temp':
		hp.write_map(outputmap+'.fits',[offsetmap,uncmap],overwrite=True)
		hp.mollview(offsetmap)
		plt.savefig(outputmap)
		plt.clf()
		plt.close()
		hp.mollview(uncmap)
		plt.savefig(outputmap.replace('.','_unc.'))
		plt.clf()
		plt.close()

	return offsetmap

def compare_polang(prefix='mfi', date='201905',datestr='may2019',use_variance=True,indirectory='',newformat=False,use_polang_err=True, polang_err_threshold = 5.0,planckmap='npipe',applyoffsets='wmap',outdirectory='',mapprefix='',inputmask='/Users/mpeel/Documents/maps/quijote_masks/mask_quijote_ncp_lowdec_nside512.fits',staticmask=False,outextraext='',simnoise=False):
	nside = 512
	npix = hp.nside2npix(nside)
	nside_out=64
	npix_out=hp.nside2npix(nside_out)

	# Map location
	if newformat:
		if simnoise:
			maps = ['64_60.0smoothed_QUIJOTEMFI1'+mapprefix+'_11.0_2021simnoise_mKCMBunits.fits','64_60.0smoothed_QUIJOTEMFI1'+mapprefix+'_13.0_2021simnoise_mKCMBunits.fits','64_60.0smoothed_QUIJOTEMFI2'+mapprefix+'_17.0_2021simnoise_mKCMBunits.fits','64_60.0smoothed_QUIJOTEMFI2'+mapprefix+'_19.0_2021simnoise_mKCMBunits.fits','64_60.0smoothed_QUIJOTEMFI3'+mapprefix+'_11.0_2021simnoise_mKCMBunits.fits','64_60.0smoothed_QUIJOTEMFI3'+mapprefix+'_13.0_2021simnoise_mKCMBunits.fits','64_60.0smoothed_QUIJOTEMFI4'+mapprefix+'_17.0_2021simnoise_mKCMBunits.fits','64_60.0smoothed_QUIJOTEMFI4'+mapprefix+'_19.0_2021simnoise_mKCMBunits.fits']
		else:
			maps = ['64_60.0smoothed_QUIJOTEMFI1'+mapprefix+'_11.0_2021_mKCMBunits.fits','64_60.0smoothed_QUIJOTEMFI1'+mapprefix+'_13.0_2021_mKCMBunits.fits','64_60.0smoothed_QUIJOTEMFI2'+mapprefix+'_17.0_2021_mKCMBunits.fits','64_60.0smoothed_QUIJOTEMFI2'+mapprefix+'_19.0_2021_mKCMBunits.fits','64_60.0smoothed_QUIJOTEMFI3'+mapprefix+'_11.0_2021_mKCMBunits.fits','64_60.0smoothed_QUIJOTEMFI3'+mapprefix+'_13.0_2021_mKCMBunits.fits','64_60.0smoothed_QUIJOTEMFI4'+mapprefix+'_17.0_2021_mKCMBunits.fits','64_60.0smoothed_QUIJOTEMFI4'+mapprefix+'_19.0_2021_mKCMBunits.fits']
		varmaps = ['', '', '', '', '', '', '', '']
	else:
		maps = [prefix+'_'+datestr+'_mapsmth_11.0_1.fits',prefix+'_'+datestr+'_mapsmth_13.0_1.fits',prefix+'_'+datestr+'_mapsmth_17.0_2.fits',prefix+'_'+datestr+'_mapsmth_19.0_2.fits',prefix+'_'+datestr+'_mapsmth_11.0_3.fits',prefix+'_'+datestr+'_mapsmth_13.0_3.fits',prefix+'_'+datestr+'_mapsmth_17.0_4.fits',prefix+'_'+datestr+'_mapsmth_19.0_4.fits']
		varmaps = ['', '', 'std_H2_17_sm1deg_nside64.fits', 'std_H2_19_sm1deg_nside64.fits', 'std_H3_11_sm1deg_nside64.fits', 'std_H3_13_sm1deg_nside64.fits', 'std_H4_17_sm1deg_nside64.fits', 'std_H4_19_sm1deg_nside64.fits']
	print(maps)
	exit()
	# Directories
	if indirectory == '':
		indirectory = '/Users/mpeel/Documents/maps/quijote_'+date+'/reform/'
		outdirectory = '/Users/mpeel/Documents/maps/quijote_'+date+'/analyse/'
	if outdirectory == '':
		outdirectory = indirectory+'/analyse'+mapprefix+outextraext
		if not newformat:
			indirectory = indirectory + '/reform/'
	if use_polang_err == False:
		outdirectory += '_all/'
	else:
		outdirectory += '_'+str(polang_err_threshold)+'/'
	print(outdirectory)
	ensure_dir(outdirectory)

	logfile = open(outdirectory+"_results.txt", "w")
	logfile.write(outdirectory+"\n")
	nsidemask = np.asarray(nside_mask(nside_out,8))
	hp.mollview(nsidemask,cmap=plt.get_cmap('jet'))
	plt.savefig(outdirectory+'_nside_mask.png')
	plt.clf()
	plt.close()

	nummaps = len(maps)
	freqs = [11,13,16.7,18.7,11.1,12.9,17,19]
	instrument = ['MFI1','MFI1','MFI2','MFI2','MFI3','MFI3','MFI4','MFI4']
	index = -3.0
	wmap_freq = 22.8
	planck_freq = 28.4
	print(inputmask)
	commonmask = hp.read_map(inputmask,field=None)
	# commonmask = hp.read_map('/Users/mpeel/Documents/maps/quijote_masks/mask_quijote_ncp_satband_nside512.fits',field=None)
	# commonmask = hp.read_map(outdirectory+'mfi_commonmask.fits',field=None)
	commonmask = hp.ud_grade(commonmask,nside_out,order_in='RING',order_out='RING')

	# Make a quick mask to get rid of low signal-to-noise pixels
	use_threshold = False
	threshold = 0.8
	threshold2 = 2.5
	use_sn = False
	sn_ratio = 5.0

	i = 4 # Use the 11GHz map for this
	mapdata = hp.read_map(indirectory+maps[i],field=None)
	qmap_311 = hp.ud_grade(mapdata[1],nside_out,order_in='RING',order_out='RING')*commonmask
	umap_311 = hp.ud_grade(mapdata[2],nside_out,order_in='RING',order_out='RING')*commonmask
	if use_variance:
		if 'nov2019' in datestr or 'mar2020' in datestr or 'apr2020' in datestr or 'nov2020' in datestr:
			if 'FG_QJT' in maps[i]:
				var = hp.read_map(indirectory+maps[i].replace('FG_QJT','FG_wei_QJT'),field=None)
				# print(np.shape(var))
				var_q_311 = 1.0/var[1]
				var_u_311 = 1.0/var[2]
				var_q_311 = hp.ud_grade(var_q_311, nside_out, power=2)
				var_u_311 = hp.ud_grade(var_u_311, nside_out, power=2)
			elif varmaps[i] != '':
				var = hp.read_map('/Volumes/Toshiba5TB2/maps/quijote_201911/noise/'+varmaps[i],field=None)
				var_q_311 = var[1].copy()**2
				var_u_311 = var[2].copy()**2
			else:
				var_q_311 = np.ones(hp.pixelfunc.nside2npix(64))
				var_u_311 = np.ones(hp.pixelfunc.nside2npix(64))
		else:
			if newformat:
				var_q_311 = mapdata[4]
				var_u_311 = mapdata[6]
			else:
				var_q_311 = hp.read_map(indirectory+maps[i].replace('60.0s','60.00s').replace('_mKCMBunits','_weight_1_variance'),field=None)
				var_u_311 = hp.read_map(indirectory+maps[i].replace('60.0s','60.00s').replace('_mKCMBunits','_weight_2_variance'),field=None)
			# var_q_311 = hp.read_map(indirectory+maps[i].replace('60.0s','60.00s').replace('_mKCMBunits','_weight_1_variance'),field=None)
			# var_u_311 = hp.read_map(indirectory+maps[i].replace('60.0s','60.00s').replace('_mKCMBunits','_weight_2_variance'),field=None)
			var_q_311 = hp.ud_grade(var_q_311, nside_out, power=0)
			var_u_311 = hp.ud_grade(var_u_311, nside_out, power=0)

		var_q_311[var_q_311 < 0.0] = 10000.0
		var_u_311[var_u_311 < 0.0] = 10000.0
		var_q_311[var_q_311==0.0] = 10000.0
		var_u_311[var_u_311==0.0] = 10000.0
		var_q_311[~np.isfinite(var_q_311)] = 10000.0
		var_u_311[~np.isfinite(var_u_311)] = 10000.0


	polang_map_11 = calc_polang(qmap_311[:],umap_311[:])
	hp.mollview(polang_map_11,cmap=plt.get_cmap('jet'))
	plt.savefig(outdirectory+'polang_map_mfi11.png')
	plt.clf()
	plt.close()
	if use_variance:
		polang_unc_map_11 = calc_polang_unc(qmap_311,umap_311,np.sqrt(var_q_311),np.sqrt(var_u_311))
		hp.mollview(polang_unc_map_11,cmap=plt.get_cmap('jet'),max=polang_err_threshold)
		plt.savefig(outdirectory+'polang_map_mfi11_unc.png')
		tempmap = polang_map_11.copy()
		tempmap[polang_unc_map_11 > polang_err_threshold] = 0.0
		hp.mollview(tempmap,cmap=plt.get_cmap('jet'))
		plt.savefig(outdirectory+'polang_map_mfi11_threshold.png')
		plt.clf()
		plt.close()

	# define the mask
	if staticmask:
		quickmask = commonmask.copy()
	else:
		quickmask = np.zeros(npix_out)
		if use_threshold:
			quickmask[np.sqrt(qmap_311**2+umap_311**2) > threshold] = 1
			quickmask[np.sqrt(qmap_311**2+umap_311**2) > threshold2] = 0
		elif use_sn:
			quickmask[np.sqrt((qmap_311/np.sqrt(var_q_311))**2+(umap_311/np.sqrt(var_u_311))**2) > sn_ratio] = 1
		elif use_polang_err:
			quickmask[polang_unc_map_11 < polang_err_threshold] = 1
		else:
			quickmask[commonmask==1] = 1
		quickmask[commonmask == 0] = 0

	# othermask = hp.read_map('/Users/mpeel/Documents/maps/quijote_masks/mask_quijote_ncp_lowdec_nside512.fits')
	# othermask = hp.ud_grade(othermask,nside_out,order_in='RING',order_out='RING')
	# othermask[othermask != 1.0] = 0
	# quickmask = othermask
	# nsidemask = np.multiply(nsidemask,othermask)
	hp.mollview(nsidemask,cmap=plt.get_cmap('jet'))
	plt.savefig(outdirectory+'_nside_mask2.png')
	plt.clf()
	plt.close()



	# Combined Planck+WMAP
	# mapdata = hp.read_map('/Users/mpeel/Documents/maps/wmap_planck_pol/weighted_both_debwk_feb2015_tqu.fits',field=None)
	# mapdata1 = hp.pixelfunc.reorder(-mapdata[1],n2r=True)
	# mapdata2 = hp.pixelfunc.reorder(-mapdata[2],n2r=True)
	if planckmap == '2015':
		mapdata = hp.read_map('/Users/mpeel/Documents/maps/weighted_wmap_planck_qu_tqu_v1.5_noise_v1.0/64_60.0smoothed_wmap9_planck2015_tqu_v1.5_noise_v1.0_-3.0_combine.fits',field=None)
	elif planckmap == '2018':
		mapdata = hp.read_map('/Users/mpeel/Documents/maps/weighted_wmap_planck_qu_tqu_v1.5_noise_v1.0/64_60.0smoothed_wmap9_planck2018_tqu_v1.5_noise_v1.0_-3.0_combine.fits',field=None)
	else:
		# npipe
		mapdata = hp.read_map('/Users/mpeel/Documents/maps/weighted_wmap_planck_qu_tqu_v1.5_noise_v1.0/64_60.0smoothed_wmap9_planck2020_tqu_v1.5_noise_v1.0_-3.0_combine.fits',field=None)
	mapdata1 = mapdata[1]
	mapdata2 = mapdata[2]
	qmap = hp.ud_grade(mapdata1,nside_out,order_in='RING',order_out='RING')*commonmask
	umap = hp.ud_grade(mapdata2,nside_out,order_in='RING',order_out='RING')*commonmask
	qmap_var = hp.ud_grade(mapdata[3],nside_out,order_in='RING',order_out='RING')*commonmask
	umap_var = hp.ud_grade(mapdata[5],nside_out,order_in='RING',order_out='RING')*commonmask
	hp.mollview(qmap,norm='hist')
	plt.savefig(outdirectory+'polang_map_planckwmap_Q.png')
	hp.mollview(umap,norm='hist')
	plt.savefig(outdirectory+'polang_map_planckwmap_U.png')
	plt.clf()
	plt.close()
	polang_map_planckwmap = calc_polang(qmap[:],umap[:])
	polang_map_planckwmap_unc = calc_polang_unc(qmap[:], umap[:], qmap_var[:], umap_var[:])

	# WMAP
	mapdata = hp.read_map('/Users/mpeel/Documents/maps/wmap9_tqu_v1.5_noise_v1.0_10k/64_60.0smoothed_wmap9beamNoise_22.8_512_2013_mKCMBunits.fits',field=None)
	qmap_wmap = hp.ud_grade(mapdata[1],nside_out,order_in='RING',order_out='RING')*commonmask
	umap_wmap = hp.ud_grade(mapdata[2],nside_out,order_in='RING',order_out='RING')*commonmask
	qmap_wmap_var = hp.ud_grade(mapdata[4],nside_out,order_in='RING',order_out='RING')*commonmask
	umap_wmap_var = hp.ud_grade(mapdata[6],nside_out,order_in='RING',order_out='RING')*commonmask
	hp.mollview(qmap_wmap,norm='hist')
	plt.savefig(outdirectory+'polang_map_wmapK_Q.png')
	hp.mollview(umap_wmap,norm='hist')
	plt.savefig(outdirectory+'polang_map_wmapK_U.png')
	plt.clf()
	plt.close()
	polang_map_wmap = calc_polang(qmap_wmap[:],umap_wmap[:])
	polang_map_wmap_unc = calc_polang_unc(qmap_wmap[:], umap_wmap[:], qmap_wmap_var[:], umap_wmap_var[:])

	threshold3 = 2.0 # mK
	if not staticmask:
		quickmask[qmap_wmap*(freqs[4]/wmap_freq)**index < -threshold3] = 0
		quickmask[qmap_wmap*(freqs[4]/wmap_freq)**index > threshold3] = 0
		quickmask[umap_wmap*(freqs[4]/wmap_freq)**index < -threshold3] = 0
		quickmask[umap_wmap*(freqs[4]/wmap_freq)**index > threshold3] = 0

	# Save the mask
	hp.write_map(outdirectory+'_quickmask.fits',quickmask,overwrite=True)
	hp.mollview(quickmask,norm='hist')
	plt.savefig(outdirectory+'_quickmask.png')
	plt.clf()
	plt.close()
	# exit()
	testmask = quickmask.copy()

	# mapdata = hp.read_map('/Users/mpeel/Documents/maps/wmap9_planck2018_tqu/512_60.0smoothed_PlanckR3fullbeam_28.4_1024_2018_mKCMBunits.fits',field=None)
	if planckmap == '2015':
		mapdata = hp.read_map('/Users/mpeel/Documents/maps/planck2015_tqu_v1.5_noise_v1.0/64_60.0smoothed_PlanckR2fullbeambpcorrNoise_28.4_256_2015_mKCMBunits.fits',field=None)
	elif planckmap == '2018':
		mapdata = hp.read_map('/Users/mpeel/Documents/maps/planck2018_tqu_v1.5_noise_v1.0_10k/64_60.0smoothed_PlanckR3fullbeamNoise_28.4_1024_2018_mKCMBunits.fits',field=None)
	else:
		# npipe
		mapdata = hp.read_map('/Users/mpeel/Documents/maps/planck2020_tqu_v1.5_noise_v1.0_10k/64_60.0smoothed_PlanckR4fullbeamnodpNoise_28.4_1024_2020_mKCMBunits.fits',field=None)
	# mapdata1 = hp.pixelfunc.reorder(mapdata[1],n2r=True)
	# mapdata2 = hp.pixelfunc.reorder(mapdata[2],n2r=True)
	qmap_planck = hp.ud_grade(mapdata[1],nside_out,order_in='RING',order_out='RING')*commonmask
	umap_planck = hp.ud_grade(mapdata[2],nside_out,order_in='RING',order_out='RING')*commonmask
	qmap_planck_var = hp.ud_grade(mapdata[4],nside_out,order_in='RING',order_out='RING')*commonmask
	umap_planck_var = hp.ud_grade(mapdata[6],nside_out,order_in='RING',order_out='RING')*commonmask
	hp.mollview(qmap_planck,norm='hist')
	plt.savefig(outdirectory+'polang_map_planck30_Q.png')
	hp.mollview(umap_planck,norm='hist')
	plt.savefig(outdirectory+'polang_map_planck30_U.png')
	plt.clf()
	plt.close()
	polang_map_planck30 = calc_polang(qmap_planck[:],umap_planck[:])
	polang_map_planck30_unc = calc_polang_unc(qmap_planck[:], umap_planck[:], qmap_planck_var[:], umap_planck_var[:])

	plot_tt(qmap_wmap[quickmask==1],qmap_planck[quickmask==1],outdirectory+'tt_wmap_planck_Q.png')
	plot_tt(umap_wmap[quickmask==1],umap_planck[quickmask==1],outdirectory+'tt_wmap_planck_U.png')

	# offsetmap = run_offset_map(qmap_wmap*(freqs[4]/wmap_freq)**index,qmap_311,np.sqrt(var_q_311),nsidemask=nsidemask,nside=8,outputtt=outdirectory+'ttplots/tt_311_wmap_Q_sigma_',outputmap=outdirectory+'_offsetQ311wmap.png')
	# print(np.average(offsetmap[offsetmap != 0.0]))
	# offsetmap = run_offset_map(umap_wmap*(freqs[4]/wmap_freq)**index,umap_311,np.sqrt(var_u_311),nsidemask=nsidemask,nside=8,outputtt=outdirectory+'ttplots/tt_311_wmap_U_sigma_',outputmap=outdirectory+'_offsetU311wmap.png')
	# print(np.average(offsetmap[offsetmap != 0.0]))


	offset_values_Q = np.zeros(nummaps)
	offset_values_U = np.zeros(nummaps)
	offset_values_Q_err = np.zeros(nummaps)
	offset_values_U_err = np.zeros(nummaps)
	offset_values_Q_planck = np.zeros(nummaps)
	offset_values_U_planck = np.zeros(nummaps)
	offset_values_Q_planck_err = np.zeros(nummaps)
	offset_values_U_planck_err = np.zeros(nummaps)

	# Get the offset between 11GHz and Planck/WMAP
	if use_variance:
		fit,fiterr=plot_tt(qmap_planck[quickmask==1]*(freqs[4]/planck_freq)**index,qmap_311[quickmask==1],outdirectory+'tt_311_planck_Q_sigma.png',sigma=np.sqrt(var_q_311[quickmask==1]),sigma_x=np.sqrt(qmap_planck_var[quickmask==1])*(freqs[4]/planck_freq)**index)
	else:
		fit,fiterr=plot_tt(qmap_planck[quickmask==1]*(freqs[4]/planck_freq)**index,qmap_311[quickmask==1],outdirectory+'tt_311_planck_Q_sigma.png')
	# print(fit,fiterr)
	offset_values_Q_planck[i] = fit[1]
	offset_values_Q_planck_err[i] = fiterr[1]
	if applyoffsets == 'planck' and fit[1] != 0.0:
		qmap_311 = qmap_311 - fit[1]

	if use_variance:
		fit,fiterr=plot_tt(umap_planck[quickmask==1]*(freqs[4]/planck_freq)**index,umap_311[quickmask==1],outdirectory+'tt_311_planck_U_sigma.png',sigma=np.sqrt(var_u_311[quickmask==1]),sigma_x=np.sqrt(umap_planck_var[quickmask==1])*(freqs[4]/planck_freq)**index)
	else:
		fit,fiterr=plot_tt(umap_planck[quickmask==1]*(freqs[4]/planck_freq)**index,umap_311[quickmask==1],outdirectory+'tt_311_planck_U_sigma.png')
	# print(fit,fiterr)
	offset_values_U_planck[i] = fit[1]
	offset_values_U_planck_err[i] = fiterr[1]
	if applyoffsets == 'planck' and fit[1] != 0.0:
		umap_311 = umap_311 - fit[1]

	if use_variance:
		fit,fiterr=plot_tt(qmap_wmap[quickmask==1]*(freqs[4]/wmap_freq)**index,qmap_311[quickmask==1],outdirectory+'tt_311_wmap_Q_sigma.png',sigma=np.sqrt(var_q_311[quickmask==1]),sigma_x=np.sqrt(qmap_wmap_var[quickmask==1])*(freqs[4]/wmap_freq)**index)
	else:
		fit,fiterr=plot_tt(qmap_wmap[quickmask==1]*(freqs[4]/wmap_freq)**index,qmap_311[quickmask==1],outdirectory+'tt_311_wmap_Q_sigma.png')
	# print(fit,fiterr)
	offset_values_Q[i] = fit[1]
	offset_values_Q_err[i] = fiterr[1]
	if applyoffsets == 'wmap' and fit[1] != 0.0:
		qmap_311 = qmap_311 - fit[1]
	if use_variance:
		fit,fiterr=plot_tt(umap_wmap[quickmask==1]*(freqs[4]/wmap_freq)**index,umap_311[quickmask==1],outdirectory+'tt_311_wmap_U_sigma.png',sigma=np.sqrt(var_u_311[quickmask==1]),sigma_x=np.sqrt(umap_wmap_var[quickmask==1])*(freqs[4]/wmap_freq)**index)
	else:
		fit,fiterr=plot_tt(umap_wmap[quickmask==1]*(freqs[4]/wmap_freq)**index,umap_311[quickmask==1],outdirectory+'tt_311_wmap_U_sigma.png')
	# print(fit,fiterr)
	offset_values_U[i] = fit[1]
	offset_values_U_err[i] = fiterr[1]
	if applyoffsets == 'wmap' and fit[1] != 0.0:
		umap_311 = umap_311 - fit[1]
	polang_map_11 = calc_polang(qmap_311[:],umap_311[:])
	if use_variance:
		polang_unc_map_11 = calc_polang_unc(qmap_311[:],umap_311[:],np.sqrt(var_q_311),np.sqrt(var_u_311))

	combhist = np.zeros((nummaps,89))
	diff_to_11 = np.zeros(nummaps)
	diff_to_planck = np.zeros(nummaps)
	diff_to_wmap = np.zeros(nummaps)
	diff_to_planckwmap = np.zeros(nummaps)
	diff_to_11_2 = np.zeros(nummaps)
	diff_to_planck_2 = np.zeros(nummaps)
	diff_to_wmap_2 = np.zeros(nummaps)
	diff_to_planckwmap_2 = np.zeros(nummaps)
	diff_to_11_std = np.zeros(nummaps)
	diff_to_planck_std = np.zeros(nummaps)
	diff_to_wmap_std = np.zeros(nummaps)
	diff_to_planckwmap_std = np.zeros(nummaps)
	diff_to_11_2_std = np.zeros(nummaps)
	diff_to_planck_2_std = np.zeros(nummaps)
	diff_to_wmap_2_std = np.zeros(nummaps)
	diff_to_planckwmap_2_std = np.zeros(nummaps)
	# Start running through, skip MFI horn 1
	for i in range(2,nummaps):
		# print(maps[i])
		skipthis = False
		try:
			mapdata = hp.read_map(indirectory+maps[i],field=None)
		except:
			print('No map found! Assuming hp.UNSEEN')
			skipthis = True
			mapdata = np.ones((7,len(qmap_311)))*hp.UNSEEN
		if skipthis:
			continue
		mapdata[mapdata < -1e10] = hp.UNSEEN
		qmap = hp.ud_grade(mapdata[1],nside_out,order_in='RING',order_out='RING')*commonmask
		umap = hp.ud_grade(mapdata[2],nside_out,order_in='RING',order_out='RING')*commonmask
		if i > 3 and not staticmask:
			quickmask[qmap < -1e10] = 0
			quickmask[umap < -1e10] = 0
			qmap[qmap < -1e10] = hp.UNSEEN
			umap[umap < -1e10] = hp.UNSEEN
			testmask[qmap == hp.UNSEEN] = 0
			testmask[umap == hp.UNSEEN] = 0

		hp.mollview(qmap,norm='hist')
		plt.savefig(outdirectory+'test_'+maps[i]+'_Q.png')
		hp.mollview(umap,norm='hist')
		plt.savefig(outdirectory+'test_'+maps[i]+'_U.png')
		plt.clf()
		plt.close()
		# Get the variance maps
		if use_variance:
			if 'nov2019' in datestr or 'mar2020' in datestr or 'apr2020' in datestr or 'nov2020' in datestr:
				if 'FG_QJT' in maps[i]:
					var = hp.read_map(indirectory+maps[i].replace('FG_QJT','FG_wei_QJT'),field=None)
					var_q = 1.0/var[1]
					var_u = 1.0/var[2]
					var_q = hp.ud_grade(var_q, nside_out, power=0)
					var_u = hp.ud_grade(var_u, nside_out, power=0)
				elif varmaps[i] != '':
					var = hp.read_map('/Volumes/Toshiba5TB2/maps/quijote_201911/noise/'+varmaps[i],field=None)
					var_q = var[1].copy()**2
					var_u = var[2].copy()**2
				else:
					var_q = np.ones(hp.pixelfunc.nside2npix(64))
					var_u = np.ones(hp.pixelfunc.nside2npix(64))
			else:
				if newformat:
					var_q = mapdata[4]
					var_u = mapdata[6]
				else:
					var_q = hp.read_map(indirectory+maps[i].replace('60.0s','60.00s').replace('_mKCMBunits','_weight_1_variance'),field=None)
					var_u = hp.read_map(indirectory+maps[i].replace('60.0s','60.00s').replace('_mKCMBunits','_weight_2_variance'),field=None)
				var_q = hp.ud_grade(var_q, nside_out, power=0)
				var_u = hp.ud_grade(var_u, nside_out, power=0)
			var_q[var_q < 0.0] = 10000.0
			var_u[var_u < 0.0] = 10000.0
			var_q[var_q==0.0] = 10000.0
			var_u[var_u==0.0] = 10000.0
			var_q[~np.isfinite(var_q)] = 10000.0
			var_u[~np.isfinite(var_u)] = 10000.0

			# offsetmap = run_offset_map(qmap_311,qmap,np.sqrt(var_q),sigmamap_x=np.sqrt(var_q_311),nsidemask=nsidemask,nside=8,outputtt=outdirectory+'ttplots/tt_311_'+maps[i]+'_Q_3_',outputmap=outdirectory+'_offset_Q311_'+maps[i]+'.png')

			# Check for offsets vs. the 11GHz maps
			# print(len(qmap_311))
			# print(len(qmap))
			# print(len(var_q))
			# print(len(var_q_311))
			fit,fiterr=plot_tt(qmap_311[quickmask==1],qmap[quickmask==1],outdirectory+'tt_311_'+maps[i]+'_Q_3.png',sigma=np.sqrt(var_q[quickmask==1]),sigma_x=np.sqrt(var_q_311[quickmask==1]))
		else:
			fit,fiterr=plot_tt(qmap_311[quickmask==1],qmap[quickmask==1],outdirectory+'tt_311_'+maps[i]+'_Q_3.png')
		# print(fit,fiterr)

		if (applyoffsets == 'wmap' or applyoffsets == 'planck') and fit[1] != 0.0:
			qmap = qmap-fit[1]
		if i != 4:
			offset_values_Q[i] = fit[1]
			offset_values_Q_err[i] = fiterr[1]

		if use_variance:
			fit,fiterr=plot_tt(umap_311[quickmask==1],umap[quickmask==1],outdirectory+'tt_311_'+maps[i]+'_U_3.png',sigma=np.sqrt(var_u[quickmask==1]),sigma_x=np.sqrt(var_u_311[quickmask==1]))

			# offsetmap = run_offset_map(umap_311,umap,np.sqrt(var_u),sigmamap_x=np.sqrt(var_u_311),nsidemask=nsidemask,nside=8,outputtt=outdirectory+'ttplots/tt_311_'+maps[i]+'_U_3_',outputmap=outdirectory+'_offset_U311_'+maps[i]+'.png')

		else:
			fit,fiterr=plot_tt(umap_311[quickmask==1],umap[quickmask==1],outdirectory+'tt_311_'+maps[i]+'_U_3.png')
		# print(fit,fiterr)
		if (applyoffsets == 'wmap' or applyoffsets == 'planck') and fit[1] != 0.0:
			umap = umap-fit[1]
		if i != 4:
			offset_values_U[i] = fit[1]
			offset_values_U_err[i] = fiterr[1]
		polang_map = calc_polang(qmap,umap)
		if use_variance:
			polang_unc_map = calc_polang_unc(qmap,umap,np.sqrt(var_q),np.sqrt(var_u))
		polang2 = polang_map.copy()
		polang2[quickmask==0] = 0.0
		polang2[polang2 > 180.0] = 0.0
		hp.mollview(polang2)
		plt.savefig(outdirectory+'test_polang_map_'+maps[i]+'.png')
		plt.clf()
		plt.close()



		# print(np.median(polang_map[quickmask==1]))
		hist = np.histogram(polang_map[quickmask == 1], bins=np.arange(0.0,90.0,1.0))
		combhist[i][:] = hist[0][:]
		plt.xlim(0,90)
		plt.plot(hist[1][:-1],hist[0])
		plt.savefig(outdirectory+'polang_'+maps[i]+'.pdf')
		plt.clf()
		plt.close()

		# hp.gnomview(qmap,reso=15.0)
		# plt.savefig(outdirectory+'galcen_Q_'+maps[i]+'.pdf')
		# hp.gnomview(umap,reso=15.0)
		# plt.savefig(outdirectory+'galcen_U_'+maps[i]+'.pdf')
		# hp.gnomview(polang_map,reso=15.0)
		# plt.savefig(outdirectory+'galcen_ang_'+maps[i]+'.pdf')

		polmap_temp = np.sqrt(qmap**2+umap**2)
		polmap_temp[quickmask==0] = hp.UNSEEN
		qmap[quickmask==0] = hp.UNSEEN
		umap[quickmask==0] = hp.UNSEEN
		# hp.mollview(qmap,norm='hist')
		# plt.savefig(outdirectory+'polang_map_'+maps[i]+'_Q.pdf')
		# hp.mollview(umap,norm='hist')
		# plt.savefig(outdirectory+'polang_map_'+maps[i]+'_U.pdf')
		# plt.clf()
		hp.mollview(polmap_temp,norm='hist')
		plt.savefig(outdirectory+'polang_map_'+maps[i]+'_P.pdf')

		if i == 2:
			polang_map_17 = polang_map.copy()
			polmap_17 = polmap_temp.copy()
			Q_17 = qmap.copy()
			U_17 = umap.copy()
			Q_17_unc = np.sqrt(var_q)
			U_17_unc = np.sqrt(var_u)
			if use_variance:
				polang_unc_map_17 = polang_unc_map.copy()
		elif i == 3:
			polang_map_19 = polang_map.copy()
			polmap_19 = polmap_temp.copy()
			Q_19 = qmap.copy()
			U_19 = umap.copy()
			Q_19_unc = np.sqrt(var_q)
			U_19_unc = np.sqrt(var_u)
			if use_variance:
				polang_unc_map_19 = polang_unc_map.copy()
		elif i == 5:
			polang_diff_temp = polang_map - polang_map_11
			polang_diff_temp[quickmask==0] = hp.UNSEEN
			hp.gnomview(polang_diff_temp,reso=15.0,min=-20,max=10,cmap=plt.get_cmap('jet'))
			plt.savefig(outdirectory+'galcen_diff_11_13.pdf')
			plt.clf()
			plt.close()

		elif i == 6:
			try:
				logfile.write("Difference at 17GHz:\n")
				logfile.write(str(np.median(dodiff(polang_map[quickmask==1],polang_map_17[quickmask==1])))+"\n")
				logfile.write('In amplitude median: ' + str(np.median(polmap_temp[quickmask==1]/polmap_17[quickmask==1]))+"\n")
				logfile.write('In amplitude sum: ' + str(np.sum(polmap_temp[quickmask==1])/np.sum(polmap_17[quickmask==1]))+"\n")
				logfile.write('In Q median: ' + str(np.median(qmap[quickmask==1]/Q_17[quickmask==1]))+"\n")
				logfile.write('In Q sum: ' + str(np.sum(qmap[quickmask==1])/np.sum(Q_17[quickmask==1]))+"\n")
				logfile.write('In U median: ' + str(np.median(umap[quickmask==1]/U_17[quickmask==1]))+"\n")
				logfile.write('In U sum: ' + str(np.sum(umap[quickmask==1])/np.sum(U_17[quickmask==1]))+"\n")

				logfile.write('In Q median: ' + str(np.sqrt(np.median(qmap[quickmask==1]**2/Q_17[quickmask==1]**2)))+"\n")
				logfile.write('In Q sum: ' + str(np.sqrt(np.sum(qmap[quickmask==1]**2)/np.sum(Q_17[quickmask==1]**2)))+"\n")
				logfile.write('In U median: ' + str(np.sqrt(np.median(umap[quickmask==1]**2/U_17[quickmask==1]**2)))+"\n")
				logfile.write('In U sum: ' + str(np.sqrt(np.sum(umap[quickmask==1]**2)/np.sum(U_17[quickmask==1]**2)))+"\n")

				fit,fiterr=plot_tt(Q_17[quickmask==1],qmap[quickmask==1],outdirectory+'tt_Qdiff_17.png',sigma=np.sqrt(var_q[quickmask==1]),sigma_x=Q_17_unc[quickmask==1])
				fit,fiterr=plot_tt(U_17[quickmask==1],umap[quickmask==1],outdirectory+'tt_Udiff_17.png',sigma=np.sqrt(var_u[quickmask==1]),sigma_x=U_17_unc[quickmask==1])
				plt.clf()
				plt.close()
				# if use_variance:
					# print(np.average(polang_map[quickmask==1]-polang_map_17[quickmask==1],weights=1.0/(polang_unc_map[quickmask==1]**2.0+polang_unc_map_17[quickmask==1]**2.0)))
				hist = np.histogram(dodiff(polang_map[quickmask==1],polang_map_17[quickmask==1]), bins=np.arange(0.0,90.0,1.0))
				plt.xlim(0,90)
				plt.title('Difference between 17GHz channels')
				plt.xlabel('Angle difference [deg]')
				plt.ylabel('Count')
				plt.plot(hist[1][:-1],hist[0])
				plt.savefig(outdirectory+'polang_diff17.png')
				plt.clf()
				plt.close()
			except:
				null = 0
		elif i == 7:
			try:
				logfile.write("Difference at 19GHz:\n")
				logfile.write(str(np.median(dodiff(polang_map[quickmask==1],polang_map_19[quickmask==1])))+"\n")
				logfile.write('In amplitude median: ' + str(np.median(polmap_temp[quickmask==1]/polmap_19[quickmask==1]))+"\n")
				logfile.write('In amplitude sum: ' + str(np.sum(polmap_temp[quickmask==1])/np.sum(polmap_19[quickmask==1]))+"\n")

				logfile.write('In Q median: ' + str(np.median(qmap[quickmask==1]/Q_19[quickmask==1]))+"\n")
				logfile.write('In Q sum: ' + str(np.sum(qmap[quickmask==1])/np.sum(Q_19[quickmask==1]))+"\n")
				logfile.write('In U median: ' + str(np.median(umap[quickmask==1]/U_19[quickmask==1]))+"\n")
				logfile.write('In U sum: ' + str(np.sum(umap[quickmask==1])/np.sum(U_19[quickmask==1]))+"\n")

				logfile.write('In Q median: ' + str(np.sqrt(np.median(qmap[quickmask==1]**2/Q_19[quickmask==1]**2)))+"\n")
				logfile.write('In Q sum: ' + str(np.sqrt(np.sum(qmap[quickmask==1]**2)/np.sum(Q_19[quickmask==1]**2)))+"\n")
				logfile.write('In U median: ' + str(np.sqrt(np.median(umap[quickmask==1]**2/U_19[quickmask==1]**2)))+"\n")
				logfile.write('In U sum: ' + str(np.sqrt(np.sum(umap[quickmask==1]**2)/np.sum(U_19[quickmask==1]**2)))+"\n")
				plt.clf()
				plt.close()
				fit,fiterr=plot_tt(Q_19[quickmask==1],qmap[quickmask==1],outdirectory+'tt_Qdiff_19.png',sigma=np.sqrt(var_q[quickmask==1]),sigma_x=Q_19_unc[quickmask==1])
				fit,fiterr=plot_tt(U_19[quickmask==1],umap[quickmask==1],outdirectory+'tt_Udiff_19.png',sigma=np.sqrt(var_u[quickmask==1]),sigma_x=U_19_unc[quickmask==1])


				# if use_variance:
					# print(np.average(polang_map[quickmask==1]-polang_map_19[quickmask==1],weights=1.0/(polang_unc_map[quickmask==1]**2.0+polang_unc_map_19[quickmask==1]**2.0)))
				hist = np.histogram(dodiff(polang_map[quickmask==1],polang_map_19[quickmask==1]), bins=np.arange(0.0,90.0,1.0))
				plt.xlim(0,90)
				plt.title('Difference between 19GHz channels')
				plt.xlabel('Angle difference [deg]')
				plt.ylabel('Count')
				plt.plot(hist[1][:-1],hist[0])
				plt.savefig(outdirectory+'polang_diff19.png')
				plt.clf()
				plt.close()
			except:
				null = 0

		diff_to_11[i] = np.median(dodiff(polang_map[quickmask==1],polang_map_11[quickmask==1]))
		diff_to_11_std[i] = calc_std_over_n(dodiff(polang_map[quickmask==1],polang_map_11[quickmask==1]))
		# print(diff_to_11[i])
		if use_variance:
			try:
				vals = np.average(dodiff(polang_map[quickmask==1],polang_map_11[quickmask==1]),weights=1.0/(polang_unc_map[quickmask==1]**2.0+polang_unc_map_11[quickmask==1]**2.0),returned=True)
				diff_to_11_2[i] = vals[0]
				diff_to_11_2_std[i] = np.sqrt(1.0/vals[1])
			except:
				vals = [0.0, 0.0]
		if i != 4:
			hist = np.histogram(dodiff(polang_map[quickmask==1],polang_map_11[quickmask==1]), bins=np.arange(0.0,90.0,1.0))
			plt.xlim(0,90)
			plt.title('Difference between 11GHz and ' + str(freqs[i]) + 'GHz (' + instrument[i] + ')')
			plt.xlabel('Angle difference [deg]')
			plt.ylabel('Count')
			plt.plot(hist[1][:-1],hist[0])
			plt.savefig(outdirectory+'polang_diff11_'+maps[i]+'.png')
			plt.clf()
			plt.close()

		# print('Compare to Planck+WMAP:')
		diff_to_planckwmap[i] = np.median(dodiff(polang_map[quickmask==1],polang_map_planckwmap[quickmask==1]))
		diff_to_planckwmap_std[i] = calc_std_over_n(dodiff(polang_map[quickmask==1],polang_map_planckwmap[quickmask==1]))
		# print(diff_to_planckwmap[i])
		if use_variance:
			try:
				vals = np.average(dodiff(polang_map[quickmask==1],polang_map_planckwmap[quickmask==1]),weights=1.0/(polang_unc_map[quickmask==1]**2.0+polang_map_planckwmap_unc[quickmask==1]**2),returned=True)
				diff_to_planckwmap_2[i] = vals[0]
				diff_to_planckwmap_2_std[i] = np.sqrt(1.0/vals[1])
			except:
				vals = [0.0, 0.0]
		hist = np.histogram(dodiff(polang_map[quickmask==1],polang_map_planckwmap[quickmask==1]), bins=np.arange(0.0,90.0,1.0))
		plt.xlim(0,90)
		plt.title('Difference between WMAP+Planck and ' + str(freqs[i]) + 'GHz (' + instrument[i] + ')')
		plt.xlabel('Angle difference [deg]')
		plt.ylabel('Count')
		plt.plot(hist[1][:-1],hist[0])
		plt.savefig(outdirectory+'polang_diffplanckwmap_'+maps[i]+'.png')
		plt.clf()
		plt.close()

		# print('Compare to WMAP:')
		diff_to_wmap[i] = np.median(dodiff(polang_map[quickmask==1],polang_map_wmap[quickmask==1]))
		diff_to_wmap_std[i] = calc_std_over_n(dodiff(polang_map[quickmask==1],polang_map_wmap[quickmask==1]))
		# print(diff_to_wmap[i])
		if use_variance:
			try:
				vals = np.average(dodiff(polang_map[quickmask==1],polang_map_wmap[quickmask==1]),weights=1.0/(polang_unc_map[quickmask==1]**2.0+polang_map_wmap_unc[quickmask==1]**2),returned=True)
				diff_to_wmap_2[i] = vals[0]
				diff_to_wmap_2_std[i] = np.sqrt(1.0/vals[1])
			except:
				vals = [0.0, 0.0]
		hist = np.histogram(dodiff(polang_map[quickmask==1],polang_map_wmap[quickmask==1]), bins=np.arange(0.0,90.0,1.0))
		plt.xlim(0,90)
		plt.title('Difference between WMAPK and ' + str(freqs[i]) + 'GHz (' + instrument[i] + ')')
		plt.xlabel('Angle difference [deg]')
		plt.ylabel('Count')
		plt.plot(hist[1][:-1],hist[0])
		plt.savefig(outdirectory+'polang_diffwmap_'+maps[i]+'.png')
		plt.clf()
		plt.close()

		# print('Compare to Planck:')
		diff_to_planck[i] = np.median(dodiff(polang_map[quickmask==1],polang_map_planck30[quickmask==1]))
		diff_to_planck_std[i] = calc_std_over_n(dodiff(polang_map[quickmask==1],polang_map_planck30[quickmask==1]))
		# print(diff_to_planck[i])
		if use_variance:
			try:
				vals = np.average(dodiff(polang_map[quickmask==1],polang_map_planck30[quickmask==1]),weights=1.0/(polang_unc_map[quickmask==1]**2.0+polang_map_planck30_unc[quickmask==1]**2.0),returned=True)
				diff_to_planck_2[i] = vals[0]
				diff_to_planck_2_std[i] = np.sqrt(1.0/vals[1])
			except:
				vals = [0.0, 0.0]
		hist = np.histogram(dodiff(polang_map[quickmask==1],polang_map_planck30[quickmask==1]), bins=np.arange(0.0,90.0,1.0))
		plt.xlim(0,90)
		plt.title('Difference between Planck30 and ' + str(freqs[i]) + 'GHz (' + instrument[i] + ')')
		plt.xlabel('Angle difference [deg]')
		plt.ylabel('Count')
		plt.plot(hist[1][:-1],hist[0])
		plt.savefig(outdirectory+'polang_diffplanck_'+maps[i]+'.png')
		plt.clf()
		plt.close()

	# print(combhist)
	logfile.write('Offsets in Q and U (x1000):\n')
	np.set_printoptions(formatter={'float': '{: 0.5f}'.format})
	logfile.write(str(offset_values_Q[2:]*1000)+"\n")
	logfile.write(str(offset_values_U[2:]*1000)+"\n")
	logfile.write('Uncertainties (x1000):\n')
	logfile.write(str(offset_values_Q_err[2:]*1000)+"\n")
	logfile.write(str(offset_values_U_err[2:]*1000)+"\n")
	logfile.write('Planck to 311:\n')
	logfile.write("{:.2f}".format(offset_values_Q_planck[4]*1000)+"\n")
	logfile.write("{:.2f}".format(offset_values_U_planck[4]*1000)+"\n")
	logfile.write("{:.2f}".format(offset_values_Q_planck_err[4]*1000)+"\n")
	logfile.write("{:.2f}".format(offset_values_U_planck_err[4]*1000)+"\n")
	np.set_printoptions(formatter={'float': '{: 0.3f}'.format})
	logfile.write('Diff to 11:'+"\n")
	logfile.write(str(diff_to_11[2:])+"\n")
	logfile.write(str(diff_to_11_std[2:])+"\n")
	if use_variance:
		logfile.write(str(diff_to_11_2[2:])+"\n")
		logfile.write(str(diff_to_11_2_std[2:])+"\n")

	logfile.write('Diff to Planck+WMAP:\n')
	logfile.write(str(diff_to_planckwmap[2:])+"\n")
	logfile.write(str(diff_to_planckwmap_std[2:])+"\n")
	if use_variance:
		logfile.write(str(diff_to_planckwmap_2[2:])+"\n")
		logfile.write(str(diff_to_planckwmap_2_std[2:])+"\n")

	logfile.write('Diff to Planck:\n')
	logfile.write(str(diff_to_planck[2:])+"\n")
	logfile.write(str(diff_to_planck_std[2:])+"\n")
	if use_variance:
		logfile.write(str(diff_to_planck_2[2:])+"\n")
		logfile.write(str(diff_to_planck_2_std[2:])+"\n")

	logfile.write('Diff to WMAP:\n')
	logfile.write(str(diff_to_wmap[2:])+"\n")
	logfile.write(str(diff_to_wmap_std[2:])+"\n")
	if use_variance:
		logfile.write(str(diff_to_wmap_2[2:])+"\n")
		logfile.write(str(diff_to_wmap_2_std[2:])+"\n")

	logfile.write('Difference between WMAPK and Planck30\n')
	logfile.write("{:.2f}".format(np.median(dodiff(polang_map_wmap[quickmask==1],polang_map_planck30[quickmask==1]))) + "+-" + "{:.2f}".format(calc_std_over_n(dodiff(polang_map_wmap[quickmask==1],polang_map_planck30[quickmask==1])))+"\n")
	# print(len(polang_map_wmap[quickmask==1]))
	if use_variance:
		try:
			vals = np.average(dodiff(polang_map_wmap[quickmask==1],polang_map_planck30[quickmask==1]),weights=1.0/(polang_map_wmap_unc[quickmask==1]**2.0+polang_map_planck30_unc[quickmask==1]**2.0),returned=True)
			logfile.write('Weighted:'+"\n")
			logfile.write("{:.2f}".format(vals[0]) + "+-" + "{:.2f}".format(np.sqrt(1.0/vals[1]))+"\n")
		except:
			null = 0

	hist = np.histogram(dodiff(polang_map_wmap[quickmask==1],polang_map_planck30[quickmask==1]), bins=np.arange(0.0,90.0,1.0))
	plt.xlim(0,90)
	plt.title('Difference between WMAPK and Planck30')
	plt.xlabel('Angle difference [deg]')
	plt.ylabel('Count')
	plt.plot(hist[1][:-1],hist[0])
	plt.savefig(outdirectory+'polang_diffplanck_wmap.png')
	plt.clf()
	plt.close()
	logfile.close()

	# Save the mask
	hp.write_map(outdirectory+'_testmask.fits',testmask,overwrite=True)
	hp.mollview(testmask,norm='hist')
	plt.savefig(outdirectory+'_testmask.png')
	plt.clf()
	plt.close()

	return [diff_to_wmap[2:], diff_to_wmap_std[2:], diff_to_wmap_2[2:], diff_to_wmap_2_std[2:]]

	# EOF
