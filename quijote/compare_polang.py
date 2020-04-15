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
from astrocode.fitspectrum.astroutils import *

def linfit(x, param):
	return param[0]*x+param[1]

def compute_residuals_linfit(param, x, y):
	model = linfit(x, param)
	residual = y - model
	return residual

def linfit2(x, A, B):
	return A*x+B

def linfit3(param, x):
	return param[0]*x+param[1]

# Debiasing, code from Federica
def debias_p_mas(Q, U, sigmaQ, sigmaU):
    p = np.sqrt(Q**2+U**2)
    b = np.sqrt((Q**2+sigmaU**2+U**2+sigmaQ**2))/p
    pmas = p-b**2*((1.-np.exp(-p**2/b**2))/2*p)
    sigmap = np.sqrt(0.5*(sigmaQ**2+sigmaU**2))
    SN = (pmas/sigmap > 3.8)
    print(SN)
    return(pmas, sigmap, SN)

# Debiasing, code from Federica
def debias_p_as(Q, U, sigmaQ, sigmaU):
    p = np.sqrt(Q**2+U**2)
    sigmap = np.sqrt((sigmaQ*Q)**2+(sigmaU*U)**2/(Q**2+U**2))
    pmas = 0.5*(p+np.sqrt(p**2-2.*sigmap**2))
    SN = (pmas/sigmap > 3.8)
    print(SN)
    return(pmas, sigmap, SN)

def calc_polang(Q, U):
	return 0.5*np.arctan2(Q, U) * 180 / np.pi

def calc_polang_unc(Q, U, Qerr, Uerr):
	unc_map = np.sqrt((Qerr**2)*(-0.5*U/(Q**2.0+U**2.0))**2.0 + (Uerr**2)*(-0.5*Q/(Q**2.0+U**2.0))**2.0) * 180 / np.pi
	unc_map[~np.isfinite(unc_map)] = 1000.0
	return unc_map

def plot_tt(vals1,vals2,outputname,sigma=np.empty(0),sigma_x=np.empty(0),leastsq=False,freq1=0,freq2=0,xlabel='',ylabel='',xyline=False):
	# print(sigma)
	# Do a fit
	params = [1.0,0]
	if leastsq:
		# least-squares fit without uncertainties
		param_est, cov_x, infodict, mesg_result, ret_value = optimize.leastsq(compute_residuals_linfit, params, args=(vals1, vals2),full_output=True)
		sigma_param_est = np.sqrt(np.diagonal(cov_x))
	elif sigma_x.size:
		# Do an odr fit
		odr_model = odr.Model(linfit3)
		dataset = odr.Data(vals1, vals2, wd=1.0/sigma_x**2, we=1.0/sigma**2)
		odr_run = odr.ODR(dataset, odr_model, beta0=params)
		out = odr_run.run()
		param_est = out.beta
		sigma_param_est = out.sd_beta
	elif sigma.size:
		# Do a curve fit to use the uncertainties
		param_est, cov_x = optimize.curve_fit(linfit2, vals1, vals2, params, sigma=sigma)
		sigma_param_est = np.sqrt(np.diagonal(cov_x))
	else:
		# Do a curve fit without uncertainties
		param_est, cov_x = optimize.curve_fit(linfit2, vals1, vals2, params)
		sigma_param_est = np.sqrt(np.diagonal(cov_x))
	
	#Plot the data and the results
	if sigma_x.size:
		plt.errorbar(vals1,vals2,yerr=sigma,xerr=sigma_x,fmt='.')
	elif sigma.size:
		plt.errorbar(vals1,vals2,yerr=sigma,fmt='.')
	else:
		plt.plot(vals1,vals2,'.')
	xvals=np.arange(np.min(vals1),np.max(vals1),(np.max(vals1)-np.min(vals1))/100.0)
	if np.isfinite(sigma_param_est[0]) and np.isfinite(sigma_param_est[1]):
		mesg_fit = (
		r'$A={:5.3e}\pm{:3.2e}$'.format(
			param_est[0], sigma_param_est[0]) + ','
		r'$B={:5.3e}\pm{:3.2e}$'.format(
			param_est[1], sigma_param_est[1]))# + ','
		if freq1 != freq2:
			beta = -np.log(param_est[0])/np.log(freq1/freq2)
			s_beta = sigma_param_est[0]/param_est[0] / np.abs(np.log(freq1/freq2))
			mesg_fit = mesg_fit + ', beta = ${:5.3f}\pm{:3.2f}$'.format(beta,s_beta)
			print(beta,s_beta)
		plt.plot(xvals,linfit(xvals,param_est),'g',label="Fit: " + mesg_fit)
	else:
		plt.plot(xvals,linfit(xvals,param_est),'g')
	if xyline:
		plt.plot(vals1,vals1,label='X=Y')
	plt.legend(prop={'size':8})
	if xlabel != '':
		plt.xlabel(xlabel)
	if ylabel != '':
		plt.ylabel(ylabel)
	plt.savefig(outputname)
	plt.clf()
	return [param_est, sigma_param_est]

def calc_std_over_n(vals):
	return np.std(vals)/np.sqrt(len(vals))

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
		hp.mollview(uncmap)
		plt.savefig(outputmap.replace('.','_unc.'))
		plt.clf()

	return offsetmap


# # First, do a quick test
# X,Xerr,Y,Yerr = np.loadtxt('polarised_intensity_23ghz_30ghz_fermi_bubble.dat',unpack=True)
# print(X)
# print(Xerr)
# print(Y)
# print(Yerr)
# print(plot_tt(X,Y,'test.png'))
# print(plot_tt(X,Y,'test_1D.png',sigma=Yerr))
# res = plot_tt(X,Y,'test_2D.png',sigma=Yerr,sigma_x=Xerr)
# print(res)
# print(res[0][0])
# print(np.log(res[0][0])/np.log(28.4/22.8))
# print(res[1][0]/res[0][0]/np.log(28.4/22.8))
# exit()

def compare_polang(prefix='mfi', date='201905',datestr='may2019',use_variance=True):
	nside = 512
	npix = hp.nside2npix(nside)
	nside_out=64
	npix_out=hp.nside2npix(nside_out)
	# maps = [prefix+'_'+datestr+'_mapsmth_11.0_1.fits',prefix+'_'+datestr+'_mapsmth_13.0_1.fits',prefix+'_'+datestr+'_mapsmth_17.0_2.fits',prefix+'_'+datestr+'_mapsmth_19.0_2.fits',prefix+'_'+datestr+'_mapsmth_11.0_3.fits',prefix+'_'+datestr+'_mapsmth_13.0_3.fits',prefix+'_'+datestr+'_mapsmth_17.0_4.fits',prefix+'_'+datestr+'_mapsmth_19.0_4.fits']
	varmaps = ['', '', 'std_H2_17_sm1deg_nside64.fits', 'std_H2_19_sm1deg_nside64.fits', 'std_H3_11_sm1deg_nside64.fits', 'std_H3_13_sm1deg_nside64.fits', 'std_H4_17_sm1deg_nside64.fits', 'std_H4_19_sm1deg_nside64.fits']
	maps = [str(nside)+'_60.00smoothed_'+prefix+'1_11.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'1_13.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'2_17.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'2_19.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'3_11.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.00smoothed_'+prefix+'3_13.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'4_17.0_512_'+date+'_mKCMBunits.fits',str(nside)+'_60.0smoothed_'+prefix+'4_19.0_512_'+date+'_mKCMBunits.fits']

	indirectory = '/Users/mpeel/Documents/maps/quijote_'+date+'/reform/'
	outdirectory = '/Users/mpeel/Documents/maps/quijote_'+date+'/analyse/'


	# maps = ['FG_QJT_nominal_pluz_haze_h1_11.0_512_1.0deg.fits','FG_QJT_nominal_pluz_haze_h1_13.0_512_1.0deg.fits','FG_QJT_nominal_pluz_haze_h2_17.0_512_1.0deg.fits','FG_QJT_nominal_pluz_haze_h2_19.0_512_1.0deg.fits','FG_QJT_nominal_pluz_haze_h3_11.0_512_1.0deg.fits','FG_QJT_nominal_pluz_haze_h3_13.0_512_1.0deg.fits','FG_QJT_nominal_pluz_haze_h4_17.0_512_1.0deg.fits','FG_QJT_nominal_pluz_haze_h4_19.0_512_1.0deg.fits']
	# indirectory = '/Users/mpeel/Documents/maps/quijote_'+date+'/haze/'
	# outdirectory = '/Users/mpeel/Documents/maps/quijote_'+date+'/analyse_haze/'

	nsidemask = np.asarray(nside_mask(nside_out,8))
	hp.mollview(nsidemask,cmap=plt.get_cmap('jet'))
	plt.savefig(outdirectory+'_nside_mask.png')

	nummaps = len(maps)
	freqs = [11,13,17,19,11,13,17,19]
	instrument = ['MFI1','MFI1','MFI2','MFI2','MFI3','MFI3','MFI4','MFI4']
	normfreq = 28.4
	index = 3.0
	commonmask = hp.read_map('/Users/mpeel/Documents/maps/quijote_masks/mask_quijote_ncp_lowdec_nside512.fits',field=None)
	commonmask = hp.read_map('/Users/mpeel/Documents/maps/quijote_masks/mask_quijote_ncp_satband_nside512.fits',field=None)
	# commonmask = hp.read_map(outdirectory+'mfi_commonmask.fits',field=None)
	commonmask = hp.ud_grade(commonmask,nside_out,order_in='RING',order_out='RING')

	# Make a quick mask to get rid of low signal-to-noise pixels
	use_threshold = False
	threshold = 0.8
	threshold2 = 2.5
	use_sn = False
	sn_ratio = 5.0
	use_polang_err = True
	polang_err_threshold = 5.0

	applyoffsets = True

	i = 4 # Use the 11GHz map for this
	mapdata = hp.read_map(indirectory+maps[i],field=None)
	qmap_311 = hp.ud_grade(mapdata[1],nside_out,order_in='RING',order_out='RING')*commonmask
	umap_311 = hp.ud_grade(mapdata[2],nside_out,order_in='RING',order_out='RING')*commonmask
	if use_variance:
		if 'nov2019' in datestr or 'mar2020' in datestr or 'apr2020' in datestr:
			if 'FG_QJT' in maps[i]:
				var = hp.read_map(indirectory+maps[i].replace('FG_QJT','FG_wei_QJT'),field=None)
				print(np.shape(var))
				var_q_311 = 1.0/var[1]
				var_u_311 = 1.0/var[2]
				var_q_311 = hp.ud_grade(var_q_311, nside_out, power=2)
				var_u_311 = hp.ud_grade(var_u_311, nside_out, power=2)
			elif varmaps[i] != '':
				var = hp.read_map(indirectory+'../../quijote_201911/noise/'+varmaps[i],field=None)
				var_q_311 = var[1].copy()**2
				var_u_311 = var[2].copy()**2
			else:
				var_q_311 = np.ones(hp.pixelfunc.nside2npix(64))
				var_u_311 = np.ones(hp.pixelfunc.nside2npix(64))
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
	if use_variance:
		polang_unc_map_11 = calc_polang_unc(qmap_311,umap_311,np.sqrt(var_q_311),np.sqrt(var_u_311))
		hp.mollview(polang_unc_map_11,cmap=plt.get_cmap('jet'),max=polang_err_threshold)
		plt.savefig(outdirectory+'polang_map_mfi11_unc.png')
		tempmap = polang_map_11.copy()
		tempmap[polang_unc_map_11 > polang_err_threshold] = 0.0
		hp.mollview(tempmap,cmap=plt.get_cmap('jet'))
		plt.savefig(outdirectory+'polang_map_mfi11_threshold.png')
		plt.clf()

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



	mapdata = hp.read_map('/Users/mpeel/Documents/maps/wmap_planck_pol/weighted_both_debwk_feb2015_tqu.fits',field=None)
	mapdata1 = hp.pixelfunc.reorder(-mapdata[1],n2r=True)
	mapdata2 = hp.pixelfunc.reorder(-mapdata[2],n2r=True)
	qmap = hp.ud_grade(mapdata1,nside_out,order_in='RING',order_out='RING')*commonmask
	umap = hp.ud_grade(mapdata2,nside_out,order_in='RING',order_out='RING')*commonmask
	hp.mollview(qmap,norm='hist')
	plt.savefig(outdirectory+'polang_map_planckwmap_Q.png')
	hp.mollview(umap,norm='hist')
	plt.savefig(outdirectory+'polang_map_planckwmap_U.png')
	plt.clf()
	polang_map_planckwmap = calc_polang(qmap[:],umap[:])

	mapdata = hp.read_map('/Users/mpeel/Documents/maps/wmap9/wmap_band_smth_iqumap_r9_9yr_K_v5.fits',field=None)
	# mapdata1 = hp.pixelfunc.reorder(mapdata[1],n2r=True)
	# mapdata2 = hp.pixelfunc.reorder(mapdata[2],n2r=True)
	qmap_wmap = hp.ud_grade(mapdata[1],nside_out,order_in='RING',order_out='RING')*commonmask
	umap_wmap = hp.ud_grade(mapdata[2],nside_out,order_in='RING',order_out='RING')*commonmask
	hp.mollview(qmap_wmap,norm='hist')
	plt.savefig(outdirectory+'polang_map_wmapK_Q.png')
	hp.mollview(umap_wmap,norm='hist')
	plt.savefig(outdirectory+'polang_map_wmapK_U.png')
	plt.clf()
	polang_map_wmap = calc_polang(qmap_wmap[:],umap_wmap[:])

	threshold3 = 2.5
	quickmask[qmap_wmap*(11.0/23.0)**-3.0 < -threshold3] = 0
	quickmask[qmap_wmap*(11.0/23.0)**-3.0 > threshold3] = 0
	quickmask[umap_wmap*(11.0/23.0)**-3.0 < -threshold3] = 0
	quickmask[umap_wmap*(11.0/23.0)**-3.0 > threshold3] = 0

	hp.write_map(outdirectory+'_quickmask.fits',quickmask,overwrite=True)
	hp.mollview(quickmask,norm='hist')
	plt.savefig(outdirectory+'_quickmask.png')
	plt.clf()
	# exit()


	offset_values_Q = np.zeros(nummaps)
	offset_values_U = np.zeros(nummaps)
	offset_values_Q_err = np.zeros(nummaps)
	offset_values_U_err = np.zeros(nummaps)


	mapdata = hp.read_map('/Users/mpeel/Documents/maps/wmap9_planck2018_tqu/512_60.0smoothed_PlanckR3fullbeam_28.4_1024_2018_mKCMBunits.fits',field=None)
	# mapdata1 = hp.pixelfunc.reorder(mapdata[1],n2r=True)
	# mapdata2 = hp.pixelfunc.reorder(mapdata[2],n2r=True)
	qmap_planck = hp.ud_grade(mapdata[1],nside_out,order_in='RING',order_out='RING')*commonmask
	umap_planck = hp.ud_grade(mapdata[2],nside_out,order_in='RING',order_out='RING')*commonmask
	hp.mollview(qmap_planck,norm='hist')
	plt.savefig(outdirectory+'polang_map_planck30_Q.png')
	hp.mollview(umap_planck,norm='hist')
	plt.savefig(outdirectory+'polang_map_planck30_U.png')
	plt.clf()
	polang_map_planck30 = calc_polang(qmap_planck[:],umap_planck[:])

	plot_tt(qmap_wmap[quickmask==1],qmap_planck[quickmask==1],outdirectory+'tt_wmap_planck_Q.png')
	plot_tt(umap_wmap[quickmask==1],umap_planck[quickmask==1],outdirectory+'tt_wmap_planck_U.png')

	# offsetmap = run_offset_map(qmap_wmap*(11.0/23.0)**-3.0,qmap_311,np.sqrt(var_q_311),nsidemask=nsidemask,nside=8,outputtt=outdirectory+'ttplots/tt_311_wmap_Q_sigma_',outputmap=outdirectory+'_offsetQ311wmap.png')
	# print(np.average(offsetmap[offsetmap != 0.0]))
	# offsetmap = run_offset_map(umap_wmap*(11.0/23.0)**-3.0,umap_311,np.sqrt(var_u_311),nsidemask=nsidemask,nside=8,outputtt=outdirectory+'ttplots/tt_311_wmap_U_sigma_',outputmap=outdirectory+'_offsetU311wmap.png')
	# print(np.average(offsetmap[offsetmap != 0.0]))



	# Get the offset between 11GHz and Planck/WMAP
	if use_variance:
		fit,fiterr=plot_tt(qmap_planck[quickmask==1]*(11.0/28.4)**-3.0,qmap_311[quickmask==1],outdirectory+'tt_311_planck_Q_sigma.png',sigma=np.sqrt(var_q_311[quickmask==1]))
	else:
		fit,fiterr=plot_tt(qmap_planck[quickmask==1]*(11.0/28.4)**-3.0,qmap_311[quickmask==1],outdirectory+'tt_311_planck_Q_sigma.png')
	print(fit,fiterr)
	if use_variance:
		fit,fiterr=plot_tt(umap_planck[quickmask==1]*(11.0/28.4)**-3.0,umap_311[quickmask==1],outdirectory+'tt_311_planck_U_sigma.png',sigma=np.sqrt(var_u_311[quickmask==1]))
	else:
		fit,fiterr=plot_tt(umap_planck[quickmask==1]*(11.0/28.4)**-3.0,umap_311[quickmask==1],outdirectory+'tt_311_planck_U_sigma.png')
	print(fit,fiterr)

	if use_variance:
		fit,fiterr=plot_tt(qmap_wmap[quickmask==1]*(11.0/23.0)**-3.0,qmap_311[quickmask==1],outdirectory+'tt_311_wmap_Q_sigma.png',sigma=np.sqrt(var_q_311[quickmask==1]))
	else:
		fit,fiterr=plot_tt(qmap_wmap[quickmask==1]*(11.0/23.0)**-3.0,qmap_311[quickmask==1],outdirectory+'tt_311_wmap_Q_sigma.png')
	print(fit,fiterr)
	offset_values_Q[i] = fit[1]
	offset_values_Q_err[i] = fiterr[1]
	if applyoffsets:
		qmap_311 = qmap_311 - fit[1]
	if use_variance:
		fit,fiterr=plot_tt(umap_wmap[quickmask==1]*(11.0/23.0)**-3.0,umap_311[quickmask==1],outdirectory+'tt_311_wmap_U_sigma.png',sigma=np.sqrt(var_u_311[quickmask==1]))
	else:
		fit,fiterr=plot_tt(umap_wmap[quickmask==1]*(11.0/23.0)**-3.0,umap_311[quickmask==1],outdirectory+'tt_311_wmap_U_sigma.png')
	print(fit,fiterr)
	offset_values_U[i] = fit[1]
	offset_values_U_err[i] = fiterr[1]
	if applyoffsets:
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
	for i in range(2,nummaps):
		print(maps[i])
		mapdata = hp.read_map(indirectory+maps[i],field=None)
		mapdata[mapdata < -1e10] = hp.UNSEEN
		qmap = hp.ud_grade(mapdata[1],nside_out,order_in='RING',order_out='RING')*commonmask
		umap = hp.ud_grade(mapdata[2],nside_out,order_in='RING',order_out='RING')*commonmask
		quickmask[qmap < -1e10] = 0
		quickmask[umap < -1e10] = 0
		qmap[qmap < -1e10] = hp.UNSEEN
		umap[umap < -1e10] = hp.UNSEEN
		hp.mollview(qmap,norm='hist')
		plt.savefig(outdirectory+'test_'+maps[i]+'_Q.png')
		hp.mollview(umap,norm='hist')
		plt.savefig(outdirectory+'test_'+maps[i]+'_U.png')
		# Get the variance maps
		if use_variance:
			if 'nov2019' in datestr or 'mar2020' in datestr or 'apr2020' in datestr:
				if 'FG_QJT' in maps[i]:
					var = hp.read_map(indirectory+maps[i].replace('FG_QJT','FG_wei_QJT'),field=None)
					var_q = 1.0/var[1]
					var_u = 1.0/var[2]
					var_q = hp.ud_grade(var_q, nside_out, power=0)
					var_u = hp.ud_grade(var_u, nside_out, power=0)
				elif varmaps[i] != '':
					var = hp.read_map(indirectory+'../../quijote_201911/noise/'+varmaps[i],field=None)
					var_q = var[1].copy()**2
					var_u = var[2].copy()**2
				else:
					var_q = np.ones(hp.pixelfunc.nside2npix(64))
					var_u = np.ones(hp.pixelfunc.nside2npix(64))
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
			fit,fiterr=plot_tt(qmap_311[quickmask==1],qmap[quickmask==1],outdirectory+'tt_311_'+maps[i]+'_Q_3.png',sigma=np.sqrt(var_q[quickmask==1]),sigma_x=np.sqrt(var_q_311[quickmask==1]))
		else:
			fit,fiterr=plot_tt(qmap_311[quickmask==1],qmap[quickmask==1],outdirectory+'tt_311_'+maps[i]+'_Q_3.png')
		print(fit,fiterr)

		if applyoffsets:
			qmap = qmap-fit[1]
		if i != 4:
			offset_values_Q[i] = fit[1]
			offset_values_Q_err[i] = fiterr[1]

		if use_variance:
			fit,fiterr=plot_tt(umap_311[quickmask==1],umap[quickmask==1],outdirectory+'tt_311_'+maps[i]+'_U_3.png',sigma=np.sqrt(var_u[quickmask==1]),sigma_x=np.sqrt(var_u_311[quickmask==1]))

			# offsetmap = run_offset_map(umap_311,umap,np.sqrt(var_u),sigmamap_x=np.sqrt(var_u_311),nsidemask=nsidemask,nside=8,outputtt=outdirectory+'ttplots/tt_311_'+maps[i]+'_U_3_',outputmap=outdirectory+'_offset_U311_'+maps[i]+'.png')

		else:
			fit,fiterr=plot_tt(umap_311[quickmask==1],umap[quickmask==1],outdirectory+'tt_311_'+maps[i]+'_U_3.png')
		print(fit,fiterr)
		if applyoffsets:
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



		print(np.median(polang_map[quickmask==1]))
		hist = np.histogram(polang_map[quickmask == 1], bins=np.arange(0.0,90.0,1.0))
		combhist[i][:] = hist[0][:]
		plt.xlim(0,90)
		plt.plot(hist[1][:-1],hist[0])
		plt.savefig(outdirectory+'polang_'+maps[i]+'.pdf')
		plt.clf()

		hp.gnomview(qmap,reso=15.0)
		plt.savefig(outdirectory+'galcen_Q_'+maps[i]+'.pdf')
		hp.gnomview(umap,reso=15.0)
		plt.savefig(outdirectory+'galcen_U_'+maps[i]+'.pdf')
		hp.gnomview(polang_map,reso=15.0)
		plt.savefig(outdirectory+'galcen_ang_'+maps[i]+'.pdf')

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

		elif i == 6:
			print("Difference at 17GHz:")
			print(np.median(polang_map[quickmask==1]-polang_map_17[quickmask==1]))
			print('In amplitude median: ' + str(np.median(polmap_temp[quickmask==1]/polmap_17[quickmask==1])))
			print('In amplitude sum: ' + str(np.sum(polmap_temp[quickmask==1])/np.sum(polmap_17[quickmask==1])))
			print('In Q median: ' + str(np.median(qmap[quickmask==1]/Q_17[quickmask==1])))
			print('In Q sum: ' + str(np.sum(qmap[quickmask==1])/np.sum(Q_17[quickmask==1])))
			print('In U median: ' + str(np.median(umap[quickmask==1]/U_17[quickmask==1])))
			print('In U sum: ' + str(np.sum(umap[quickmask==1])/np.sum(U_17[quickmask==1])))

			print('In Q median: ' + str(np.sqrt(np.median(qmap[quickmask==1]**2/Q_17[quickmask==1]**2))))
			print('In Q sum: ' + str(np.sqrt(np.sum(qmap[quickmask==1]**2)/np.sum(Q_17[quickmask==1]**2))))
			print('In U median: ' + str(np.sqrt(np.median(umap[quickmask==1]**2/U_17[quickmask==1]**2))))
			print('In U sum: ' + str(np.sqrt(np.sum(umap[quickmask==1]**2)/np.sum(U_17[quickmask==1]**2))))

			plt.close()
			plt.clf()
			fit,fiterr=plot_tt(Q_17[quickmask==1],qmap[quickmask==1],outdirectory+'tt_Qdiff_17.png',sigma=np.sqrt(var_q[quickmask==1]),sigma_x=Q_17_unc[quickmask==1])
			fit,fiterr=plot_tt(U_17[quickmask==1],umap[quickmask==1],outdirectory+'tt_Udiff_17.png',sigma=np.sqrt(var_u[quickmask==1]),sigma_x=U_17_unc[quickmask==1])

			if use_variance:
				print(np.average(polang_map[quickmask==1]-polang_map_17[quickmask==1],weights=1.0/(polang_unc_map[quickmask==1]**2.0+polang_unc_map_17[quickmask==1]**2.0)))
			hist = np.histogram(polang_map[quickmask==1]-polang_map_17[quickmask==1], bins=np.arange(0.0,90.0,1.0))
			plt.xlim(0,90)
			plt.title('Difference between 17GHz channels')
			plt.xlabel('Angle difference [deg]')
			plt.ylabel('Count')
			plt.plot(hist[1][:-1],hist[0])
			plt.savefig(outdirectory+'polang_diff17.png')
			plt.clf()
		elif i == 7:
			print("Difference at 19GHz:")
			print(np.median(polang_map[quickmask==1]-polang_map_19[quickmask==1]))
			print('In amplitude median: ' + str(np.median(polmap_temp[quickmask==1]/polmap_19[quickmask==1])))
			print('In amplitude sum: ' + str(np.sum(polmap_temp[quickmask==1])/np.sum(polmap_19[quickmask==1])))

			print('In Q median: ' + str(np.median(qmap[quickmask==1]/Q_19[quickmask==1])))
			print('In Q sum: ' + str(np.sum(qmap[quickmask==1])/np.sum(Q_19[quickmask==1])))
			print('In U median: ' + str(np.median(umap[quickmask==1]/U_19[quickmask==1])))
			print('In U sum: ' + str(np.sum(umap[quickmask==1])/np.sum(U_19[quickmask==1])))

			print('In Q median: ' + str(np.sqrt(np.median(qmap[quickmask==1]**2/Q_19[quickmask==1]**2))))
			print('In Q sum: ' + str(np.sqrt(np.sum(qmap[quickmask==1]**2)/np.sum(Q_19[quickmask==1]**2))))
			print('In U median: ' + str(np.sqrt(np.median(umap[quickmask==1]**2/U_19[quickmask==1]**2))))
			print('In U sum: ' + str(np.sqrt(np.sum(umap[quickmask==1]**2)/np.sum(U_19[quickmask==1]**2))))
			plt.close()
			plt.clf()
			fit,fiterr=plot_tt(Q_19[quickmask==1],qmap[quickmask==1],outdirectory+'tt_Qdiff_19.png',sigma=np.sqrt(var_q[quickmask==1]),sigma_x=Q_19_unc[quickmask==1])
			fit,fiterr=plot_tt(U_19[quickmask==1],umap[quickmask==1],outdirectory+'tt_Udiff_19.png',sigma=np.sqrt(var_u[quickmask==1]),sigma_x=U_19_unc[quickmask==1])


			if use_variance:
				print(np.average(polang_map[quickmask==1]-polang_map_19[quickmask==1],weights=1.0/(polang_unc_map[quickmask==1]**2.0+polang_unc_map_19[quickmask==1]**2.0)))
			hist = np.histogram(polang_map[quickmask==1]-polang_map_19[quickmask==1], bins=np.arange(0.0,90.0,1.0))
			plt.xlim(0,90)
			plt.title('Difference between 19GHz channels')
			plt.xlabel('Angle difference [deg]')
			plt.ylabel('Count')
			plt.plot(hist[1][:-1],hist[0])
			plt.savefig(outdirectory+'polang_diff19.png')
			plt.clf()

		diff_to_11[i] = np.median(polang_map[quickmask==1]-polang_map_11[quickmask==1])
		diff_to_11_std[i] = calc_std_over_n(polang_map[quickmask==1]-polang_map_11[quickmask==1])
		# print(diff_to_11[i])
		if use_variance:
			vals = np.average(polang_map[quickmask==1]-polang_map_11[quickmask==1],weights=1.0/(polang_unc_map[quickmask==1]**2.0+polang_unc_map_11[quickmask==1]**2.0),returned=True)
			diff_to_11_2[i] = vals[0]
			diff_to_11_2_std[i] = np.sqrt(1.0/vals[1])
		if i != 4:
			hist = np.histogram(polang_map[quickmask==1]-polang_map_11[quickmask==1], bins=np.arange(0.0,90.0,1.0))
			plt.xlim(0,90)
			plt.title('Difference between 11GHz and ' + str(freqs[i]) + 'GHz (' + instrument[i] + ')')
			plt.xlabel('Angle difference [deg]')
			plt.ylabel('Count')
			plt.plot(hist[1][:-1],hist[0])
			plt.savefig(outdirectory+'polang_diff11_'+maps[i]+'.png')
			plt.clf()

		# print('Compare to Planck+WMAP:')
		diff_to_planckwmap[i] = np.median(polang_map[quickmask==1]-polang_map_planckwmap[quickmask==1])
		diff_to_planckwmap_std[i] = calc_std_over_n(polang_map[quickmask==1]-polang_map_planckwmap[quickmask==1])
		# print(diff_to_planckwmap[i])
		if use_variance:
			vals = np.average(polang_map[quickmask==1]-polang_map_planckwmap[quickmask==1],weights=1.0/(polang_unc_map[quickmask==1]**2.0),returned=True)
			diff_to_planckwmap_2[i] = vals[0]
			diff_to_planckwmap_2_std[i] = np.sqrt(1.0/vals[1])
		hist = np.histogram(polang_map[quickmask==1]-polang_map_planckwmap[quickmask==1], bins=np.arange(0.0,90.0,1.0))
		plt.xlim(0,90)
		plt.title('Difference between WMAP+Planck and ' + str(freqs[i]) + 'GHz (' + instrument[i] + ')')
		plt.xlabel('Angle difference [deg]')
		plt.ylabel('Count')
		plt.plot(hist[1][:-1],hist[0])
		plt.savefig(outdirectory+'polang_diffplanckwmap_'+maps[i]+'.png')
		plt.clf()

		# print('Compare to WMAP:')
		diff_to_wmap[i] = np.median(polang_map[quickmask==1]-polang_map_wmap[quickmask==1])
		diff_to_wmap_std[i] = calc_std_over_n(polang_map[quickmask==1]-polang_map_wmap[quickmask==1])
		# print(diff_to_wmap[i])
		if use_variance:
			vals = np.average(polang_map[quickmask==1]-polang_map_wmap[quickmask==1],weights=1.0/(polang_unc_map[quickmask==1]**2.0),returned=True)
			diff_to_wmap_2[i] = vals[0]
			diff_to_wmap_2_std[i] = np.sqrt(1.0/vals[1])
		hist = np.histogram(polang_map[quickmask==1]-polang_map_wmap[quickmask==1], bins=np.arange(0.0,90.0,1.0))
		plt.xlim(0,90)
		plt.title('Difference between WMAPK and ' + str(freqs[i]) + 'GHz (' + instrument[i] + ')')
		plt.xlabel('Angle difference [deg]')
		plt.ylabel('Count')
		plt.plot(hist[1][:-1],hist[0])
		plt.savefig(outdirectory+'polang_diffwmap_'+maps[i]+'.png')
		plt.clf()

		# print('Compare to Planck:')
		diff_to_planck[i] = np.median(polang_map[quickmask==1]-polang_map_planck30[quickmask==1])
		diff_to_planck_std[i] = calc_std_over_n(polang_map[quickmask==1]-polang_map_planck30[quickmask==1])
		# print(diff_to_planck[i])
		if use_variance:
			vals = np.average(polang_map[quickmask==1]-polang_map_planck30[quickmask==1],weights=1.0/(polang_unc_map[quickmask==1]**2.0),returned=True)
			diff_to_planck_2[i] = vals[0]
			diff_to_planck_2_std[i] = np.sqrt(1.0/vals[1])
		hist = np.histogram(polang_map[quickmask==1]-polang_map_planck30[quickmask==1], bins=np.arange(0.0,90.0,1.0))
		plt.xlim(0,90)
		plt.title('Difference between Planck30 and ' + str(freqs[i]) + 'GHz (' + instrument[i] + ')')
		plt.xlabel('Angle difference [deg]')
		plt.ylabel('Count')
		plt.plot(hist[1][:-1],hist[0])
		plt.savefig(outdirectory+'polang_diffplanck_'+maps[i]+'.png')
		plt.clf()

	# print(combhist)
	print('Offsets in Q and U:')
	np.set_printoptions(formatter={'float': '{: 0.5f}'.format})
	print(offset_values_Q[2:]*1000)
	print(offset_values_U[2:]*1000)
	print('Uncertainties:')
	print(offset_values_Q_err[2:]*1000)
	print(offset_values_U_err[2:]*1000)
	np.set_printoptions(formatter={'float': '{: 0.3f}'.format})
	print('Diff to 11:')
	print(diff_to_11[2:])
	print(diff_to_11_std[2:])
	if use_variance:
		print(diff_to_11_2[2:])
		print(diff_to_11_2_std[2:])

	print('Diff to Planck+WMAP:')
	print(diff_to_planckwmap[2:])
	print(diff_to_planckwmap_std[2:])
	if use_variance:
		print(diff_to_planckwmap_2[2:])
		print(diff_to_planckwmap_2_std[2:])

	print('Diff to Planck:')
	print(diff_to_planck[2:])
	print(diff_to_planck_std[2:])
	if use_variance:
		print(diff_to_planck_2[2:])
		print(diff_to_planck_2_std[2:])

	print('Diff to WMAP:')
	print(diff_to_wmap[2:])
	print(diff_to_wmap_std[2:])
	if use_variance:
		print(diff_to_wmap_2[2:])
		print(diff_to_wmap_2_std[2:])

	print('Difference between WMAPK and Planck30')
	print(np.median(polang_map_wmap[quickmask==1]-polang_map_planck30[quickmask==1]))
	print(calc_std_over_n(polang_map_wmap[quickmask==1]-polang_map_planck30[quickmask==1]))
	print(len(polang_map_wmap[quickmask==1]))
	hist = np.histogram(polang_map_wmap[quickmask==1]-polang_map_planck30[quickmask==1], bins=np.arange(0.0,90.0,1.0))
	plt.xlim(0,90)
	plt.title('Difference between WMAPK and Planck30')
	plt.xlabel('Angle difference [deg]')
	plt.ylabel('Count')
	plt.plot(hist[1][:-1],hist[0])
	plt.savefig(outdirectory+'polang_diffplanck_wmap.png')
	plt.clf()

	# EOF