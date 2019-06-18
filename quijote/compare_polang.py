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

def linfit(x, param):
	return param[0]*x+param[1]

def compute_residuals_linfit(param, x, y):
	model = linfit(x, param)
	residual = y - model
	return residual

def linfit2(x, A, B):
	return A*x+B

def calc_polang(Q, U):
	return 0.5*np.arctan2(Q, U) * 180 / np.pi

def calc_polang_unc(Q, U, Qerr, Uerr):
	unc_map = np.sqrt((Qerr**2)*(-0.5*U/(Q**2.0+U**2.0))**2.0 + (Uerr**2)*(0.5*Q/(Q**2.0+U**2.0))**2.0) * 180 / np.pi
	unc_map[~np.isfinite(unc_map)] = 1000.0
	return unc_map

def plot_tt(vals1,vals2,outputname,sigma=np.empty(0),leastsq=False):
	# print(sigma)
	# Do a fit
	params = [0.5,0]
	if leastsq:
		# least-squares fit without uncertainties
		param_est, cov_x, infodict, mesg_result, ret_value = optimize.leastsq(compute_residuals_linfit, params, args=(vals1, vals2),full_output=True)
	elif sigma.size:
		# Do a curve fit to use the uncertainties
		param_est, cov_x = optimize.curve_fit(linfit2, vals1, vals2, params, sigma=sigma)
	else:
		# Do a curve fit without uncertainties
		param_est, cov_x = optimize.curve_fit(linfit2, vals1, vals2, params)
	sigma_param_est = np.sqrt(np.diagonal(cov_x))
	
	#Plot the data and the results
	if sigma.size:
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
		plt.plot(xvals,linfit(xvals,param_est),'g',label="Fit: " + mesg_fit)
	else:
		plt.plot(xvals,linfit(xvals,param_est),'g')
	plt.legend(prop={'size':8})
	plt.savefig(outputname)
	plt.clf()
	return [param_est, sigma_param_est]

nside = 512
npix = hp.nside2npix(nside)
nside_out=64
npix_out=hp.nside2npix(nside_out)
maps = [str(nside)+'_60.00smoothed_mfi1_11.0_512_201905_mKCMBunits.fits',str(nside)+'_60.00smoothed_mfi1_13.0_512_201905_mKCMBunits.fits',str(nside)+'_60.0smoothed_mfi2_17.0_512_201905_mKCMBunits.fits',str(nside)+'_60.0smoothed_mfi2_19.0_512_201905_mKCMBunits.fits',str(nside)+'_60.00smoothed_mfi3_11.0_512_201905_mKCMBunits.fits',str(nside)+'_60.00smoothed_mfi3_13.0_512_201905_mKCMBunits.fits',str(nside)+'_60.0smoothed_mfi4_17.0_512_201905_mKCMBunits.fits',str(nside)+'_60.0smoothed_mfi4_19.0_512_201905_mKCMBunits.fits']

indirectory = '/Users/mpeel/Documents/maps/quijote_201905/smooth/'
outdirectory = '/Users/mpeel/Documents/maps/quijote_201905/analyse/'

nummaps = len(maps)
freqs = [11,13,17,19,11,13,17,19]
instrument = ['MFI1','MFI1','MFI2','MFI2','MFI3','MFI3','MFI4','MFI4']
normfreq = 28.4
index = 3.0
commonmask = hp.read_map(outdirectory+'mfi_commonmask.fits',field=None)
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
var_q_311 = hp.read_map(indirectory+maps[i].replace('60.0s','60.00s').replace('_mKCMBunits','_weight_1_variance'),field=None)
var_u_311 = hp.read_map(indirectory+maps[i].replace('60.0s','60.00s').replace('_mKCMBunits','_weight_2_variance'),field=None)
var_q_311 = hp.ud_grade(var_q_311, nside_out, power=0)
var_u_311 = hp.ud_grade(var_u_311, nside_out, power=0)

var_q_311[var_q_311 < 0.0] = 10000.0
var_u_311[var_u_311 < 0.0] = 10000.0
var_q_311[var_q_311==0.0] = 10000.0
var_u_311[var_u_311==0.0] = 10000.0
var_q_311[~np.isfinite(var_q_311)] = 10000.0
var_u_311[~np.isfinite(var_u_311)] = 10000.0


polang_map_11 = calc_polang(qmap_311[:],umap_311[:])
polang_unc_map_11 = calc_polang_unc(qmap_311,umap_311,np.sqrt(var_q_311),np.sqrt(var_u_311))
hp.mollview(polang_map_11,cmap=plt.get_cmap('jet'))
plt.savefig(outdirectory+'polang_map_mfi11.png')
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
quickmask[commonmask == 0] = 0

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
polang_map_planckwmap = (0.5*np.arctan2(qmap[:],umap[:])* 180 / np.pi)

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
polang_map_wmap = (0.5*np.arctan2(qmap_wmap[:],umap_wmap[:])* 180 / np.pi)

threshold3 = 2.5
quickmask[qmap_wmap*(11.0/23.0)**-3.0 < -threshold3] = 0
quickmask[qmap_wmap*(11.0/23.0)**-3.0 > threshold3] = 0
quickmask[umap_wmap*(11.0/23.0)**-3.0 < -threshold3] = 0
quickmask[umap_wmap*(11.0/23.0)**-3.0 > threshold3] = 0

offset_values_Q = np.zeros(nummaps)
offset_values_U = np.zeros(nummaps)


mapdata = hp.read_map('/Users/mpeel/Documents/maps/planck2018/smooth_pol/512_60.00smoothed_PlanckR3fullbeam_28.4_1024_2018_mKCMBunits.fits',field=None)
# mapdata1 = hp.pixelfunc.reorder(mapdata[1],n2r=True)
# mapdata2 = hp.pixelfunc.reorder(mapdata[2],n2r=True)
qmap_planck = hp.ud_grade(mapdata[1],nside_out,order_in='RING',order_out='RING')*commonmask
umap_planck = hp.ud_grade(mapdata[2],nside_out,order_in='RING',order_out='RING')*commonmask
hp.mollview(qmap_planck,norm='hist')
plt.savefig(outdirectory+'polang_map_planck30_Q.png')
hp.mollview(umap_planck,norm='hist')
plt.savefig(outdirectory+'polang_map_planck30_U.png')
plt.clf()
polang_map_planck30 = (0.5*np.arctan2(qmap_planck[:],umap_planck[:])* 180 / np.pi)

plot_tt(qmap_wmap[quickmask==1],qmap_planck[quickmask==1],outdirectory+'tt_wmap_planck_Q.png')
plot_tt(umap_wmap[quickmask==1],umap_planck[quickmask==1],outdirectory+'tt_wmap_planck_U.png')

# Get the offset between 11GHz and Planck/WMAP
# fit,fiterr=plot_tt(qmap_planck[quickmask==1]*(11.0/28.4)**-3.0,qmap_311[quickmask==1],outdirectory+'tt_311_planck_Q.png')
# print(fit,fiterr)
fit,fiterr=plot_tt(qmap_planck[quickmask==1]*(11.0/28.4)**-3.0,qmap_311[quickmask==1],outdirectory+'tt_311_planck_Q_sigma.png',sigma=np.sqrt(var_q_311[quickmask==1]))
print(fit,fiterr)

# fit,fiterr=plot_tt(umap_planck[quickmask==1]*(11.0/28.4)**-3.0,umap_311[quickmask==1],outdirectory+'tt_311_planck_U.png')
# print(fit,fiterr)
fit,fiterr=plot_tt(umap_planck[quickmask==1]*(11.0/28.4)**-3.0,umap_311[quickmask==1],outdirectory+'tt_311_planck_U_sigma.png',np.sqrt(var_u_311[quickmask==1]))
print(fit,fiterr)

# fit,fiterr=plot_tt(qmap_wmap[quickmask==1]*(11.0/23.0)**-3.0,qmap_311[quickmask==1],outdirectory+'tt_311_wmap_Q.png')
# print(fit,fiterr)
fit,fiterr=plot_tt(qmap_wmap[quickmask==1]*(11.0/23.0)**-3.0,qmap_311[quickmask==1],outdirectory+'tt_311_wmap_Q_sigma.png',np.sqrt(var_q_311[quickmask==1]))
print(fit,fiterr)
offset_values_Q[i] = fit[1]
if applyoffsets:
	qmap_311 = qmap_311 - fit[1]
# fit,fiterr=plot_tt(umap_wmap[quickmask==1]*(11.0/23.0)**-3.0,umap_311[quickmask==1],outdirectory+'tt_311_wmap_U.png')
# print(fit,fiterr)
fit,fiterr=plot_tt(umap_wmap[quickmask==1]*(11.0/23.0)**-3.0,umap_311[quickmask==1],outdirectory+'tt_311_wmap_U_sigma.png',np.sqrt(var_u_311[quickmask==1]))
print(fit,fiterr)
offset_values_U[i] = fit[1]
if applyoffsets:
	umap_311 = umap_311 - fit[1]
polang_map_11 = (0.5*np.arctan2(qmap_311[:],umap_311[:])* 180 / np.pi)


combhist = np.zeros((nummaps,89))
diff_to_11 = np.zeros(nummaps)
diff_to_planck = np.zeros(nummaps)
diff_to_wmap = np.zeros(nummaps)
diff_to_planckwmap = np.zeros(nummaps)
diff_to_11_2 = np.zeros(nummaps)
diff_to_planck_2 = np.zeros(nummaps)
diff_to_wmap_2 = np.zeros(nummaps)
diff_to_planckwmap_2 = np.zeros(nummaps)
for i in range(2,nummaps):
	print(maps[i])
	mapdata = hp.read_map(indirectory+maps[i],field=None)
	qmap = hp.ud_grade(mapdata[1],nside_out,order_in='RING',order_out='RING')*commonmask
	umap = hp.ud_grade(mapdata[2],nside_out,order_in='RING',order_out='RING')*commonmask

	# Get the variance maps
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
	print(min(np.sqrt(var_q[quickmask==1])))
	print(max(np.sqrt(var_q[quickmask==1])))
	print(min(np.sqrt(var_u[quickmask==1])))
	print(max(np.sqrt(var_u[quickmask==1])))
	# Check for offsets vs. the 11GHz maps
	fit,fiterr=plot_tt(qmap_311[quickmask==1],qmap[quickmask==1],outdirectory+'tt_311_'+maps[i]+'_Q.pdf',np.sqrt(var_q[quickmask==1]))
	print(fit,fiterr)
	if applyoffsets:
		qmap = qmap-fit[1]
	if i != 4:
		offset_values_Q[i] = fit[1]

	fit,fiterr=plot_tt(umap_311[quickmask==1],umap[quickmask==1],outdirectory+'tt_311_'+maps[i]+'_U.pdf',np.sqrt(var_u[quickmask==1]))
	print(fit,fiterr)
	if applyoffsets:
		umap = umap-fit[1]
	if i != 4:
		offset_values_U[i] = fit[1]
	polang_map = calc_polang(qmap,umap)
	polang_unc_map = calc_polang_unc(qmap,umap,np.sqrt(var_q),np.sqrt(var_u))

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

	qmap[quickmask==0] = hp.UNSEEN
	umap[quickmask==0] = hp.UNSEEN
	hp.mollview(qmap,norm='hist')
	plt.savefig(outdirectory+'polang_map_'+maps[i]+'_Q.pdf')
	hp.mollview(umap,norm='hist')
	plt.savefig(outdirectory+'polang_map_'+maps[i]+'_U.pdf')
	plt.clf()


	if i == 2:
		polang_map_17 = polang_map.copy()
		polang_unc_map_17 = polang_unc_map.copy()
	elif i == 3:
		polang_map_19 = polang_map.copy()
		polang_unc_map_19 = polang_unc_map.copy()
	elif i == 5:
		polang_diff_temp = polang_map - polang_map_11
		polang_diff_temp[quickmask==0] = hp.UNSEEN
		hp.gnomview(polang_diff_temp,reso=15.0,min=-20,max=10,cmap=plt.get_cmap('jet'))
		plt.savefig(outdirectory+'galcen_diff_11_13.pdf')

	elif i == 6:
		print("Difference at 17GHz:")
		print(np.median(polang_map[quickmask==1]-polang_map_17[quickmask==1]))
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
	# print(diff_to_11[i])
	diff_to_11_2[i] = np.average(polang_map[quickmask==1]-polang_map_11[quickmask==1],weights=1.0/(polang_unc_map[quickmask==1]**2.0))

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
	# print(diff_to_planckwmap[i])
	diff_to_planckwmap_2[i] = np.average(polang_map[quickmask==1]-polang_map_planckwmap[quickmask==1],weights=1.0/(polang_unc_map[quickmask==1]**2.0))
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
	# print(diff_to_wmap[i])
	diff_to_wmap_2[i] = np.average(polang_map[quickmask==1]-polang_map_wmap[quickmask==1],weights=1.0/(polang_unc_map[quickmask==1]**2.0))
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
	# print(diff_to_planck[i])
	diff_to_planck_2[i] = np.average(polang_map[quickmask==1]-polang_map_planck30[quickmask==1],weights=1.0/(polang_unc_map[quickmask==1]**2.0))
	hist = np.histogram(polang_map[quickmask==1]-polang_map_planck30[quickmask==1], bins=np.arange(0.0,90.0,1.0))
	plt.xlim(0,90)
	plt.title('Difference between Planck30 and ' + str(freqs[i]) + 'GHz (' + instrument[i] + ')')
	plt.xlabel('Angle difference [deg]')
	plt.ylabel('Count')
	plt.plot(hist[1][:-1],hist[0])
	plt.savefig(outdirectory+'polang_diffplanck_'+maps[i]+'.png')
	plt.clf()

# print(combhist)
np.set_printoptions(formatter={'float': '{: 0.5f}'.format})
print(offset_values_Q[2:])
print(offset_values_U[2:])
np.set_printoptions(formatter={'float': '{: 0.3f}'.format})
print(diff_to_11[2:])
print(diff_to_planckwmap[2:])
print(diff_to_planck[2:])
print(diff_to_wmap[2:])
print(diff_to_11_2[2:])
print(diff_to_planckwmap_2[2:])
print(diff_to_planck_2[2:])
print(diff_to_wmap_2[2:])

print('Difference between WMAPK and Planck30')
print(np.median(polang_map_wmap[quickmask==1]-polang_map_planck30[quickmask==1]))
hist = np.histogram(polang_map_wmap[quickmask==1]-polang_map_planck30[quickmask==1], bins=np.arange(0.0,90.0,1.0))
plt.xlim(0,90)
plt.title('Difference between WMAPK and Planck30')
plt.xlabel('Angle difference [deg]')
plt.ylabel('Count')
plt.plot(hist[1][:-1],hist[0])
plt.savefig(outdirectory+'polang_diffplanck_wmap.png')
plt.clf()

# EOF