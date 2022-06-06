#!/usr/bin/env python
# -*- coding: utf-8  -*-
#
# Do a quick test of polarised tt plots
#
# Version history:
#
# 29-Jul-2019  M. Peel       Started

import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from scipy import odr
from compare_polang import *
from smoothmap.smoothnoisemap import noiserealisation
from astrocode.polfunc import *

nside = 512
npix = hp.nside2npix(nside)
nside_out=64
npix_out=hp.nside2npix(nside_out)

indirectory='/Users/mpeel/Documents/maps/quijote_201905/smooth/'
outdirectory='/Users/mpeel/Documents/maps/quijote_201905/ttplots/'
map_11 = '512_60.00smoothed_mfi3_11.0_512_201905_mKCMBunits.fits'
map_11_qvar = '512_60.00smoothed_mfi3_11.0_512_201905_weight_1_variance.fits'
map_11_uvar = '512_60.00smoothed_mfi3_11.0_512_201905_weight_2_variance.fits'
map_17_qvar = '512_60.00smoothed_mfi2_17.0_512_201905_weight_1_variance.fits'
map_17_uvar = '512_60.00smoothed_mfi2_17.0_512_201905_weight_2_variance.fits'
map_11_std = '../std_H3_11_sm1deg_nside64.fits'
map_17_std = '../std_H4_17_sm1deg_nside64.fits'
map_11_noise = '../quijote_11GHz_horn3_0001_sm1deg.fits'
map_17_noise = '../quijote_17GHz_horn4_0001_sm1deg.fits'
rescalenoise = 1.0
index = 3.0
usenoisemap = False
extramask_filename = ''
# extramask_filename = '../../quijote_masks/fan_mask_extended_final_nside64.fits'
# extramask_filename = '../../quijote_masks/NPS_mask_ns256.fits'

# commonmask = hp.read_map(indirectory+'../analyse/mfi_commonmask.fits',field=None)
commonmask = hp.read_map(indirectory+'../../quijote_masks/mask_quijote_ncp_lowdec_nside512.fits',field=None)
commonmask = hp.ud_grade(commonmask,nside_out,order_in='RING',order_out='RING')

if extramask_filename != '':
	extramask = hp.read_map(indirectory+extramask_filename,field=None)
	extramask = hp.ud_grade(extramask,nside_out,order_in='RING',order_out='RING')
	commonmask = commonmask * extramask

snmask = commonmask.copy()
mapdata = hp.read_map(indirectory+map_11,field=None)
qmap_311 = hp.ud_grade(mapdata[1],nside_out,order_in='RING',order_out='RING')*commonmask
umap_311 = hp.ud_grade(mapdata[2],nside_out,order_in='RING',order_out='RING')*commonmask
if usenoisemap:
	mapdata = hp.read_map(indirectory+map_11_std,field=None)
	qmap_311_var = hp.ud_grade(mapdata[1],nside_out,order_in='RING',order_out='RING')*commonmask
	umap_311_var = hp.ud_grade(mapdata[2],nside_out,order_in='RING',order_out='RING')*commonmask
	mapdata = hp.read_map(indirectory+map_11_noise,field=None)
	noise_q = hp.ud_grade(mapdata[1],nside_out,order_in='RING',order_out='RING')*commonmask
	noise_u = hp.ud_grade(mapdata[2],nside_out,order_in='RING',order_out='RING')*commonmask
else:
	mapdata = hp.read_map(indirectory+map_11_qvar,field=None)
	qmap_311_var = np.sqrt(hp.ud_grade(mapdata,nside_out,order_in='RING',order_out='RING')*commonmask)
	mapdata = hp.read_map(indirectory+map_11_uvar,field=None)
	umap_311_var = np.sqrt(hp.ud_grade(mapdata,nside_out,order_in='RING',order_out='RING')*commonmask)
	noise_q = noiserealisation(qmap_311_var, npix_out)
	noise_u = noiserealisation(umap_311_var, npix_out)
hp.mollview(qmap_311_var)
plt.savefig(outdirectory+'Q_311_var.png')
plt.clf()
hp.mollview(umap_311_var)
plt.savefig(outdirectory+'U_311_var.png')
plt.clf()

qmap_317 = qmap_311.copy() * (11.0/17.0)**index
umap_317 = umap_311.copy() * (11.0/17.0)**index


if usenoisemap:
	mapdata = hp.read_map(indirectory+map_17_std,field=None)
	qmap_317_var = hp.ud_grade(mapdata[1],nside_out,order_in='RING',order_out='RING')*commonmask
	umap_317_var = hp.ud_grade(mapdata[2],nside_out,order_in='RING',order_out='RING')*commonmask

	mapdata = hp.read_map(indirectory+map_17_noise,field=None)
	noise_q_317 = hp.ud_grade(mapdata[1],nside_out,order_in='RING',order_out='RING')*commonmask
	noise_u_317 = hp.ud_grade(mapdata[2],nside_out,order_in='RING',order_out='RING')*commonmask
else:
	mapdata = hp.read_map(indirectory+map_17_qvar,field=None)
	qmap_317_var = np.sqrt(hp.ud_grade(mapdata,nside_out,order_in='RING',order_out='RING')*commonmask)
	mapdata = hp.read_map(indirectory+map_17_uvar,field=None)
	umap_317_var = np.sqrt(hp.ud_grade(mapdata,nside_out,order_in='RING',order_out='RING')*commonmask)

	# qmap_317_var = qmap_311_var * (11.0/17.0)**index
	# umap_317_var = umap_311_var * (11.0/17.0)**index
	noise_q_317 = noiserealisation(qmap_317_var, npix_out)
	noise_u_317 = noiserealisation(umap_317_var, npix_out)

hp.mollview(qmap_317_var)
plt.savefig(outdirectory+'Q_317_var.png')
plt.clf()
hp.mollview(umap_317_var)
plt.savefig(outdirectory+'U_317_var.png')
plt.clf()

# Add noise to the maps
qmap_311 = qmap_311+rescalenoise*noise_q
umap_311 = umap_311+rescalenoise*noise_u
qmap_317 = qmap_317+rescalenoise*noise_q_317
umap_317 = umap_317+rescalenoise*noise_u_317

P_311 = np.sqrt(qmap_311**2+umap_311**2)
P_311_sigma = np.sqrt(((qmap_311*rescalenoise*qmap_311_var)**2+(umap_311*rescalenoise*umap_311_var)**2)/((qmap_311)**2+(umap_311)**2))
commonmask[P_311 == hp.UNSEEN] = 0
commonmask[~np.isfinite(P_311)] = 0

P_317 = np.sqrt(qmap_317**2+umap_317**2)
P_317_sigma = np.sqrt(((qmap_317*rescalenoise*qmap_317_var)**2+(umap_317*rescalenoise*umap_317_var)**2)/((qmap_317)**2+(umap_317)**2))
commonmask[P_317 == hp.UNSEEN] = 0
commonmask[~np.isfinite(P_317)] = 0

P_311[commonmask==0] = hp.UNSEEN
P_311_sigma[commonmask==0] = hp.UNSEEN
hp.mollview(P_311)
plt.savefig(outdirectory+'P_11.png')
plt.clf()
P_311_deb_as = [P_311.copy(), P_311.copy()]
P_311_deb_as[0][commonmask==1], P_311_deb_as[1][commonmask==1] = debias_p_as(qmap_311[commonmask==1], umap_311[commonmask==1], rescalenoise*qmap_311_var[commonmask==1], rescalenoise*umap_311_var[commonmask==1])
hp.mollview(P_311_deb_as[0])
plt.savefig(outdirectory+'P_11_deb_as.png')
plt.clf()
P_311_deb_mas = [P_311.copy(), P_311.copy()]
P_311_deb_mas[0][commonmask==1], P_311_deb_mas[1][commonmask==1] = debias_p_mas(qmap_311[commonmask==1], umap_311[commonmask==1], rescalenoise*qmap_311_var[commonmask==1], rescalenoise*umap_311_var[commonmask==1])
hp.mollview(P_311_deb_mas[0])
plt.savefig(outdirectory+'P_11_deb_mas.png')
plt.clf()

P_317[commonmask==0] = hp.UNSEEN
P_317_sigma[commonmask==0] = hp.UNSEEN
hp.mollview(P_317)
plt.savefig(outdirectory+'P_17.png')
plt.clf()
P_317_deb_as = [P_317.copy(), P_317.copy()]
P_317_deb_as[0][commonmask==1], P_317_deb_as[1][commonmask==1] = debias_p_as(qmap_317[commonmask==1], umap_317[commonmask==1], rescalenoise*qmap_317_var[commonmask==1], rescalenoise*umap_317_var[commonmask==1])
hp.mollview(P_317_deb_as[0])
plt.savefig(outdirectory+'P_17_deb_as.png')
plt.clf()
P_317_deb_mas = [P_317.copy(), P_317.copy()]
P_317_deb_mas[0][commonmask==1], P_317_deb_mas[1][commonmask==1] = debias_p_mas(qmap_317[commonmask==1], umap_317[commonmask==1], rescalenoise*qmap_317_var[commonmask==1], rescalenoise*umap_317_var[commonmask==1])
hp.mollview(P_317_deb_mas[0])
plt.savefig(outdirectory+'P_17_deb_mas.png')
plt.clf()

# Cut out bright pixels
# commonmask[P_311 > 2.0] = 0

print('## NO SN CUT')
print('P, biased')
plot_tt(P_311[commonmask==1],P_317[commonmask==1],outdirectory+'ttplot_11_17.png',freq1=11.0,freq2=17.0,sigma=P_317_sigma[commonmask==1],sigma_x=P_311_sigma[commonmask==1])
plt.clf()

print('P, deb_as')
plot_tt(P_311_deb_as[0][commonmask==1],P_317_deb_as[0][commonmask==1],outdirectory+'ttplot_11_17_deb_as.png',freq1=11.0,freq2=17.0,sigma=P_317_deb_as[1][commonmask==1],sigma_x=P_311_deb_as[1][commonmask==1])
plt.clf()

print('P, deb_mas')
plot_tt(P_311_deb_mas[0][commonmask==1],P_317_deb_mas[0][commonmask==1],outdirectory+'ttplot_11_17_deb_mas.png',freq1=11.0,freq2=17.0,sigma=P_317_deb_mas[1][commonmask==1],sigma_x=P_311_deb_mas[1][commonmask==1])
plt.clf()

print('Q')
plot_tt(qmap_311[commonmask==1],qmap_317[commonmask==1],outdirectory+'ttplot_11_17_Q.png',freq1=11.0,freq2=17.0,sigma=rescalenoise*qmap_317_var[commonmask==1],sigma_x=rescalenoise*qmap_311_var[commonmask==1])
plt.clf()
print('U')
plot_tt(umap_311[commonmask==1],umap_317[commonmask==1],outdirectory+'ttplot_11_17_U.png',freq1=11.0,freq2=17.0,sigma=rescalenoise*umap_317_var[commonmask==1],sigma_x=rescalenoise*umap_311_var[commonmask==1])
plt.clf()

print('## SN 1.0 CUT')
sigma = 1.0
snmask[np.abs(qmap_311) < sigma * rescalenoise * qmap_311_var] = 0.0
snmask[np.abs(umap_311) < sigma * rescalenoise * umap_311_var] = 0.0
snmask[np.abs(qmap_317) < sigma * rescalenoise * qmap_317_var] = 0.0
snmask[np.abs(umap_317) < sigma * rescalenoise * umap_317_var] = 0.0

print('P, biased')
plot_tt(P_311[snmask==1],P_317[snmask==1],outdirectory+'ttplot_11_17_SN.png',freq1=11.0,freq2=17.0,sigma=P_317_sigma[snmask==1],sigma_x=P_311_sigma[snmask==1])
plt.clf()
print('P, deb_as')
plot_tt(P_311_deb_as[0][snmask==1],P_317_deb_as[0][snmask==1],outdirectory+'ttplot_11_17_SN_deb_as.png',freq1=11.0,freq2=17.0,sigma=P_317_deb_as[1][snmask==1],sigma_x=P_311_deb_as[1][snmask==1])
plt.clf()

print('P, deb_mas')
plot_tt(P_311_deb_mas[0][snmask==1],P_317_deb_mas[0][snmask==1],outdirectory+'ttplot_11_17_SN_deb_as.png',freq1=11.0,freq2=17.0,sigma=P_317_deb_mas[1][snmask==1],sigma_x=P_311_deb_mas[1][snmask==1])
plt.clf()

print('Q')
plot_tt(qmap_311[snmask==1],qmap_317[snmask==1],outdirectory+'ttplot_11_17_Q_SN.png',freq1=11.0,freq2=17.0,sigma=rescalenoise*qmap_317_var[snmask==1],sigma_x=rescalenoise*qmap_311_var[snmask==1])
plt.clf()
print('U')
plot_tt(umap_311[snmask==1],umap_317[snmask==1],outdirectory+'ttplot_11_17_U_SN.png',freq1=11.0,freq2=17.0,sigma=rescalenoise*umap_317_var[snmask==1],sigma_x=rescalenoise*umap_311_var[snmask==1])
plt.clf()

P_311[snmask==0] = hp.UNSEEN
hp.mollview(P_311)
plt.savefig(outdirectory+'P_11_cut.png')
plt.clf()

P_317[snmask==0] = hp.UNSEEN
hp.mollview(P_317)
plt.savefig(outdirectory+'P_17_cut.png')
plt.clf()

print('## SN 3.8 CUT')
sigma = 3.8
snmask[np.abs(qmap_311) < sigma * rescalenoise * qmap_311_var] = 0.0
snmask[np.abs(umap_311) < sigma * rescalenoise * umap_311_var] = 0.0
snmask[np.abs(qmap_317) < sigma * rescalenoise * qmap_317_var] = 0.0
snmask[np.abs(umap_317) < sigma * rescalenoise * umap_317_var] = 0.0

print('P, biased')
plot_tt(P_311[snmask==1],P_317[snmask==1],outdirectory+'ttplot_11_17_SN10.png',freq1=11.0,freq2=17.0,sigma=P_317_sigma[snmask==1],sigma_x=P_311_sigma[snmask==1])
plt.clf()
print('P, deb_as')
plot_tt(P_311_deb_as[0][snmask==1],P_317_deb_as[0][snmask==1],outdirectory+'ttplot_11_17_SN10_deb_as.png',freq1=11.0,freq2=17.0,sigma=P_317_deb_as[1][snmask==1],sigma_x=P_311_deb_as[1][snmask==1])
plt.clf()

print('P, deb_mas')
plot_tt(P_311_deb_mas[0][snmask==1],P_317_deb_mas[0][snmask==1],outdirectory+'ttplot_11_17_SN10_deb_as.png',freq1=11.0,freq2=17.0,sigma=P_317_deb_mas[1][snmask==1],sigma_x=P_311_deb_mas[1][snmask==1])
plt.clf()

print('Q')
plot_tt(qmap_311[snmask==1],qmap_317[snmask==1],outdirectory+'ttplot_11_17_Q_SN10.png',freq1=11.0,freq2=17.0,sigma=rescalenoise*qmap_317_var[snmask==1],sigma_x=rescalenoise*qmap_311_var[snmask==1])
plt.clf()
print('U')
plot_tt(umap_311[snmask==1],umap_317[snmask==1],outdirectory+'ttplot_11_17_U_SN10.png',freq1=11.0,freq2=17.0,sigma=rescalenoise*umap_317_var[snmask==1],sigma_x=rescalenoise*umap_311_var[snmask==1])
plt.clf()

# EOF
