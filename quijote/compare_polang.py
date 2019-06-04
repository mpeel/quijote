#!/usr/bin/env python
# -*- coding: utf-8  -*-
#
# Do a quick analysis of the MFI maps
# 
# Version history:
#
# 31-May-2019  M. Peel       Started

import healpy as hp
import numpy as np
import matplotlib.pyplot as plt

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
commonmask = hp.read_map(outdirectory+'commonmask.fits',field=None)
commonmask = hp.ud_grade(commonmask,nside_out,order_in='RING',order_out='RING')

# Make a quick mask to get rid of low signal-to-noise pixels
threshold = 1.0#0.6

i = 4 # Use the 11GHz map for this
mapdata = hp.read_map(indirectory+maps[i],field=None)
qmap = hp.ud_grade(mapdata[1],nside_out,order_in='RING',order_out='RING')*commonmask
umap = hp.ud_grade(mapdata[2],nside_out,order_in='RING',order_out='RING')*commonmask
polang_map_11 = (0.5*np.arctan2(qmap[:],umap[:])* 180 / np.pi)
quickmask = np.zeros(npix_out)
quickmask[np.sqrt(qmap**2+umap**2) > threshold] = 1
# quickmask[qmap < -threshold] = 1
# quickmask[qmap > threshold] = 1
# quickmask[umap < -threshold] = 1
# quickmask[umap > threshold] = 1

mapdata = hp.read_map('/Users/mpeel/Documents/maps/wmap_planck_pol/weighted_both_debwk_feb2015_tqu.fits',field=None)
mapdata1 = hp.pixelfunc.reorder(-mapdata[1],n2r=True)
mapdata2 = hp.pixelfunc.reorder(-mapdata[2],n2r=True)
qmap = hp.ud_grade(mapdata1,nside_out,order_in='RING',order_out='RING')*commonmask
umap = hp.ud_grade(mapdata2,nside_out,order_in='RING',order_out='RING')*commonmask
hp.mollview(qmap,norm='hist')
plt.savefig(outdirectory+'polang_map_planckwmap_Q.pdf')
hp.mollview(umap,norm='hist')
plt.savefig(outdirectory+'polang_map_planckwmap_U.pdf')
polang_map_planckwmap = (0.5*np.arctan2(qmap[:],umap[:])* 180 / np.pi)

mapdata = hp.read_map('/Users/mpeel/Documents/maps/wmap9/wmap_band_smth_iqumap_r9_9yr_K_v5.fits',field=None)
# mapdata1 = hp.pixelfunc.reorder(mapdata[1],n2r=True)
# mapdata2 = hp.pixelfunc.reorder(mapdata[2],n2r=True)
qmap = hp.ud_grade(mapdata[1],nside_out,order_in='RING',order_out='RING')*commonmask
umap = hp.ud_grade(mapdata[2],nside_out,order_in='RING',order_out='RING')*commonmask
hp.mollview(qmap,norm='hist')
plt.savefig(outdirectory+'polang_map_wmapK_Q.pdf')
hp.mollview(umap,norm='hist')
plt.savefig(outdirectory+'polang_map_wmapK_U.pdf')
polang_map_wmap = (0.5*np.arctan2(qmap[:],umap[:])* 180 / np.pi)


mapdata = hp.read_map('/Users/mpeel/Documents/maps/planck2018/smooth_pol/512_60.00smoothed_PlanckR3fullbeam_28.4_1024_2018_mKCMBunits.fits',field=None)
# mapdata1 = hp.pixelfunc.reorder(mapdata[1],n2r=True)
# mapdata2 = hp.pixelfunc.reorder(mapdata[2],n2r=True)
qmap = hp.ud_grade(mapdata[1],nside_out,order_in='RING',order_out='RING')*commonmask
umap = hp.ud_grade(mapdata[2],nside_out,order_in='RING',order_out='RING')*commonmask
hp.mollview(qmap,norm='hist')
plt.savefig(outdirectory+'polang_map_planck30_Q.pdf')
hp.mollview(umap,norm='hist')
plt.savefig(outdirectory+'polang_map_planck30_U.pdf')
polang_map_planck30 = (0.5*np.arctan2(qmap[:],umap[:])* 180 / np.pi)


combhist = np.zeros((nummaps,89))
diff_to_11 = np.zeros(nummaps)
diff_to_planck = np.zeros(nummaps)
diff_to_wmap = np.zeros(nummaps)
diff_to_planckwmap = np.zeros(nummaps)
for i in range(2,nummaps):
	mapdata = hp.read_map(indirectory+maps[i],field=None)
	qmap = hp.ud_grade(mapdata[1],nside_out,order_in='RING',order_out='RING')*commonmask
	umap = hp.ud_grade(mapdata[2],nside_out,order_in='RING',order_out='RING')*commonmask
	polang_map = (0.5*np.arctan2(qmap[:],umap[:])* 180 / np.pi)

	print(np.median(polang_map[quickmask==1]))
	hist = np.histogram(polang_map[quickmask == 1], bins=np.arange(0.0,90.0,1.0))
	combhist[i][:] = hist[0][:]
	plt.xlim(0,90)
	plt.plot(hist[1][:-1],hist[0])
	plt.savefig(outdirectory+'polang_'+maps[i]+'.pdf')
	plt.clf()

	qmap[quickmask==0] = hp.UNSEEN
	umap[quickmask==0] = hp.UNSEEN
	hp.mollview(qmap,norm='hist')
	plt.savefig(outdirectory+'polang_map_'+maps[i]+'_Q.pdf')
	hp.mollview(umap,norm='hist')
	plt.savefig(outdirectory+'polang_map_'+maps[i]+'_U.pdf')
	plt.clf()

	if i == 2:
		polang_map_17 = polang_map.copy()
	elif i == 3:
		polang_map_19 = polang_map.copy()
	elif i == 6:
		print("Difference at 17GHz:")
		print(np.median(polang_map[quickmask==1]-polang_map_17[quickmask==1]))
		hist = np.histogram(polang_map[quickmask==1]-polang_map_17[quickmask==1], bins=np.arange(0.0,90.0,1.0))
		plt.xlim(0,90)
		plt.title('Difference between 17GHz channels')
		plt.xlabel('Angle difference [deg')
		plt.ylabel('Count')
		plt.plot(hist[1][:-1],hist[0])
		plt.savefig(outdirectory+'polang_diff17.png')
		plt.clf()
	elif i == 7:
		print("Difference at 19GHz:")
		print(np.median(polang_map[quickmask==1]-polang_map_19[quickmask==1]))
		hist = np.histogram(polang_map[quickmask==1]-polang_map_19[quickmask==1], bins=np.arange(0.0,90.0,1.0))
		plt.xlim(0,90)
		plt.title('Difference between 19GHz channels')
		plt.xlabel('Angle difference [deg')
		plt.ylabel('Count')
		plt.plot(hist[1][:-1],hist[0])
		plt.savefig(outdirectory+'polang_diff19.png')
		plt.clf()

	diff_to_11[i] = np.median(polang_map[quickmask==1]-polang_map_11[quickmask==1])

	if i != 4:
		hist = np.histogram(polang_map[quickmask==1]-polang_map_11[quickmask==1], bins=np.arange(0.0,90.0,1.0))
		plt.xlim(0,90)
		plt.title('Difference between 11GHz and ' + str(freqs[i]) + 'GHz (' + instrument[i] + ')')
		plt.xlabel('Angle difference [deg')
		plt.ylabel('Count')
		plt.plot(hist[1][:-1],hist[0])
		plt.savefig(outdirectory+'polang_diff11_'+maps[i]+'.png')
		plt.clf()

	diff_to_planckwmap[i] = np.median(polang_map[quickmask==1]-polang_map_planckwmap[quickmask==1])
	hist = np.histogram(polang_map[quickmask==1]-polang_map_planckwmap[quickmask==1], bins=np.arange(0.0,90.0,1.0))
	plt.xlim(0,90)
	plt.title('Difference between WMAP+Planck and ' + str(freqs[i]) + 'GHz (' + instrument[i] + ')')
	plt.xlabel('Angle difference [deg')
	plt.ylabel('Count')
	plt.plot(hist[1][:-1],hist[0])
	plt.savefig(outdirectory+'polang_diffplanckwmap_'+maps[i]+'.png')
	plt.clf()

	diff_to_wmap[i] = np.median(polang_map[quickmask==1]-polang_map_wmap[quickmask==1])
	hist = np.histogram(polang_map[quickmask==1]-polang_map_wmap[quickmask==1], bins=np.arange(0.0,90.0,1.0))
	plt.xlim(0,90)
	plt.title('Difference between WMAPK and ' + str(freqs[i]) + 'GHz (' + instrument[i] + ')')
	plt.xlabel('Angle difference [deg')
	plt.ylabel('Count')
	plt.plot(hist[1][:-1],hist[0])
	plt.savefig(outdirectory+'polang_diffwmap_'+maps[i]+'.png')
	plt.clf()

	diff_to_planck[i] = np.median(polang_map[quickmask==1]-polang_map_planck30[quickmask==1])
	hist = np.histogram(polang_map[quickmask==1]-polang_map_planck30[quickmask==1], bins=np.arange(0.0,90.0,1.0))
	plt.xlim(0,90)
	plt.title('Difference between Planck30 and ' + str(freqs[i]) + 'GHz (' + instrument[i] + ')')
	plt.xlabel('Angle difference [deg')
	plt.ylabel('Count')
	plt.plot(hist[1][:-1],hist[0])
	plt.savefig(outdirectory+'polang_diffplanck_'+maps[i]+'.png')
	plt.clf()

# print(combhist)
print(diff_to_11)
print(diff_to_planckwmap)
print(diff_to_planck)
print(diff_to_wmap)

print('Difference between WMAPK and Planck30')
print(np.median(polang_map_wmap[quickmask==1]-polang_map_planck30[quickmask==1]))
hist = np.histogram(polang_map_wmap[quickmask==1]-polang_map_planck30[quickmask==1], bins=np.arange(0.0,90.0,1.0))
plt.xlim(0,90)
plt.title('Difference between WMAPK and Planck30')
plt.xlabel('Angle difference [deg')
plt.ylabel('Count')
plt.plot(hist[1][:-1],hist[0])
plt.savefig(outdirectory+'polang_diffplanck_wmap.png')
plt.clf()


# combine_q /= 6.0
# combine_u /= 6.0

# hp.write_map(outdirectory+'commonmask.fits',commonmask,overwrite=True)
# hp.mollview(commonmask)
# plt.savefig(outdirectory+'commonmask.pdf')

# hp.write_map(outdirectory+'combine_q.fits',combine_q*commonmask,overwrite=True)
# hp.mollview(combine_q*commonmask,min=-1,max=1)
# plt.savefig(outdirectory+'combine_q.pdf')
# hp.write_map(outdirectory+'combine_u.fits',combine_u*commonmask,overwrite=True)
# hp.mollview(combine_u*commonmask,min=-1,max=1)
# plt.savefig(outdirectory+'combine_u.pdf')
# hp.write_map(outdirectory+'combine_P.fits',np.sqrt(combine_q**2+combine_u**2)*commonmask,overwrite=True)
# hp.mollview(np.sqrt(combine_q**2+combine_u**2)*commonmask,min=0,max=2.0,cmap=plt.get_cmap('jet'))
# plt.savefig(outdirectory+'combine_P.pdf')


# commonmask2 = hp.ud_grade(commonmask,256,order_in='RING',order_out='RING')

# hp.mollview(mapdata*1000.0*commonmask2,min=0,max=0.03,cmap=plt.get_cmap('jet'))
# plt.savefig(outdirectory+'combine_P_planck.pdf')
