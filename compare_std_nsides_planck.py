#!/usr/bin/env python
# -*- coding: utf-8  -*-
#
# Do a quick comparison of the simulated standard deviations with different nsides
# 
# Version history:
#
# 03-Feb-2020  M. Peel       Started

import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
indirectory = '/Users/mpeel/Documents/maps/wmap9_planck2018_tqu_noise/'
outdirectory = indirectory

map_512 = hp.read_map(indirectory+'60.00smoothed_PlanckR3fullbeam_28.4_1024_2018_mKCMBunits_7_variance_512.fits',field=None)
map_256 = hp.read_map(indirectory+'60.00smoothed_PlanckR3fullbeam_28.4_1024_2018_mKCMBunits_7_variance_256.fits',field=None)
map_64 = hp.read_map(indirectory+'60.00smoothed_PlanckR3fullbeam_28.4_1024_2018_mKCMBunits_7_variance_64.fits',field=None)

map_512 = np.sqrt(map_512)
map_256 = np.sqrt(map_256)
map_64 = np.sqrt(map_64)

mask_512 = np.ones(12*512*512)
mask_256 = hp.ud_grade(mask_512,nside_out=256)
mask_64 = hp.ud_grade(mask_512,nside_out=64)

power=0
map_512_64 = np.sqrt(hp.ud_grade(map_512**2, nside_out=64,power=power))
map_512_256 = np.sqrt(hp.ud_grade(map_512**2, nside_out=256,power=power))
map_256_64 = np.sqrt(hp.ud_grade(map_256**2, nside_out=64,power=power))

hp.mollview(map_512)
plt.savefig(outdirectory+'std_512.pdf')
hp.mollview(map_256)
plt.savefig(outdirectory+'std_256.pdf')
hp.mollview(map_64)
plt.savefig(outdirectory+'std_64.pdf')

hp.mollview(map_512_64)
plt.savefig(outdirectory+'std_512_64.pdf')

hp.mollview(map_512_256)
plt.savefig(outdirectory+'std_512_256.pdf')

hp.mollview(map_256_64)
plt.savefig(outdirectory+'std_256_64.pdf')

print('std')
print(np.std(map_512[mask_512==1]))
print(np.std(map_256[mask_256==1]))
print(np.std(map_64[mask_64==1]))

print('median')
print(np.median(map_512[mask_512==1]))
print(np.median(map_256[mask_256==1]))
print(np.median(map_64[mask_64==1]))

print('512_to_64')
print(np.std(map_512_64[mask_64==1]))
print(np.median(map_512_64[mask_64==1]))
print('512_to_256')
print(np.std(map_512_256[mask_256==1]))
print(np.median(map_512_256[mask_256==1]))

print('256_to_64')
print(np.std(map_256_64[mask_64==1]))
print(np.median(map_256_64[mask_64==1]))
