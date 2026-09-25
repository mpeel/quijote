import healpy as hp
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import numpy as np

def conv_nobs_variance_map(inputmap, sigma_0):
	newmap = sigma_0**2 / inputmap
	return newmap

# QUIJOTE MFI11
orig_map_filename = '/Volumes/Main5TB/maps/quijote_202103/reform/mfi_mar2021_11.0_3.fits'
smoothed_map_filename = '/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/512_60.0smoothed_QUIJOTEMFI3_11.0_2021_mKCMBunits.fits'
smoothed_map_64_filename = '/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/64_60.0smoothed_QUIJOTEMFI3_11.0_2021_mKCMBunits.fits'
rescale_orig_map = 1.0
sigma_P = 0

# WMAP K
orig_map_filename = '/Volumes/Main5TB/maps/wmap9/wmap_band_iqumap_r9_9yr_K_v5.fits'
smoothed_map_filename = '/Users/mpeel/Documents/maps/wmap9_tqu_v1.5_noise_v1.0_10k/512_60.0smoothed_wmap9beamNoise_22.8_512_2013_mKCMBunits.fits'
smoothed_map_64_filename = '/Users/mpeel/Documents/maps/wmap9_tqu_v1.5_noise_v1.0_10k/64_60.0smoothed_wmap9beamNoise_22.8_512_2013_mKCMBunits.fits'
rescale_orig_map = 1.0
sigma_P = 1.435

# Planck 30
# orig_map_filename = '/Volumes/Main5TB/maps/planck2020/LFI_SkyMap_030_1024_R4.00_full.fits'
# smoothed_map_filename = '/Users/mpeel/Documents/maps/planck2020_tqu_v1.5_noise_v1.0_10k/512_60.0smoothed_PlanckR4fullbeamnodpNoise_28.4_1024_2020_mKCMBunits.fits'
# smoothed_map_64_filename = '/Users/mpeel/Documents/maps/planck2020_tqu_v1.5_noise_v1.0_10k/64_60.0smoothed_PlanckR4fullbeamnodpNoise_28.4_1024_2020_mKCMBunits.fits'
# rescale_orig_map = 1e3
# sigma_P = 0

orig_map, orig_map_hdr = hp.read_map(orig_map_filename,field=None,h=True)
smoothed_map, smoothed_map_hdr = hp.read_map(smoothed_map_filename,field=None,h=True)
smoothed_map_64, smoothed_map_64_hdr = hp.read_map(smoothed_map_64_filename,field=None,h=True)
print(len(orig_map))
print(len(smoothed_map))
print(orig_map_hdr)
print(smoothed_map_hdr)

# inputfits = fits.open(orig_map_filename)
# cols = inputfits[1].columns
# col_names = cols.names
# nmaps = len(cols)
# maps = []
# for i in range(0,nmaps):
# 	maps.append(inputfits[1].data.field(i))
# 	print(len(maps[i]))
# print(len(maps[0]))
# # Check to see whether we have nested data, and switch to ring if that is the case.
# try:
# 	if (inputfits[1].header['ORDERING'] == 'NESTED'):
# 		maps = hp.reorder(maps,n2r=True)
# except:
# 	null = 0
# newheader = inputfits[1].header.copy(strip=False)
# inputfits.close()
# orig_map = maps

if sigma_P == 0:
	orig_map_64_qvar = hp.ud_grade(orig_map[4]*rescale_orig_map*rescale_orig_map,nside_out=64,power=2)
else:
	orig_map_64_qvar = hp.ud_grade(conv_nobs_variance_map(orig_map[3],sigma_P),nside_out=64,power=2)


# for i in range(0,5):
# 	print(str(i) + "	" + str(orig_map[4][i]) + "	" + str(smoothed_map[4][i]) + "	" + str(smoothed_map[4][i]/orig_map[4][i]))
# print(np.median(smoothed_map[4]/orig_map[4]))

for i in range(0,5):
	print(str(i) + "	" + str(orig_map_64_qvar[i]) + "	" + str(smoothed_map_64[4][i]) + "	" + str(smoothed_map[4][i]/orig_map_64_qvar[i]) + "	" + str(((64.0/512.0)**2)*smoothed_map[4][i]/orig_map_64_qvar[i]))
print(np.median(smoothed_map_64[4]/orig_map_64_qvar))
print(((64.0/512.0)**2)*np.median(smoothed_map_64[4]/orig_map_64_qvar))

hp.mollview(orig_map_64_qvar)
plt.savefig('test_orig_map_64_qvar.png')
plt.clf()
hp.mollview(orig_map[5])
plt.savefig('test_orig_map_qvar.png')
plt.clf()
# hp.mollview(maps[5])
# plt.savefig('test_orig_map2_qvar.png')
# plt.clf()
hp.mollview(smoothed_map[4])
plt.savefig('test_smoothed_map_qvar.png')
plt.clf()
