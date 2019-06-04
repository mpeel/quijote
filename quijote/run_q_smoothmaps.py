from astrocode.fitspectrum.smoothmap import smoothmap
import numpy as np
import healpy as hp
import astropy.io.fits as fits

output_resolution = 60.0
output_nside = [512]#, 256, 128, 64]

directory = '/Users/mpeel/Documents/maps/quijote_201905/reform/'
outdirectory = '/Users/mpeel/Documents/maps/quijote_201905/smooth/'

mfi_l,mfi_111,mfi_113,mfi_217,mfi_219,mfi_311,mfi_313,mfi_417,mfi_419 = np.loadtxt('/Users/mpeel/Documents/QUIJOTE/MFI/quijote_mfi_wfs.dat',unpack=True)

numnside = len(output_nside)
for i in range(0,numnside):
	smoothmap(directory,outdirectory,'mfi_map_11.0_1.fits',str(output_nside[i])+'_60.00smoothed_mfi1_11.0_512_201905_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=np.sqrt(mfi_111),usehealpixfits=True,taper=True,lmin_taper=300,lmax_taper=450)
	smoothmap(directory,outdirectory,'mfi_map_11.0_3.fits',str(output_nside[i])+'_60.00smoothed_mfi3_11.0_512_201905_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=np.sqrt(mfi_311),usehealpixfits=True,taper=True,lmin_taper=300,lmax_taper=450)
	smoothmap(directory,outdirectory,'mfi_map_13.0_1.fits',str(output_nside[i])+'_60.00smoothed_mfi1_13.0_512_201905_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=np.sqrt(mfi_113),usehealpixfits=True,taper=True,lmin_taper=300,lmax_taper=450)
	smoothmap(directory,outdirectory,'mfi_map_13.0_3.fits',str(output_nside[i])+'_60.00smoothed_mfi3_13.0_512_201905_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=np.sqrt(mfi_313),usehealpixfits=True,taper=True,lmin_taper=300,lmax_taper=450)
	smoothmap(directory,outdirectory,'mfi_map_17.0_2.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_mfi2_17.0_512_201905_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=np.sqrt(mfi_217),usehealpixfits=True,taper=True,lmin_taper=350,lmax_taper=600)
	smoothmap(directory,outdirectory,'mfi_map_17.0_4.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_mfi4_17.0_512_201905_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=np.sqrt(mfi_417),usehealpixfits=True,taper=True,lmin_taper=350,lmax_taper=600)
	smoothmap(directory,outdirectory,'mfi_map_19.0_2.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_mfi2_19.0_512_201905_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=np.sqrt(mfi_219),usehealpixfits=True,taper=True,lmin_taper=350,lmax_taper=600)
	smoothmap(directory,outdirectory,'mfi_map_19.0_4.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_mfi4_19.0_512_201905_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=np.sqrt(mfi_419),usehealpixfits=True,taper=True,lmin_taper=350,lmax_taper=600)

# EOF