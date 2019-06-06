from astrocode.fitspectrum.smoothnoisemap import smoothnoisemap
import numpy as np
import healpy as hp
import astropy.io.fits as fits

output_resolution = [60.0]
output_nside = [512]#, 256, 128, 64]

directory = '/Users/mpeel/Documents/maps/quijote_201905/reform/'
outdirectory = '/Users/mpeel/Documents/maps/quijote_201905/smooth/'

mfi_l,mfi_111,mfi_113,mfi_217,mfi_219,mfi_311,mfi_313,mfi_417,mfi_419 = np.loadtxt('/Users/mpeel/Documents/QUIJOTE/MFI/quijote_mfi_wfs.dat',unpack=True)
prefixes = ['mfi_may2019']#,'mfi_may2019_half1','mfi_may2019_half2']
newprefixes = ['mfi']#,'half1mfi','half2mfi']
newdate = '201905'
numrealisations=100#0

mapnumbers = [0,1,2]

for k in range(0,len(prefixes)):
	numres = len(output_resolution)
	for i in range(0,numres):
		resolution = "%.2f" % output_resolution[i]
		for m in range(0,len(mapnumbers)):
			smoothnoisemap(directory, outdirectory, str(output_nside[i])+'_'+resolution+'smoothed_'+newprefixes[k]+'2_17.0_512_'+newdate+'_nhits_'+str(mapnumbers[m]), prefixes[k]+'_nhits_17.0_2.fits',mapnumber=mapnumbers[m],numrealisations=numrealisations,fwhm=output_resolution[i],sigma_0=1.0*np.sqrt(25),nside=output_nside,windowfunction=np.sqrt(mfi_217),usehealpixfits=True)

			smoothnoisemap(directory, outdirectory, str(output_nside[i])+'_'+resolution+'smoothed_'+newprefixes[k]+'2_17.0_512_'+newdate+'_weight_'+str(mapnumbers[m]), prefixes[k]+'_weights_17.0_2.fits',mapnumber=mapnumbers[m],numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=np.sqrt(mfi_217),usehealpixfits=True)
