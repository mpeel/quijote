from astrocode.fitspectrum.smoothnoisemap import smoothnoisemap
import numpy as np
import healpy as hp
import astropy.io.fits as fits

output_resolution = [60.0]
output_nside = [512]#, 256, 128, 64]

# directory = '/Users/mpeel/Documents/maps/quijote_201905/reform/'
# outdirectory = '/Users/mpeel/Documents/maps/quijote_201905/smooth/'


#mfi_l,mfi_111,mfi_113,mfi_217,mfi_219,mfi_311,mfi_313,mfi_417,mfi_419 = np.loadtxt('/Users/mpeel/Documents/QUIJOTE/MFI/quijote_mfi_wfs.dat',unpack=True)

mfi_l,mfi_111,mfi_113,mfi_217,mfi_219,mfi_311,mfi_313,mfi_417,mfi_419 = np.loadtxt('/net/nas/proyectos/quijote/mikepeel/quijote_mfi_wfs.dat',unpack=True)


# July 2019 release
directory = '/net/nas/proyectos/quijote/mikepeel/quijote_201907/reform/'
outdirectory = '/net/nas/proyectos/quijote/mikepeel/quijote_201907/smooth/'
prefixes = ['mfi_jul2019','mfi_jul2019_half1','mfi_jul2019_half2']
newprefixes = ['mfi','half1mfi','half2mfi']
newdate = '201907'

# May 2019 release
# directory = '/net/nas/proyectos/quijote/mikepeel/quijote_201905/reform/'
# outdirectory = '/net/nas/proyectos/quijote/mikepeel/quijote_201905/smooth/'
# prefixes = ['mfi_may2019','mfi_may2019_half1','mfi_may2019_half2']
# newprefixes = ['mfi','half1mfi','half2mfi']
# newdate = '201905'
numrealisations=1000

mapnumbers = [0,1,2]

for k in range(0,len(prefixes)):
	numres = len(output_resolution)
	for i in range(0,numres):
		resolution = "%.2f" % output_resolution[i]
		for m in range(0,len(mapnumbers)):
			# smoothnoisemap(directory, outdirectory, str(output_nside[i])+'_'+resolution+'smoothed_'+newprefixes[k]+'2_17.0_512_'+newdate+'_nhits_'+str(mapnumbers[m]), prefixes[k]+'_nhits_17.0_2.fits',mapnumber=mapnumbers[m],numrealisations=numrealisations,fwhm=output_resolution[i],sigma_0=1.0*np.sqrt(25),nside=output_nside,windowfunction=np.sqrt(mfi_217),usehealpixfits=True)

			# smoothnoisemap(directory, outdirectory, str(output_nside[i])+'_'+resolution+'smoothed_'+newprefixes[k]+'2_17.0_512_'+newdate+'_nhits_'+str(mapnumbers[m]), prefixes[k]+'_nhits_17.0_2.fits',mapnumber=mapnumbers[m],numrealisations=numrealisations,fwhm=output_resolution[i],sigma_0=1.0*np.sqrt(25),nside=output_nside,windowfunction=np.sqrt(mfi_217),usehealpixfits=True,taper_gauss=True)

			smoothnoisemap(directory, outdirectory, str(output_nside[i])+'_'+resolution+'smoothed_'+newprefixes[k]+'2_17.0_512_'+newdate+'_weight_'+str(mapnumbers[m]), prefixes[k]+'_weights_17.0_2.fits',mapnumber=mapnumbers[m],numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=np.sqrt(mfi_217),usehealpixfits=True,taper_gauss=True,taper_gauss_sigma=0.630364*60,normalise=False)

			smoothnoisemap(directory, outdirectory, str(output_nside[i])+'_'+resolution+'smoothed_'+newprefixes[k]+'2_19.0_512_'+newdate+'_weight_'+str(mapnumbers[m]), prefixes[k]+'_weights_19.0_2.fits',mapnumber=mapnumbers[m],numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=np.sqrt(mfi_219),usehealpixfits=True,taper_gauss=True,taper_gauss_sigma=0.630364*60,normalise=False)

			smoothnoisemap(directory, outdirectory, str(output_nside[i])+'_'+resolution+'smoothed_'+newprefixes[k]+'3_11.0_512_'+newdate+'_weight_'+str(mapnumbers[m]), prefixes[k]+'_weights_11.0_3.fits',mapnumber=mapnumbers[m],numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=np.sqrt(mfi_311),usehealpixfits=True,taper_gauss=True,taper_gauss_sigma=0.840227*60,normalise=False)

			smoothnoisemap(directory, outdirectory, str(output_nside[i])+'_'+resolution+'smoothed_'+newprefixes[k]+'3_13.0_512_'+newdate+'_weight_'+str(mapnumbers[m]), prefixes[k]+'_weights_13.0_3.fits',mapnumber=mapnumbers[m],numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=np.sqrt(mfi_313),usehealpixfits=True,taper_gauss=True,taper_gauss_sigma=0.845983*60,normalise=False)

			smoothnoisemap(directory, outdirectory, str(output_nside[i])+'_'+resolution+'smoothed_'+newprefixes[k]+'4_17.0_512_'+newdate+'_weight_'+str(mapnumbers[m]), prefixes[k]+'_weights_17.0_4.fits',mapnumber=mapnumbers[m],numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=np.sqrt(mfi_417),usehealpixfits=True,taper_gauss=True,taper_gauss_sigma=0.651673*60,normalise=False)

			smoothnoisemap(directory, outdirectory, str(output_nside[i])+'_'+resolution+'smoothed_'+newprefixes[k]+'4_19.0_512_'+newdate+'_weight_'+str(mapnumbers[m]), prefixes[k]+'_weights_19.0_4.fits',mapnumber=mapnumbers[m],numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=np.sqrt(mfi_417),usehealpixfits=True,taper_gauss=True,taper_gauss_sigma=0.651673*60,normalise=False)

# EOF