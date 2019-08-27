from astrocode.fitspectrum.smoothmap import smoothmap
import numpy as np
import healpy as hp
import astropy.io.fits as fits

output_resolution = 60.0
output_nside = [512]#, 256, 128, 64]

mfi_l,mfi_111,mfi_113,mfi_217,mfi_219,mfi_311,mfi_313,mfi_417,mfi_419 = np.loadtxt('/Users/mpeel/Documents/QUIJOTE/MFI/quijote_mfi_wfs.dat',unpack=True)


# July 2019 release
directory = '/Users/mpeel/Documents/maps/quijote_201907/reform/'
outdirectory = '/Users/mpeel/Documents/maps/quijote_201907/smooth/'
prefixes = ['mfi_jul2019','mfi_jul2019_half1','mfi_jul2019_half2']
newprefixes = ['mfi','half1mfi','half2mfi']
newdate = '201907'


# May 2019 release
# directory = '/Users/mpeel/Documents/maps/quijote_201905/reform/'
# outdirectory = '/Users/mpeel/Documents/maps/quijote_201905/smooth/'
# prefixes = ['mfi_may2019','mfi_may2019_half1','mfi_may2019_half2']
# newprefixes = ['mfi','half1mfi','half2mfi']
# newdate = '201905'

# October 2018 release
# directory = '/Users/mpeel/Documents/maps/quijote_201810/reform/'
# outdirectory = '/Users/mpeel/Documents/maps/quijote_201810/smooth/'
# prefixes = ['mfi_oct2018']
# newprefixes = ['mfi']
# newdate = '201810'

# April 2019 pre-release multi-epoch
# directory = '/Users/mpeel/Documents/maps/quijote_201904/reform/'
# outdirectory = '/Users/mpeel/Documents/maps/quijote_201904/smooth/'
# prefixes = ['mfi_apr2019_30_2','mfi_apr2019_35_6','mfi_apr2019_40_2','mfi_apr2019_40_5','mfi_apr2019_50_2','mfi_apr2019_50_5','mfi_apr2019_50_6','mfi_apr2019_60_1','mfi_apr2019_60_2','mfi_apr2019_60_5','mfi_apr2019_60_6','mfi_apr2019_65_1','mfi_apr2019_65_2','mfi_apr2019_65_6','mfi_apr2019_70_6']
# newprefixes = ['mfi302','mfi356','mfi402','mfi405','mfi502','mfi505','mfi506','mfi601','mfi602','mfi605','mfi606','mfi651','mfi652','mfi656','mfi706']
# newdate = '201904'

# March 2018 pre-release multi-epoch
# directory = '/Users/mpeel/Documents/maps/quijote_201803/reform/'
# outdirectory = '/Users/mpeel/Documents/maps/quijote_201803/smooth/'
# prefixes = ['mfi_mar2018_30_2','mfi_mar2018_35_6','mfi_mar2018_40_2','mfi_mar2018_40_5','mfi_mar2018_50_2','mfi_mar2018_50_5','mfi_mar2018_50_6','mfi_mar2018_60_1','mfi_mar2018_60_2','mfi_mar2018_60_5','mfi_mar2018_60_6','mfi_mar2018_65_1','mfi_mar2018_65_2','mfi_mar2018_65_6','mfi_mar2018_70_6']
# newprefixes = ['mfi302','mfi356','mfi402','mfi405','mfi502','mfi505','mfi506','mfi601','mfi602','mfi605','mfi606','mfi651','mfi652','mfi656','mfi706']
# newdate = '201803'

# Gaussian smoothing
for k in range(0,len(prefixes)):
	numnside = len(output_nside)
	for i in range(0,numnside):
		smoothmap(directory,outdirectory,prefixes[k]+'_map_11.0_1.fits',str(output_nside[i])+'_60.00smoothed_'+newprefixes[k]+'1_11.0_512_'+newdate+'_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=np.sqrt(mfi_111),usehealpixfits=True,taper_gauss=True,taper_gauss_sigma=0.886835*60,minmapvalue=-100,maxmapvalue=100,minmaxmaps=[1,2])
		smoothmap(directory,outdirectory,prefixes[k]+'_map_11.0_3.fits',str(output_nside[i])+'_60.00smoothed_'+newprefixes[k]+'3_11.0_512_'+newdate+'_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=np.sqrt(mfi_311),usehealpixfits=True,taper_gauss=True,taper_gauss_sigma=0.840227*60,minmapvalue=-100,maxmapvalue=100,minmaxmaps=[1,2])
		smoothmap(directory,outdirectory,prefixes[k]+'_map_13.0_1.fits',str(output_nside[i])+'_60.00smoothed_'+newprefixes[k]+'1_13.0_512_'+newdate+'_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=np.sqrt(mfi_113),usehealpixfits=True,taper_gauss=True,taper_gauss_sigma=0.892011*60,minmapvalue=-100,maxmapvalue=100,minmaxmaps=[1,2])
		smoothmap(directory,outdirectory,prefixes[k]+'_map_13.0_3.fits',str(output_nside[i])+'_60.00smoothed_'+newprefixes[k]+'3_13.0_512_'+newdate+'_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=np.sqrt(mfi_313),usehealpixfits=True,taper_gauss=True,taper_gauss_sigma=0.845983*60,minmapvalue=-100,maxmapvalue=100,minmaxmaps=[1,2])
		smoothmap(directory,outdirectory,prefixes[k]+'_map_17.0_2.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_'+newprefixes[k]+'2_17.0_512_'+newdate+'_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=np.sqrt(mfi_217),usehealpixfits=True,taper_gauss=True,taper_gauss_sigma=0.630364*60,minmapvalue=-100,maxmapvalue=100,minmaxmaps=[1,2])
		smoothmap(directory,outdirectory,prefixes[k]+'_map_17.0_4.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_'+newprefixes[k]+'4_17.0_512_'+newdate+'_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=np.sqrt(mfi_417),usehealpixfits=True,taper_gauss=True,taper_gauss_sigma=0.651673*60,minmapvalue=-100,maxmapvalue=100,minmaxmaps=[1,2])
		smoothmap(directory,outdirectory,prefixes[k]+'_map_19.0_2.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_'+newprefixes[k]+'2_19.0_512_'+newdate+'_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=np.sqrt(mfi_219),usehealpixfits=True,taper_gauss=True,taper_gauss_sigma=0.630364*60,minmapvalue=-100,maxmapvalue=100,minmaxmaps=[1,2])
		smoothmap(directory,outdirectory,prefixes[k]+'_map_19.0_4.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_'+newprefixes[k]+'4_19.0_512_'+newdate+'_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=np.sqrt(mfi_419),usehealpixfits=True,taper_gauss=True,taper_gauss_sigma=0.651673*60,minmapvalue=-100,maxmapvalue=100,minmaxmaps=[1,2])


# This was cosine smoothing
# for k in range(0,len(prefixes)):
# 	numnside = len(output_nside)
# 	for i in range(0,numnside):
# 		smoothmap(directory,outdirectory,prefixes[k]+'_map_11.0_1.fits',str(output_nside[i])+'_60.00smoothed_'+newprefixes[k]+'1_11.0_512_'+newdate+'_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=np.sqrt(mfi_111),usehealpixfits=True,taper=True,lmin_taper=300,lmax_taper=450,minmapvalue=-100,maxmapvalue=100,minmaxmaps=[1,2])
# 		smoothmap(directory,outdirectory,prefixes[k]+'_map_11.0_3.fits',str(output_nside[i])+'_60.00smoothed_'+newprefixes[k]+'3_11.0_512_'+newdate+'_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=np.sqrt(mfi_311),usehealpixfits=True,taper=True,lmin_taper=300,lmax_taper=450,minmapvalue=-100,maxmapvalue=100,minmaxmaps=[1,2])
# 		smoothmap(directory,outdirectory,prefixes[k]+'_map_13.0_1.fits',str(output_nside[i])+'_60.00smoothed_'+newprefixes[k]+'1_13.0_512_'+newdate+'_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=np.sqrt(mfi_113),usehealpixfits=True,taper=True,lmin_taper=300,lmax_taper=450,minmapvalue=-100,maxmapvalue=100,minmaxmaps=[1,2])
# 		smoothmap(directory,outdirectory,prefixes[k]+'_map_13.0_3.fits',str(output_nside[i])+'_60.00smoothed_'+newprefixes[k]+'3_13.0_512_'+newdate+'_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=np.sqrt(mfi_313),usehealpixfits=True,taper=True,lmin_taper=300,lmax_taper=450,minmapvalue=-100,maxmapvalue=100,minmaxmaps=[1,2])
# 		smoothmap(directory,outdirectory,prefixes[k]+'_map_17.0_2.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_'+newprefixes[k]+'2_17.0_512_'+newdate+'_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=np.sqrt(mfi_217),usehealpixfits=True,taper=True,lmin_taper=350,lmax_taper=600,minmapvalue=-100,maxmapvalue=100,minmaxmaps=[1,2])
# 		smoothmap(directory,outdirectory,prefixes[k]+'_map_17.0_4.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_'+newprefixes[k]+'4_17.0_512_'+newdate+'_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=np.sqrt(mfi_417),usehealpixfits=True,taper=True,lmin_taper=350,lmax_taper=600,minmapvalue=-100,maxmapvalue=100,minmaxmaps=[1,2])
# 		smoothmap(directory,outdirectory,prefixes[k]+'_map_19.0_2.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_'+newprefixes[k]+'2_19.0_512_'+newdate+'_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=np.sqrt(mfi_219),usehealpixfits=True,taper=True,lmin_taper=350,lmax_taper=600,minmapvalue=-100,maxmapvalue=100,minmaxmaps=[1,2])
# 		smoothmap(directory,outdirectory,prefixes[k]+'_map_19.0_4.fits',str(output_nside[i])+'_'+str(output_resolution)+'smoothed_'+newprefixes[k]+'4_19.0_512_'+newdate+'_mKCMBunits.fits', output_resolution,nside_out=output_nside[i],windowfunction=np.sqrt(mfi_419),usehealpixfits=True,taper=True,lmin_taper=350,lmax_taper=600,minmapvalue=-100,maxmapvalue=100,minmaxmaps=[1,2])

# EOF