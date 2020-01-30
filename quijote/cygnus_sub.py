# Quickly read in the QUIJOTE maps, focus on the Cygnus region, and try subtracting the intensity map from the polarisation to see what amplitude is needed.
# Mike Peel		30-Jan-2020		Started
import numpy as np
import healpy as hp
import scipy as sp
import matplotlib.pyplot as plt

def difference(P, x,y):
	return np.linalg.norm(y-P*x)

basedir="/Users/mpeel/Documents/maps/quijote_201911/"
outdir = basedir+'analyse/'
mapdir = basedir+'reform/'

inputmapnames = ['../../wmap9/wmap_band_smth_iqumap_r9_9yr_Ka_v5.fits']#'mfi_nov2019_mapsmth_11.0_3.fits', 'mfi_nov2019_mapsmth_13.0_3.fits', 'mfi_nov2019_mapsmth_17.0_4.fits', 'mfi_nov2019_mapsmth_19.0_4.fits', 'mfi_nov2019_mapsmth_17.0_2.fits', 'mfi_nov2019_mapsmth_19.0_2.fits']
region_long = [70, 90]
region_lat = [-4, 4]

for inputmapname in inputmapnames:
	inputmap = hp.read_map(mapdir+inputmapname,field=None)
	# hp.mollview(inputmap)
	# plt.savefig(outdir+inputmapname.replace('.fits','')+'_cygnus.pdf')
	# plt.clf()

	print(np.shape(inputmap))
	numpix = len(inputmap[0])
	nside = hp.pixelfunc.npix2nside(numpix)
	print(numpix)
	print(nside)
	positions = hp.pixelfunc.pix2ang(nside,range(0,numpix))
	phi = positions[1]*(180.0/np.pi)
	theta = 90.0-(positions[0]*180.0/np.pi)
	for i in range(0,3):
		inputmap[i][phi > region_long[1]] = hp.UNSEEN
		inputmap[i][phi < region_long[0]] = hp.UNSEEN
		inputmap[i][theta > region_lat[1]] = hp.UNSEEN
		inputmap[i][theta < region_lat[0]] = hp.UNSEEN
		hp.mollview(inputmap[i])
		plt.savefig(outdir+inputmapname.replace('.fits','')+'_cygnus_'+str(i)+'.pdf')
		plt.clf()
		hp.gnomview(inputmap[i],rot=[np.mean(region_long),np.mean(region_lat)],reso=1.3,xsize=1000,ysize=500)
		plt.savefig(outdir+inputmapname.replace('.fits','')+'_cygnus_'+str(i)+'_zoom.pdf')
		plt.clf()


	a = 0.0
	b = 0.0
	result = sp.optimize.minimize_scalar(difference, args=(inputmap[0][inputmap[0][:] != hp.UNSEEN], inputmap[1][inputmap[1][:] != hp.UNSEEN]))
	print(result)
	submap = inputmap[1]
	submap[inputmap[0] != hp.UNSEEN] = submap[inputmap[0] != hp.UNSEEN] -result.x*inputmap[0][inputmap[0] != hp.UNSEEN]
	hp.gnomview(submap,rot=[np.mean(region_long),np.mean(region_lat)],reso=1.3,xsize=1000,ysize=500,title=str(result.x))
	plt.savefig(outdir+inputmapname.replace('.fits','')+'_cygnus_Q_zoom.pdf')
	plt.clf()

	result = sp.optimize.minimize_scalar(difference, args=(inputmap[0][inputmap[0][:] != hp.UNSEEN], inputmap[2][inputmap[1][:] != hp.UNSEEN]))
	print(result)

	submap = inputmap[2]
	submap[inputmap[0] != hp.UNSEEN] = submap[inputmap[0] != hp.UNSEEN] -result.x*inputmap[0][inputmap[0] != hp.UNSEEN]
	hp.gnomview(submap,rot=[np.mean(region_long),np.mean(region_lat)],reso=1.3,xsize=1000,ysize=500,title=str(result.x))
	plt.savefig(outdir+inputmapname.replace('.fits','')+'_cygnus_U_zoom.pdf')
	plt.clf()
