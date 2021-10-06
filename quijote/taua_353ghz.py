import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
from astrocode.fitspectrum.astroutils import *


def calc_polang(Q, U):
	return 0.5*np.arctan2(Q, U) * 180 / np.pi

def calc_polang2(Q, U):
	return 0.5*np.arctan2(-U, Q) * 180 / np.pi

def calc_polang_unc(Q, U, Qerr, Uerr):
	unc_map = np.sqrt((Qerr**2)*(-0.5*U/(Q**2.0+U**2.0))**2.0 + (Uerr**2)*(-0.5*Q/(Q**2.0+U**2.0))**2.0) * 180 / np.pi
	unc_map[~np.isfinite(unc_map)] = 1000.0
	return unc_map

# print(calc_polang(10.0,0.0))
# print(calc_polang2(10.0,0.0))
# exit()

frequencies = ['100','143','217','353']
# frequencies = ['353']
for freq in frequencies:
	maps = hp.read_map('/Users/mpeel/Documents/maps/planck2020/HFI_SkyMap_'+freq+'_2048_R4.00_full.fits',field=(1,2))
	hp.mollview(maps[0],norm='hist')
	plt.savefig('p'+freq+'_Q.png')
	plt.clf()
	# maps[1] = -maps[1]
	hp.mollview(maps[1],norm='hist')
	plt.savefig('p'+freq+'_U.png')
	plt.clf()
	p353pol = np.sqrt(maps[0]*maps[0]+maps[1]*maps[1])
	hp.mollview(p353pol,norm='hist')
	plt.savefig('p'+freq+'pol.png')
	plt.clf()

	lon = 184.55745
	lat = -05.78436
	radii = np.arange(1.0,60.0,1.0)
	aper_outer_radius1 = 60.0
	aper_outer_radius2 = 80.0
	fd_q = np.zeros(len(radii))
	fd_u = np.zeros(len(radii))
	fd_p = np.zeros(len(radii))
	fd_q_err = np.zeros(len(radii))
	fd_u_err = np.zeros(len(radii))
	fd_p_err = np.zeros(len(radii))
	units='KCMB'
	aper_inner_radius = 10.0
	fd_q_5, fd_q_err_5, fd_bg = haperflux(maps[0], float(freq), 5.0, lon, lat, aper_inner_radius, aper_outer_radius1, aper_outer_radius2, units,quickplot='p'+freq+'_taua_Q.png')
	fd_u_5, fd_u_err_5, fd_bg = haperflux(maps[1], float(freq), 5.0, lon, lat, aper_inner_radius, aper_outer_radius1, aper_outer_radius2, units,quickplot='p'+freq+'_taua_U.png')
	haperflux(p353pol, float(freq), 5.0, lon, lat, aper_inner_radius, aper_outer_radius1, aper_outer_radius2, units,quickplot='p'+freq+'_taua.png')
	for i in range(0,len(radii)):
		aper_inner_radius = radii[i]
		fd_q[i], fd_q_err[i], fd_bg = haperflux(maps[0], float(freq), 5.0, lon, lat, aper_inner_radius, aper_outer_radius1, aper_outer_radius2, units,quickplot='')
		fd_u[i], fd_u_err[i], fd_bg = haperflux(maps[1], float(freq), 5.0, lon, lat, aper_inner_radius, aper_outer_radius1, aper_outer_radius2, units,quickplot='')
		fd_p[i], fd_p_err[i], fd_bg = haperflux(p353pol, float(freq), 5.0, lon, lat, aper_inner_radius, aper_outer_radius1, aper_outer_radius2, units,quickplot='')

	plt.errorbar(radii, fd_q, yerr=fd_q_err, label='Q')
	plt.errorbar(radii, fd_u, yerr=fd_u_err, label='U')
	plt.errorbar(radii, fd_p, yerr=fd_p_err, label='P')
	plt.xlabel('Distance (arcmin)')
	plt.ylabel('Flux density (Jy)')
	plt.legend()
	plt.savefig('p'+freq+'_tau_radii.png')
	plt.clf()
	# angle = calc_polang(fd_q, fd_u)
	angle = calc_polang2(fd_q, fd_u)
	angle[angle>0] = -180.0+angle[angle>0]
	angle_err = calc_polang_unc(fd_q, fd_u,fd_q_err,fd_u_err)
	plt.errorbar(radii, angle, yerr=angle_err)
	plt.xlabel('Distance (arcmin)')
	plt.ylabel('Polarisation angle')
	plt.savefig('p'+freq+'_tau_angle_radii.png')

	angle = calc_polang2(fd_q-fd_q_5, fd_u-fd_u_5)
	angle[angle>0] = -180.0+angle[angle>0]
	angle[angle>90] = np.nan
	angle[angle<-180] = np.nan
	angle_err = calc_polang_unc(fd_q-fd_q_5, fd_u-fd_u_5,fd_q_err,fd_u_err)
	angle_err[angle_err>90] = 0.0
	angle_err[angle_err<-180] = 0.0
	plt.errorbar(radii, angle, yerr=angle_err)
	plt.xlim(6,60)
	plt.xlabel('Distance (arcmin)')
	plt.ylabel('Polarisation angle')
	plt.savefig('p'+freq+'_tau_angle_radii_notaua.png')