#!/usr/bin/env python
# -*- coding: utf-8  -*-
#
# Have a quick look at the flux density of Perseus over time
# 
# Version history:
#
# 30-Apr-2020  M. Peel       Started

import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
from astrocode.fitspectrum.astroutils import *

res_arcmin = 60.0
lon = 159.73
lat = -18.60
aper_inner_radius = 100.0
aper_outer_radius1 = 100.0
aper_outer_radius2 = 140.0
noise_model = 1
units='mK_CMB'

basedir = '/Users/mpeel/Documents/maps/quijote_202103/reform/'
maps = ['mfi_mar2021_period1_13.0_3.fits','mfi_mar2021_period2_13.0_3.fits','mfi_mar2021_period5_13.0_3.fits','mfi_mar2021_period6_13.0_3.fits']
# for filename in filenames:
fd = []
fd_err = []
for filename in maps:
	mapdata = hp.read_map(basedir + filename)
	vals = haperflux(mapdata, 13.0, res_arcmin, lon, lat, aper_inner_radius, aper_outer_radius1, aper_outer_radius2, units, column=0, dopol=False, nested=False, noise_model=noise_model, abratio=1.0, angle=0.0, silent=False)
	fd.append(vals[0])
	fd_err.append(vals[1])
plt.errorbar(range(0,len(fd)),fd,yerr=fd_err)
plt.savefig('test_perseus_periods.png')