#!/usr/bin/env python
# -*- coding: utf-8  -*-
#
# Do a comparison of calculated MFI bandpasses
# 
# Version history:
#
# 02-Jun-2020  M. Peel       Started
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits

def read_quijote_mfi_bandpass(filename):
	freq = []
	bp1 = []
	bp2 = []
	bp3 = []
	bp4 = []
	bp5 = []
	bp6 = []
	bp7 = []
	bp8 = []
	with open(filename) as f:
		for line in f:
			if '(' not in line:
				val = line.strip().split()
				freq.append(float(val[0]))
				if float(val[1]) > 0:
					bp1.append(float(val[1]))
				else:
					bp1.append(0.0)
				if float(val[2]) > 0:
					bp2.append(float(val[2]))
				else:
					bp2.append(0.0)
				if float(val[3]) > 0:
					bp3.append(float(val[3]))
				else:
					bp3.append(0.0)
				if float(val[4]) > 0:
					bp4.append(float(val[4]))
				else:
					bp4.append(0.0)
				if float(val[5]) > 0:
					bp5.append(float(val[5]))
				else:
					bp5.append(0.0)
				if float(val[6]) > 0:
					bp6.append(float(val[6]))
				else:
					bp6.append(0.0)
				if float(val[7]) > 0:
					bp7.append(float(val[7]))
				else:
					bp7.append(0.0)
				if float(val[8]) > 0:
					bp8.append(float(val[8]))
				else:
					bp8.append(0.0)
	bp1 = bp1 / np.sum(bp1)
	bp2 = bp2 / np.sum(bp2)
	bp3 = bp3 / np.sum(bp3)
	bp4 = bp4 / np.sum(bp4)
	bp5 = bp5 / np.sum(bp5)
	bp6 = bp6 / np.sum(bp6)
	bp7 = bp7 / np.sum(bp7)
	bp8 = bp8 / np.sum(bp8)
	return np.asarray([freq, bp1, bp2, bp3, bp4, bp5, bp6, bp7, bp8])


indirectory = '/Users/mpeel/Documents/maps/quijote_mfi/'
plotdir = indirectory+'plots/'
files = ['Pol3_bandpass.dat','2020_05_21_Pol3_bandpass.dat','2020_05_22_0deg_Pol3_bandpass.dat','2020_05_22_22p5deg_Pol3_bandpass.dat']
replace = ['Pol3','Pol2','Pol4']
# indirectory = '/Users/mpeel/Documents/QUIJOTE/MFI/2020-05_MFI_measurements/_plots/'
plotdir = indirectory+'plots2/'
files = ['Pol3_bandpass.dat','Pol3_bandpass_mod_0deg_hornang_0.dat']
files = ['Pol3_bandpass_sum_mod_0deg_hornang_0.dat','Pol3_bandpass_sum_mod_0deg_hornang_90.dat','Pol3_bandpass_sum_mod_22_5deg_hornang_0.dat','Pol3_bandpass_sum_mod_22_5deg_hornang_90.dat','Pol3_bandpass_sum_mod_45deg_hornang_0.dat','Pol3_bandpass_sum_mod_45deg_hornang_90.dat','Pol3_bandpass_sum_mod_67_5deg_hornang_0.dat','Pol3_bandpass_sum_mod_67_5deg_hornang_90.dat','Pol3_bandpass.dat']
replace = ['Pol3']


numdatasets = len(files)
numchannels = 8

for k in range(0,len(replace)):
	datasets = []
	for file in files:
		newset = read_quijote_mfi_bandpass(indirectory+file.replace(replace[0],replace[k]))
		datasets.append(newset)

	for i in range(1,numchannels+1):
		for j in range(0,numdatasets):
			label = files[j].split('band')[1].replace('pass','').replace('dat','')
			label = label.replace('_sum_mod','mod')
			print(label)
			if label == '.':
				plt.plot(datasets[j][0],datasets[j][i]/(datasets[j][0][1]-datasets[j][0][0]),'b-',label=label,lw=0.5)
			else:
				plt.plot(datasets[j][0],datasets[j][i]/(datasets[j][0][1]-datasets[j][0][0]),'-',label=label,lw=0.5)

		plt.legend()
		plt.ylabel('Amplitude (renormalised)')
		plt.xlabel('Frequency (GHz)')
		plt.savefig(plotdir+replace[k]+'_'+str(i)+'.png')
		plt.clf()
