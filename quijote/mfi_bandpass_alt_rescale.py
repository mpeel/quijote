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
import csv
def read_quijote_mfi_bandpass_alt(filename):
	freq = []
	bp1 = []
	bp2 = []
	bp3 = []
	bp4 = []
	print(filename)
	with open(filename, mode='r', encoding='cp1252') as f:
		for line in f:
			if '(' not in line and '!' not in line and '#' not in line:
				val = line.strip().split()
				freq.append(float(val[3])/1e9)
				if float(val[4]) > 0:
					bp1.append(float(val[4]))
				else:
					bp1.append(0.0)
				if float(val[5]) > 0:
					bp2.append(float(val[5]))
				else:
					bp2.append(0.0)
				if float(val[6]) > 0:
					bp3.append(float(val[6]))
				else:
					bp3.append(0.0)
				if float(val[7]) > 0:
					bp4.append(float(val[7]))
				else:
					bp4.append(0.0)
	return np.asarray([freq, bp1, bp2, bp3, bp4])

indirectory = '/Users/mpeel/Documents/QUIJOTE/MFI/2020-06_bandpass/'
plotdir = indirectory+'plots/'

reference = ''
comparison = ''

reference = read_quijote_mfi_bandpass_alt(indirectory+'2020-08-06/12-22-39_horn2_14_22GHz_ch13_16_mod_67_5deg_hornang_-45.dat')
comparison = read_quijote_mfi_bandpass_alt(indirectory+'2020-08-06/12-15-04_horn2_14_22GHz_ch13_16_mod_67_5deg_hornang_-45.dat')

vals = []
for i in range(0,4):
	ratio = np.asarray(reference[i+1]/comparison[i+1])
	mask = np.ones(len(ratio))
	mask[comparison[i+1] > 10.0] = 0
	mask[comparison[i+1] < np.median(comparison[i+1])+1] = 0
	# print(np.shape(ratio))
	print(np.median(ratio[mask==1]))
	vals.append(np.median(ratio[mask==1]))
	plt.plot(reference[i+1])
	plt.plot(comparison[i+1])
	plt.plot(reference[i+1]*mask)
	plt.savefig('bandpass_test_'+str(i+1)+'.png')
	plt.clf()

reference = read_quijote_mfi_bandpass_alt(indirectory+'2020-08-06/12-06-44_horn2_14_22GHz_ch13_16_mod_45deg_hornang_-45.dat')
comparison = read_quijote_mfi_bandpass_alt(indirectory+'2020-08-06/11-58-25_horn2_14_22GHz_ch13_16_mod_45deg_hornang_-45.dat')

for i in range(0,4):
	ratio = np.asarray(reference[i+1]/comparison[i+1])
	mask = np.ones(len(ratio))
	mask[comparison[i+1] > 10.0] = 0
	mask[comparison[i+1] < np.median(comparison[i+1])+0.4] = 0
	# print(np.shape(ratio))
	print(np.median(ratio[mask==1]))
	vals.append(np.median(ratio[mask==1]))
	plt.plot(reference[i+1])
	plt.plot(comparison[i+1])
	plt.plot(reference[i+1]*mask)
	plt.savefig('bandpass2_test_'+str(i+1)+'.png')
	plt.clf()

print(np.mean(vals))
print(np.median(vals))
# EOF