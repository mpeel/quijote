#!/usr/bin/env python
# -*- coding: utf-8  -*-
#
# Reduce the MFI polarisation tests from May 2020
# 
# Version history:
#
# 13-May-2020  M. Peel       Started

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits

indirectory = '/Users/mpeel/Documents/QUIJOTE/MFI/2020-05_MFI_measurements/'
files = ['horn3_poltest-200512-1356-000.tod','horn3_linearity-200512-1420-000.tod','horn3_linearity-200512-1426-000.tod','HORN3_lintest-200515-1237-000.tod']

files = ['horn2_lintest-200515-1257-000.tod','HORN2_poltest-200515-1248-000.tod','HORN2_poltest-200515-1250-000.tod','HORN2_poltest-200515-1252-000.tod','horn3_linearity-200512-1420-000.tod','horn3_linearity-200512-1426-000.tod','HORN3_lintest-200515-1237-000.tod','HORN3_lintest-200515-1238-000.tod','horn3_poltest-200512-1356-000.tod','horn3_poltest-200512-1402-000.tod','horn3_poltest-200512-135921-000.tod','horn3_poltest-200512-135950-000.tod','HORN3_poltest-200515-1227-000.tod','HORN3_poltest-200515-1228-000.tod','HORN3_poltest-200515-1230-000.tod','HORN3_poltest-200515-1232-000.tod','HORN3_poltest-200515-1234-000.tod','HORN4_lintes-200515-1139-000.tod','HORN4_poltest-200515-1124-000.tod','HORN4_poltest-200515-1128-000.tod','HORN4_poltest-200515-1131-000.tod','HORN4_poltest-200515-1133-000.tod']
horns = [2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4]
# horns = [3,3,3,3]
numchannels = 8

for i in range(0,len(horns)):
	start = numchannels*(horns[i]-1)
	end = numchannels*horns[i]
	inputfits = fits.open(indirectory+files[i])
	# print(inputfits[1].header)
	data = inputfits[1].data.field('data')[0][:,start:end]
	print(np.shape(data))
	data_diff = np.diff(data,axis=0)
	print(np.shape(data_diff))
	# print(data)
	# print(np.shape(data[0]))
	# print(np.shape(data[:,1]))
	# print(numchannels*(horns[i]-1))
	# print(numchannels*horns[i])

	tempmask = np.ones(len(data_diff[:,0]))
	print(len(tempmask))
	print(np.median(data_diff[:,0]))
	print(np.std(data_diff[:,0]))
	numsigma = 5
	tempmask[data_diff[:,0] > np.median(data_diff[:,0]) + numsigma*np.std(data_diff[:,0])] = 0
	tempmask[data_diff[:,0] < np.median(data_diff[:,0]) - numsigma*np.std(data_diff[:,0])] = 0
	smooth = 150
	for j in range(0,int(len(tempmask)/smooth)):
		if np.sum(tempmask[j*smooth:j*smooth+smooth]) != smooth:
			tempmask[j*smooth:j*smooth+smooth] = 0

	mask = tempmask.copy()

	# Find out whether we're going up or down in elevation
	state = -1
	trip = 1
	count = 0
	for j in range(0,len(mask)):
		if tempmask[j] == 0:
			# print('hi')
			if trip == 0:
				# if state < 0:
				state = count
				# else:
				# 	state = -count
				trip = 1
				count += 1
		else:
			mask[j] = state
			trip = 0
	mask = np.pad(mask,(0,1),'constant')
	print(len(mask))
	print(len(data[:,0]))
	data[mask==0,:] = float('NaN')
	for j in range(0,len(data[0])):
	# for j in range(0,1):
		print(j)
		plt.plot(data[:,j])
		# plt.plot(data_diff[:,j])
		plt.plot(mask/100)
		# plt.plot(tempmask)
	plt.savefig(indirectory+files[i]+'.png')
	plt.clf()
	if count != 0:
		print(count)
		vals = np.zeros((count-1,len(data[0])))
		print(np.shape(vals))
		for k in range(0,len(data[0])):
			for j in range(1,count):
				vals[j-1,k] = np.nanmedian(data[mask==j,k])
				print(str(j) + " - " + str(np.sum(mask==j)) + " - " + str(vals[j-1,k]))
		diff = np.diff(vals,axis=0)
		plt.plot(diff[::2])
		plt.plot(diff[1::2],'.-')
		plt.xlabel('Switch #')
		plt.ylabel('Diode amplitude [V]')
		# plt.legend()
		plt.savefig(indirectory+files[i]+'_diode.png')
