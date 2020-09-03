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
	# with open(usequarry, mode='r') as infile:
		# reader = csv.reader(infile)
	# data = np.loadtxt(filename,comments=['#','!'])
	# print(data)
	# exit()
	with open(filename, mode='r', encoding='cp1252') as f:
		# reader = csv.reader(f)
		# lines = {rows[0] for rows in reader}
		for line in f:
			# print(line)
			if '(' not in line and '!' not in line and '#' not in line:
				val = line.strip().split()
				# print(val)
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
	# print(freq)
	# bp1 = bp1 - np.median(bp1)
	# bp2 = bp2 - np.median(bp2)
	# bp3 = bp3 - np.median(bp3)
	# bp4 = bp4 - np.median(bp4)
	# bp1 = bp1 / np.sum(bp1)
	# bp2 = bp2 / np.sum(bp2)
	# bp3 = bp3 / np.sum(bp3)
	# bp4 = bp4 / np.sum(bp4)
	return np.asarray([freq, bp1, bp2, bp3, bp4])

horns = [2,3,4]
indirectory = '/Users/mpeel/Documents/QUIJOTE/MFI/2020-06_bandpass/'
plotdir = indirectory+'plots/'

for dohorn in horns:
	if dohorn == 3:
		files_1 = ['2020-06-24/11-14-48_horn3_9_15GHz_ch17_20_mod_0deg_hornang_0.dat','2020-06-24/11-48-17_horn3_9_15GHz_ch17_20_mod_22_5deg_hornang_0.dat','2020-06-24/11-58-05_horn3_9_15GHz_ch17_20_mod_45deg_hornang_0.dat','2020-06-24/12-26-26_horn3_9_15GHz_ch17_20_mod_67_5deg_hornang_0.dat','2020-06-24/14-00-39_horn3_9_15GHz_ch17_20_mod_0deg_hornang_45.dat','2020-06-24/13-36-56_horn3_9_15GHz_ch17_20_mod_22_5deg_hornang_45.dat','2020-06-24/13-29-03_horn3_9_15GHz_ch17_20_mod_45deg_hornang_45.dat','2020-06-24/13-05-51_horn3_9_15GHz_ch17_20_mod_67_5deg_hornang_45.dat','2020-06-24/14-08-31_horn3_9_15GHz_ch17_20_mod_0deg_hornang_90.dat','2020-06-24/14-34-09_horn3_9_15GHz_ch17_20_mod_22_5deg_hornang_90.dat','2020-06-24/14-42-18_horn3_9_15GHz_ch17_20_mod_45deg_hornang_90.dat','2020-06-24/15-15-07_horn3_9_15GHz_ch17_20_mod_67_5deg_hornang_90.dat','2020-06-24/16-17-49_horn3_9_15GHz_ch17_20_mod_0deg_hornang135.dat','2020-06-24/15-54-58_horn3_9_15GHz_ch17_20_mod_22_5deg_hornang135.dat','2020-06-24/15-46-56_horn3_9_15GHz_ch17_20_mod_45deg_hornang135.dat','2020-06-24/15-22-52_horn3_9_15GHz_ch17_20_mod_67_5deg_hornang135.dat']
		files_2 = ['2020-06-24/11-32-31_horn3_9_15GHz_ch21_24_mod_0deg_hornang_0.dat','2020-06-24/11-40-39_horn3_9_15GHz_ch21_24_mod_22_5deg_hornang_0.dat','2020-06-24/12-07-07_horn3_9_15GHz_ch21_24_mod_45deg_hornang_0.dat','2020-06-24/12-15-29_horn3_9_15GHz_ch21_24_mod_67_5deg_hornang_0.dat','2020-06-24/13-52-44_horn3_9_15GHz_ch21_24_mod_0deg_hornang_45.dat','2020-06-24/13-44-39_horn3_9_15GHz_ch21_24_mod_22_5deg_hornang_45.dat','2020-06-24/13-21-45_horn3_9_15GHz_ch21_24_mod_45deg_hornang_45.dat','2020-06-24/13-13-22_horn3_9_15GHz_ch21_24_mod_67_5deg_hornang_45.dat','2020-06-24/14-16-48_horn3_9_15GHz_ch21_24_mod_0deg_hornang_90.dat','2020-06-24/14-26-18_horn3_9_15GHz_ch21_24_mod_22_5deg_hornang_90.dat','2020-06-24/14-50-45_horn3_9_15GHz_ch21_24_mod_45deg_hornang_90.dat','2020-06-24/15-07-07_horn3_9_15GHz_ch21_24_mod_67_5deg_hornang_90.dat','2020-06-24/16-10-19_horn3_9_15GHz_ch21_24_mod_0deg_hornang135.dat','2020-06-24/16-02-26_horn3_9_15GHz_ch21_24_mod_22_5deg_hornang135.dat','2020-06-24/15-38-53_horn3_9_15GHz_ch21_24_mod_45deg_hornang135.dat','2020-06-24/15-30-45_horn3_9_15GHz_ch21_24_mod_67_5deg_hornang135.dat']
		replace = ['horn3']
		rescale_1 = np.ones(len(files_1))
		rescale_2 = np.ones(len(files_2))
	elif dohorn == 2:
		files_1 =['2020-08-06/11-17-00_horn2_14_22GHz_ch9_12_mod_0deg_hornang_-45.dat','2020-08-06/11-41-55_horn2_14_22GHz_ch9_12_mod_22_5deg_hornang_-45.dat','2020-08-06/11-50-40_horn2_14_22GHz_ch9_12_mod_45deg_hornang_-45.dat','2020-08-06/12-30-59_horn2_14_22GHz_ch9_12_mod_67_5deg_hornang_-45.dat','2020-08-06/12-38-27_horn2_14_22GHz_ch9_12_mod_67_5deg_hornang_0.dat','2020-08-06/13-03-00_horn2_14_22GHz_ch9_12_mod_45deg_hornang_0.dat','2020-08-06/13-11-15_horn2_14_22GHz_ch9_12_mod_22_5deg_hornang_0.dat','2020-08-06/13-35-25_horn2_14_22GHz_ch9_12_mod_0deg_hornang_0.dat','2020-08-06/13-43-00_horn2_14_22GHz_ch9_12_mod_0deg_hornang_45.dat','2020-08-06/14-07-49_horn2_14_22GHz_ch9_12_mod_22_5deg_hornang_45.dat','2020-08-06/14-16-03_horn2_14_22GHz_ch9_12_mod_45deg_hornang_45.dat','2020-08-06/14-39-47_horn2_14_22GHz_ch9_12_mod_67.5deg_hornang_45.dat','2020-08-06/14-47-40_horn2_14_22GHz_ch9_12_mod_67.5deg_hornang_90.dat','2020-08-06/15-11-01_horn2_14_22GHz_ch9_12_mod_45deg_hornang_90.dat','2020-08-06/15-19-08_horn2_14_22GHz_ch9_12_mod_22.5deg_hornang_90.dat','2020-08-06/15-42-33_horn2_14_22GHz_ch9_12_mod_0deg_hornang_90.dat']
		files_2 = ['2020-08-06/11-25-14_horn2_14_22GHz_ch13_16_mod_0deg_hornang_-45.dat','2020-08-06/11-33-51_horn2_14_22GHz_ch13_16_mod_22_5deg_hornang_-45.dat','2020-08-06/12-06-44_horn2_14_22GHz_ch13_16_mod_45deg_hornang_-45.dat','2020-08-06/12-22-39_horn2_14_22GHz_ch13_16_mod_67_5deg_hornang_-45.dat','2020-08-06/12-46-13_horn2_14_22GHz_ch13_16_mod_67_5deg_hornang_0.dat','2020-08-06/12-54-42_horn2_14_22GHz_ch13_16_mod_45deg_hornang_0.dat','2020-08-06/13-19-16_horn2_14_22GHz_ch13_16_mod_22_5deg_hornang_0.dat','2020-08-06/13-27-36_horn2_14_22GHz_ch13_16_mod_0deg_hornang_0.dat','2020-08-06/13-51-00_horn2_14_22GHz_ch13_16_mod_0deg_hornang_45.dat','2020-08-06/13-59-48_horn2_14_22GHz_ch13_16_mod_22_5deg_hornang_45.dat','2020-08-06/14-24-00_horn2_14_22GHz_ch13_16_mod_45deg_hornang_45.dat','2020-08-06/14-32-06_horn2_14_22GHz_ch13_16_mod_67.5deg_hornang_45.dat','2020-08-06/14-55-22_horn2_14_22GHz_ch13_16_mod_67.5deg_hornang_90.dat','2020-08-06/15-03-18_horn2_14_22GHz_ch13_16_mod_45deg_hornang_90.dat','2020-08-06/15-26-48_horn2_14_22GHz_ch13_16_mod_22.5deg_hornang_90.dat','2020-08-06/15-34-56_horn2_14_22GHz_ch13_16_mod_0deg_hornang_90.dat']
		rescale_1 = np.ones(len(files_1))
		rescale_2 = np.ones(len(files_2))
		rescalefactor = 0.858
		rescale_1[0] = rescalefactor
		rescale_1[1] = rescalefactor
		rescale_2[0] = rescalefactor
		rescale_2[1] = rescalefactor
	elif dohorn == 4:
		files_1 = ['2020-08-07/10-32-04_horn4_14_22GHz_ch25_28_mod_0deg_hornang_-45.dat','2020-08-07/10-57-05_horn4_14_22GHz_ch25_28_mod_22_5deg_hornang_-45.dat','2020-08-07/11-05-25_horn4_14_22GHz_ch25_28_mod_45deg_hornang_-45.dat','2020-08-07/11-33-01_horn4_14_22GHz_ch25_28_mod_67_5deg_hornang_-45.dat','2020-08-07/11-40-42_horn4_14_22GHz_ch25_28_mod_67_5deg_hornang_0.dat','2020-08-07/12-07-22_horn4_14_22GHz_ch25_28_mod_45deg_hornang_0.dat','2020-08-07/12-15-58_horn4_14_22GHz_ch25_28_mod_22_5deg_hornang_0.dat','2020-08-07/12-40-29_horn4_14_22GHz_ch25_28_mod_0deg_hornang_0.dat','2020-08-07/12-49-43_horn4_14_22GHz_ch25_28_mod_0deg_hornang_45.dat','2020-08-07/13-16-53_horn4_14_22GHz_ch25_28_mod_22_5deg_hornang_45.dat','2020-08-07/13-25-45_horn4_14_22GHz_ch25_28_mod_45deg_hornang_45.dat','2020-08-07/13-50-44_horn4_14_22GHz_ch25_28_mod_67_5deg_hornang_45.dat','2020-08-07/13-58-42_horn4_14_22GHz_ch25_28_mod_67_5deg_hornang_90.dat','2020-08-07/14-26-59_horn4_14_22GHz_ch25_28_mod_45deg_hornang_90.dat','2020-08-07/14-36-31_horn4_14_22GHz_ch25_28_mod_22_5deg_hornang_90.dat','2020-08-07/15-02-10_horn4_14_22GHz_ch25_28_mod_0deg_hornang_90.dat']
		files_2 = ['2020-08-07/10-40-26_horn4_14_22GHz_ch29_32_mod_0deg_hornang_-45.dat','2020-08-07/10-49-09_horn4_14_22GHz_ch29_32_mod_22_5deg_hornang_-45.dat','2020-08-07/11-13-20_horn4_14_22GHz_ch29_32_mod_45deg_hornang_-45.dat','2020-08-07/11-21-43_horn4_14_22GHz_ch29_32_mod_67_5deg_hornang_-45.dat','2020-08-07/11-48-36_horn4_14_22GHz_ch29_32_mod_67_5deg_hornang_0.dat','2020-08-07/11-58-52_horn4_14_22GHz_ch29_32_mod_45deg_hornang_0.dat','2020-08-07/12-23-52_horn4_14_22GHz_ch29_32_mod_22_5deg_hornang_0.dat','2020-08-07/12-32-23_horn4_14_22GHz_ch29_32_mod_0deg_hornang_0.dat','2020-08-07/12-57-51_horn4_14_22GHz_ch29_32_mod_0deg_hornang_45.dat','2020-08-07/13-08-19_horn4_14_22GHz_ch29_32_mod_22_5deg_hornang_45.dat','2020-08-07/13-34-06_horn4_14_22GHz_ch29_32_mod_45deg_hornang_45.dat','2020-08-07/13-42-37_horn4_14_22GHz_ch29_32_mod_67_5deg_hornang_45.dat','2020-08-07/14-10-35_horn4_14_22GHz_ch29_32_mod_67_5deg_hornang_90.dat','2020-08-07/14-19-05_horn4_14_22GHz_ch29_32_mod_45deg_hornang_90.dat','2020-08-07/14-45-16_horn4_14_22GHz_ch29_32_mod_22_5deg_hornang_90.dat','2020-08-07/14-54-06_horn4_14_22GHz_ch29_32_mod_0deg_hornang_90.dat']
		rescale_1 = np.ones(len(files_1))
		rescale_2 = np.ones(len(files_2))

	numdatasets = len(files_1)
	numchannels = 8
	nummodangs = 4
	numhornangs = 4

	datasets = []
	for i in range(0,numdatasets):
		newset1 = read_quijote_mfi_bandpass_alt(indirectory+files_1[i])
		newset2 = read_quijote_mfi_bandpass_alt(indirectory+files_2[i])
		newset1[1:] *= rescale_1[i]
		newset2[1:] *= rescale_2[i]
		# print(np.shape(newset1))
		# print(np.shape(newset1[1]))
		newset2 = np.asarray([newset2[1],newset2[2],newset2[3],newset2[4]])
		# print(np.shape(newset2))
		datasets.append(np.concatenate((newset1,newset2)))
	# print(np.shape(datasets))
	for i in range(1,numchannels+1):
		# print(i)
		for j in range(0,numdatasets):
			# print(' ' + str(j))
			plt.plot(datasets[j][0],datasets[j][i]/(datasets[j][0][1]-datasets[j][0][0]),'-',label=files_1[j],lw=0.5)
		# plt.legend()
		plt.ylabel('Amplitude (renormalised)')
		plt.xlabel('Frequency (GHz)')
		plt.savefig(plotdir+'horn'+str(dohorn)+'_'+str(i)+'.png')
		plt.clf()

	for j in range(0,numdatasets):
		runname = files_1[j].split('mo')[1].replace('d_','mod_')
		np.savetxt(plotdir+'Pol'+str(dohorn)+'_bandpass_nochange_'+runname,np.transpose([datasets[j][0],datasets[j][4],datasets[j][3],datasets[j][1],datasets[j][2],datasets[j][8],datasets[j][7],datasets[j][5],datasets[j][6]]),header="freq(GHz)  ch1  ch2  ch3  ch4  ch5  ch6  ch7  ch8",fmt="%f",comments='')

	newdataset = np.zeros((8,9,len(datasets[0][0])))
	for i in range(0,4):
		# print(i)
		newdataset[i] = datasets[i] + datasets[i+nummodangs]
	for i in range(8,12):
		# print(i)
		newdataset[i-4] = datasets[i] + datasets[i+nummodangs]
	for j in range(0,8):
		newdataset[j][0] = newdataset[j][0] / 2.0
		
	for i in range(1,numchannels+1):
		# print(i)
		for j in range(0,8):
			newdataset[j][i] = newdataset[j][i]/np.sum(newdataset[j][i])
			if j >= 4:
				label = files_1[j+4]
			else:
				label = files_1[j]
			# print(' ' + str(j))
			plt.plot(newdataset[j][0],newdataset[j][i]/(newdataset[j][0][1]-newdataset[j][0][0]),'-',label=label.split('mod_')[1].replace('.dat','').replace('hornang','horn'),lw=0.5)
		plt.legend()
		plt.ylabel('Amplitude (renormalised)')
		plt.xlabel('Frequency (GHz)')
		plt.savefig(plotdir+'comb_horn'+str(dohorn)+'_'+str(i)+'.png')
		plt.clf()

	for j in range(0,8):
		if j >= 4:
			runname = files_1[j+4].split('mo')[1].replace('d_','sum_mod_')
		else:
			runname = files_1[j].split('mo')[1].replace('d_','sum_mod_')
		np.savetxt(plotdir+'Pol'+str(dohorn)+'_bandpass_'+runname,np.transpose([newdataset[j][0],newdataset[j][4],newdataset[j][3],newdataset[j][1],newdataset[j][2],newdataset[j][8],newdataset[j][7],newdataset[j][5],newdataset[j][6]]),header="freq(GHz)  ch1  ch2  ch3  ch4  ch5  ch6  ch7  ch8",fmt="%f",comments='')
		# np.savetxt(plotdir+'Pol3_bandpass_'+runname,np.transpose([newdataset[j][0],newdataset[j][2],newdataset[j][1],newdataset[j][3],newdataset[j][4],newdataset[j][6],newdataset[j][5],newdataset[j][7],newdataset[j][8]]),header="freq(GHz)  ch1  ch2  ch3  ch4  ch5  ch6  ch7  ch8",fmt="%f",comments='')
