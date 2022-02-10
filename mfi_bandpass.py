#!/usr/bin/env python
# -*- coding: utf-8  -*-
#
# Reduce the MFI polarisation tests from May 2020
# 
# Version history:
#
# 13-May-2020  M. Peel       Started
# 26-May-2020  M. Peel		 Reworking to fit offsets and frequencies
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits

# Work out the frequencies - using the centre frequency of the bin, hence the extra deltafreq/2.
def get_frequency_range(start, end, numsamples):
	deltafreq = (end-start)/numsamples
	return np.arange(start,end,deltafreq) + deltafreq/2.0

def get_offset_wide_from_narrow(bandpass,freqs,prevbandpass,prevfreqs):
	# frequencies = get_frequency_range(freqs[0],freqs[1],len(bandpass))
	prevfrequencies = get_frequency_range(prevfreqs[0],prevfreqs[1],len(prevbandpass))
	prevfrequencies = prevfrequencies[0:len(prevbandpass)]
	deltafreq = (freqs[1]-freqs[0])/len(bandpass)
	deltaprevfreq = (prevfreqs[1]-prevfreqs[0])/len(prevbandpass)
	stretch = deltafreq/deltaprevfreq
	resample_size = int(len(prevbandpass)/stretch)
	resample_frequencies = get_frequency_range(prevfreqs[0],prevfreqs[1],resample_size)
	resample_bandpass = np.interp(resample_frequencies,prevfrequencies,prevbandpass)
	# plt.plot(prevfrequencies,prevbandpass,'b-')
	# plt.plot(resample_frequencies,resample_bandpass,'r-')
	# plt.savefig('/Users/mpeel/Documents/QUIJOTE/MFI/2020-05_MFI_measurements/_plots/test.png')
	# plt.clf()
	# exit()
	# print(len(prevbandpass))
	# print(stretch)
	# print(resample_size)
	# exit()

	# We now have a narrow bandpass with known frequencies, and a wider bandpass with unknown offset
	# Both on the same frequency step scale
	# Run through and calculate the optimal offset
	# We want to extend the bandpass array in case we have an offset close to the end
	print(np.shape(bandpass))
	bandpass2 = np.concatenate((bandpass,bandpass))
	diffs = np.zeros(len(bandpass))
	print(np.shape(resample_bandpass))
	print(np.shape(bandpass2))
	for k in range(0,len(bandpass)):
		diffs[k] = np.std(bandpass2[k:k+resample_size]-resample_bandpass[:])
	offset = np.argmin(diffs)
	# That is the offset in the array. Now adjust it for the different starting frequencies.
	offset_in_freq = int((prevfreqs[0]-freqs[0])/deltafreq)
	print(offset_in_freq)
	offset -= offset_in_freq
	# We might have gone negative, in that case add the length of the bandpass to correct.
	if offset < 0:
		offset += len(bandpass)
	print(offset)
	return offset


# def autocorr(x):
#     n = x.size
#     norm = (x - np.mean(x))
#     result = np.correlate(norm, norm, mode='same')
#     acorr = result[n//2 + 1:] / (x.var() * np.arange(n-1, n//2, -1))
#     lag = np.abs(acorr).argmax() + 1
#     r = acorr[lag-1]        
#     if np.abs(r) > 0.5:
#       print('Appears to be autocorrelated with r = {}, lag = {}'. format(r, lag))
#     else: 
#       print('Appears to be not autocorrelated')
#     return r, lag

indirectory = '/Users/mpeel/Documents/QUIJOTE/MFI/2020-05_MFI_measurements/'
plotdir = indirectory+'_plots/'
# runname = '2020_05_21'
# files = ['horn2_vna_bandpass_cal_17-19GHz-200521-1156.tod','horn2_vna_bandpass_cal-200521-1148-000.tod','horn3_vna_bandpass_cal_10-14GHz-200521-1100-000.tod','horn3_vna_bandpass-200521-1043-000.tod','horn4_vna_bandpass_cal_17-19GHz-200521-1210-000.tod','horn4_vna_bandpass_cal-200521-1204-000.tod']
runname = '2020_05_22_0deg'
files = ['horn2_vna_bandpass_17_19GHz_0-200522-1104-000.tod','horn2_vna_bandpass_wide_0-200522-1058-000.tod','horn3_vna_bandpass_10_14GHz_0-200522-1210-000.tod','horn3_vna_bandpass_wide_0-200522-1206-000.tod','horn4_vna_bandpass_17_19GHz_0-200522-1120-000.tod','horn4_vna_bandpass_wide_0-200522-1116-000.tod','horn4_vna_bandpass_3dB_17_19GHz_0-200522-1149-000.tod','horn4_vna_bandpass_3dB_wide_0-200522-1152-000.tod']
# runname = '2020_05_22_22p5deg'
# files = ['horn2_vna_bandpass_17_19GHz_22_5-200522-1106-000.tod','horn2_vna_bandpass_wide_22_5-200522-1109-000.tod','horn3_vna_bandpass_10_14GHz_22_5-200522-1211-000.tod','horn3_vna_bandpass_wide_22_5-200522-1214-000.tod','horn4_vna_bandpass_17_19GHz_22.5-200522-1122-000.tod','horn4_vna_bandpass_wide_22.5-200522-1125-000.tod','horn4_vna_bandpass_3dB_17_19GHz_22.5-200522-1147-000.tod','horn4_vna_bandpass_3dB_wide_22.5-200522-1137-000.tod']
horns = [2,2,3,3,4,4,4,4]
rectype = ['n','w','n','w','n','w','n','w']
freqrange = [[17.0,19.0],[9.0,21.0],[10.0,14.0],[9.0,21.0],[17.0,19.0],[9.0,21.0],[17.0,19.0],[9.0,21.0]]
# horns = [3,3,3,3]
numchannels = 8
# nominallength = 7000
nominallength = 4000

runname = '2020_06_18'
files = ['horn3_bandpass_11_13GHz_0deg_cal_horn-200618-1533-000.tod','horn3_bandpass_9_15GHz_0deg_cal_horn-200618-1516-000.tod','horn3_bandpass_11_13GHz_22_5deg_cal_horn-200618-1532-000.tod','horn3_bandpass_9_15GHz_22_5deg_cal_horn-200618-1519-000.tod','horn3_bandpass_11_13GHz_45deg_cal_horn-200618-1531-000.tod','horn3_bandpass_9_15GHz_45deg_cal_horn-200618-1521-000.tod','horn3_bandpass_11_13GHz_67_5deg_cal_horn-200618-1529-000.tod','horn3_bandpass_9_15GHz_67_5deg_cal_horn-200618-1524-000.tod']
horns = [3,3,3,3,3,3,3,3]
rectype = ['n','w','n','w','n','w','n','w']
freqrange = [[11,13],[9,15],[11,13],[9,15],[11,13],[9,15],[11,13],[9,15]]
nominallength = 7000


runname = '2020_06_19'
files = ['BEM3_11_13GHz_config1-200619-1403-000.tod','BEM3_9_15GHz_config1-200619-1359-000.tod','BEM3_11_13GHz_config2-200619-1406-000.tod','BEM3_9_15GHz_config2-200619-1409-000.tod']
horns = [3,3,3,3]
rectype = ['n','w','n','w']
freqrange = [[11,13],[9,15],[11,13],[9,15]]
nominallength = 7000


fold = np.ones(len(horns),dtype=int)*int(nominallength)
print(fold)
# fold = [6889,7004,6938,6994,6889,6994,7000,7000]
# fold = [4000,4000,4000,400,4889,4994,4000,4000]
offset = np.zeros(len(horns),dtype=int)
num_test_samples = int(nominallength/5)
sample1 = 1
sample2 = 3


prevbandpass = [] # For later
# for i in range(0,1):
for i in range(0,len(horns)):
	start = numchannels*(horns[i]-1)
	end = numchannels*horns[i]
	inputfits = fits.open(indirectory+files[i])
	# print(inputfits[1].header)
	data = inputfits[1].data.field('data')[0][:,start:end]
	print(np.shape(data))
	# data_diff = np.diff(data,axis=0)
	# print(np.shape(data_diff))

	# Find the folding frequency starting with the input assumption, and minimising the std of the difference.
	## Previous attempt to use an fft, which wasn't accurate enough.
	# fft = np.fft.rfft(data[:,0])
	# freq = np.fft.rfftfreq(len(data[:,0]))
	# index = np.argmax(np.abs(fft[50:]))+50
	# print(index)
	# print(freq[index]*2)
	# fold[i] = int(np.abs((1.0/freq[index])))
	# print(fold[i]*2)
	# if fold[i] < 5000:
	# 	fold[i] = fold[i]*2
	# plt.plot(np.abs(fft[10:1000]))
	# plt.savefig(indirectory+'fft.png')
	# print(fft)
	## ... and to use autocorr, which didn't work.
	# r, lag = autocorr(data[:,0])
	diffs = np.zeros(2*num_test_samples+1)
	for k in range(0,2*num_test_samples+1):
		diffs[k] = np.std(data[sample1*(fold[i]-num_test_samples+k):(sample1+1)*(fold[i]-num_test_samples+k),0] - data[sample2*(fold[i]-num_test_samples+k):(sample2+1)*(fold[i]-num_test_samples+k),0])
		# plt.plot(data[sample1*(fold[i]-num_test_samples+k):(sample1+1)*(fold[i]-num_test_samples+k),0] - data[sample2*(fold[i]-num_test_samples+k):(sample2+1)*(fold[i]-num_test_samples+k),0])
		# plt.savefig(plotdir+'test_'+str(k)+'.png')
		# plt.clf()
	# print(diffs)
	# print(np.argmin(diffs))
	# And find the folding position
	fold[i] = fold[i]-num_test_samples+np.argmin(diffs)

	# Find the offset. This is nominally found by the maximum difference between samples, which should be at the start/end of a narrow scan.
	test = data[sample1*fold[i]:(sample1+1)*fold[i],0]
	# Do a bit of smoothing in case we have noise
	numsmooth = 10
	test = np.convolve(test, np.ones((numsmooth,))/numsmooth, mode='valid')
	difference = np.abs(np.diff(test))
	plt.plot(data[sample1*fold[i]:(sample1+1)*fold[i],0])
	plt.plot(difference)
	plt.savefig(plotdir+'_test_'+str(i)+'.png')
	plt.clf()
	# exit()

	offset[i] = np.argmax(difference)+1
	offset[i] = 0
	print(offset[i])
	# ... if we have a previous narrower bandpass, then we can try bootsrapping from that.
	if rectype[i] == 'w':
		offset[i] = get_offset_wide_from_narrow(data[sample1*fold[i]:(sample1+1)*fold[i],0],freqrange[i],prevbandpass[0],freqrange[i-1])

	frequencies = get_frequency_range(freqrange[i][0],freqrange[i][1],fold[i])
	frequencies = frequencies[0:fold[i]]

	comb_bandpass = np.zeros((numchannels,fold[i]))
	numsamples = len(data)
	numsets = int(numsamples/fold[i])-1
	print('Number of sets is ')
	print(numsets)
	print(fold[i])

	# This is the main loop, where we go through the multiple measurements and combine them to a single bandpass.
	for k in range(0,numchannels):
		plt.title('Horn ' +str(horns[i]) + ' channel ' + str(k))
		bandpass = np.zeros(fold[i])
		for j in range(0,numsets):
		# for j in range(0,5):
			print(j)
			measurement = data[offset[i]+j*fold[i]:offset[i]+(j+1)*fold[i],k]
			measurement -= np.median(measurement)
			plt.plot(frequencies,measurement,',')
			bandpass = np.add(bandpass, measurement)
		bandpass = bandpass/numsets
		comb_bandpass[k][:] = bandpass
		plt.plot(frequencies,bandpass,'r-',linewidth=1)
		plt.ylabel('Amplitude')
		plt.xlabel('Frequency [GHz]')
		plt.savefig(plotdir+files[i]+'_'+str(k)+'.png')
		plt.clf()

	# Save out the combined bandpass
	if rectype[i] == 'w':
		np.savetxt(plotdir+runname+'_Pol'+str(horns[i])+'_bandpass_'+str(i)+'.dat',np.transpose([frequencies,comb_bandpass[0],comb_bandpass[1],comb_bandpass[2],comb_bandpass[3],comb_bandpass[4],comb_bandpass[5],comb_bandpass[6],comb_bandpass[7]]),header="freq(GHz)  ch1  ch2  ch3  ch4  ch5  ch6  ch7  ch8",fmt="%f",comments='')

	# If we have a narrow scan, we're going to want to save it for the next loop.
	if rectype[i] == 'n':
		prevbandpass = comb_bandpass.copy()

print(fold)
print(offset)
