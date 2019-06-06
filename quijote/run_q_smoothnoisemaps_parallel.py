#!/usr/bin/env python
# Original functionality from:
# http://www.shocksolution.com/2010/04/managing-a-pool-of-mpi-processes-with-python-and-pypar/
# History since then:
# Mike Peel   23-Sep-2015   Initial version. Updated to use mpi4py; start up xvfb and call the Matlab pipeline.
# Mike Peel   24-Sep-2015   Add start/end dates, use day_reduce_mod
# Mike Peel   12-Jan-2015   Adding C-BASS South reduction
# Mike Peel   06-Jun-2019   Rework to use for quijote MC runs
from numpy import *
from mpi4py import MPI
import time
import subprocess
import sys
import os
from datetime import date, timedelta as td
from astrocode.fitspectrum.smoothnoisemap import smoothnoisemap
import numpy as np
import healpy as hp
import astropy.io.fits as fits

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
status = MPI.Status()

print("-"*78)
print(" Running on %d cores" % comm.size)
print(" Rank %d" % rank)
print("-"*78)

# Constants
MASTER_PROCESS = 0
WORK_TAG = 1
DIE_TAG = 2


output_resolution = [60.0]
output_nside = [512]#, 256, 128, 64]

# directory = '/Users/mpeel/Documents/maps/quijote_201905/reform/'
# outdirectory = '/Users/mpeel/Documents/maps/quijote_201905/smooth/'
directory = '/net/nas/proyectos/quijote/mikepeel/quijote_201905/reform/'
outdirectory = '/net/nas/proyectos/quijote/mikepeel/quijote_201905/smooth/'

# mfi_l,mfi_111,mfi_113,mfi_217,mfi_219,mfi_311,mfi_313,mfi_417,mfi_419 = np.loadtxt('/net/nas/proyectos/quijote/mikepeel/quijote_mfi_wfs.dat',unpack=True)
prefixes = ['mfi_may2019','mfi_may2019_half1','mfi_may2019_half2']
newprefixes = ['mfi','half1mfi','half2mfi']
newdate = '201905'
numrealisations=1000

mapnumbers = [0,1,2]

### Master Process ###
if rank == MASTER_PROCESS:
	num_processors = comm.size
	print("Master process found " + str(num_processors) + " worker processors.")


	# Work out how many work tasks we have
	numtasks = (len(prefixes)) * (len(output_resolution)) * (len(mapnumbers)) * (6-1)
	work_array = np.arange(0,numtasks+1)
	work_index = 0
	print(work_array)
	# Dispatch jobs to worker processes
	# work_index = 0
	num_completed = 0
	work_size = numtasks+1

	# Start all worker processes
	for i in range(1, min(num_processors, work_size+1)):
		comm.send(work_index, i, tag=WORK_TAG)
		comm.send(work_array[work_index], i)
		print("Sent work index " + str(work_index) + " to processor " + str(i))
		work_index += 1

	# Receive results from each worker, and send it new data
	for i in range(num_processors, work_size+1):
		results = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
		index = status.tag
		proc = status.source
		num_completed += 1
		comm.send(work_index, proc, tag=WORK_TAG)
		comm.send(work_array[work_index], proc)
		print("Sent work index " + str(work_index) + " to processor " + str(proc))
		work_index += 1

	# Get results from remaining worker processes
	while num_completed < work_size-1:
		results = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
		num_completed += 1

	# Shut down worker processes
	for proc in range(1, num_processors):
		print("Stopping worker process " + str(proc))
		comm.send(-1, proc, tag=DIE_TAG)

else:
	### Worker Processes ###
	continue_working = True
	while continue_working:

		work_index =  comm.recv(source=MASTER_PROCESS, tag=MPI.ANY_TAG, status=status)

		if status.tag == DIE_TAG:
			continue_working = False
		else:
			work_array = comm.recv(source=MASTER_PROCESS, tag=MPI.ANY_TAG, status=status)
			work_index = status.tag

		count = 0
		print(work_array)
		print(work_index)
		for k in range(0,len(prefixes)):
			numres = len(output_resolution)
			for i in range(0,numres):
				resolution = "%.2f" % output_resolution[i]
				for m in range(0,len(mapnumbers)):
					# smoothnoisemap(directory, outdirectory, str(output_nside[i])+'_'+resolution+'smoothed_'+newprefixes[k]+'2_17.0_512_'+newdate+'_nhits_'+str(mapnumbers[m]), prefixes[k]+'_nhits_17.0_2.fits',mapnumber=mapnumbers[m],numrealisations=numrealisations,fwhm=output_resolution[i],sigma_0=1.0*np.sqrt(25),nside=output_nside,windowfunction=np.sqrt(mfi_217),usehealpixfits=True)

					if work_array==count:
						print(k)
						print(i)
						print(m)
						print(count)
						print('17-2')
						# smoothnoisemap(directory, outdirectory, str(output_nside[i])+'_'+resolution+'smoothed_'+newprefixes[k]+'2_17.0_512_'+newdate+'_weight_'+str(mapnumbers[m]), prefixes[k]+'_weights_17.0_2.fits',mapnumber=mapnumbers[m],numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=np.sqrt(mfi_217),usehealpixfits=True)
					count += 1
					if work_array==count:
						print(k)
						print(i)
						print(m)
						print(count)
						print('19-2')
						# smoothnoisemap(directory, outdirectory, str(output_nside[i])+'_'+resolution+'smoothed_'+newprefixes[k]+'2_19.0_512_'+newdate+'_weight_'+str(mapnumbers[m]), prefixes[k]+'_weights_19.0_2.fits',mapnumber=mapnumbers[m],numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=np.sqrt(mfi_219),usehealpixfits=True)
					count += 1
					if work_array==count:
						print(k)
						print(i)
						print(m)
						print(count)
						print('11-3')
						# smoothnoisemap(directory, outdirectory, str(output_nside[i])+'_'+resolution+'smoothed_'+newprefixes[k]+'3_11.0_512_'+newdate+'_weight_'+str(mapnumbers[m]), prefixes[k]+'_weights_11.0_3.fits',mapnumber=mapnumbers[m],numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=np.sqrt(mfi_311),usehealpixfits=True)
					count += 1
					if work_array==count:
						print(k)
						print(i)
						print(m)
						print(count)
						print('13-3')
						# smoothnoisemap(directory, outdirectory, str(output_nside[i])+'_'+resolution+'smoothed_'+newprefixes[k]+'3_13.0_512_'+newdate+'_weight_'+str(mapnumbers[m]), prefixes[k]+'_weights_13.0_3.fits',mapnumber=mapnumbers[m],numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=np.sqrt(mfi_313),usehealpixfits=True)
					count += 1
					if work_array==count:
						print(k)
						print(i)
						print(m)
						print(count)
						print('17-4')
						# smoothnoisemap(directory, outdirectory, str(output_nside[i])+'_'+resolution+'smoothed_'+newprefixes[k]+'4_17.0_512_'+newdate+'_weight_'+str(mapnumbers[m]), prefixes[k]+'_weights_17.0_4.fits',mapnumber=mapnumbers[m],numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=np.sqrt(mfi_417),usehealpixfits=True)
					count += 1
					if work_array==count:
						print(k)
						print(i)
						print(m)
						print(count)
						print('19-4')
						# smoothnoisemap(directory, outdirectory, str(output_nside[i])+'_'+resolution+'smoothed_'+newprefixes[k]+'4_19.0_512_'+newdate+'_weight_'+str(mapnumbers[m]), prefixes[k]+'_weights_19.0_4.fits',mapnumber=mapnumbers[m],numrealisations=numrealisations,fwhm=output_resolution[i],nside=output_nside,windowfunction=np.sqrt(mfi_417),usehealpixfits=True)
		print('Max should be ' + str(count))
		comm.send(work_array, dest=MASTER_PROCESS, tag=work_index)
	#### while
#### if worker


MPI.Finalize()

