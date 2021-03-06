from compare_polang import *
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors

# This is needed if you want to write out lots of plots
mpl.rcParams['agg.path.chunksize'] = 10000

colours = []
for c in mcolors.BASE_COLORS:
	colours.append(c)
for c in mcolors.TABLEAU_COLORS:
	colours.append(c)
for c in mcolors.CSS4_COLORS:
	colours.append(c)

# print(colours)
# exit()

# compare_polang(prefix='mfi',date='201905')
#compare_polang(prefix='mfi302',date='201904',use_variance=False)
#compare_polang(prefix='mfi',date='201907')
# compare_polang(prefix='mfi',date='201911',datestr='nov2019')
# compare_polang(prefix='mfi',date='201911',datestr='nov2019_half2')

# compare_polang(prefix='mfi',date='202003',datestr='mar2020')
# compare_polang(prefix='mfi',date='202004',datestr='apr2020')
# compare_polang(prefix='mfi',date='202004b',datestr='apr2020b')
# compare_polang(prefix='mfi_recalib3',date='202004b',datestr='apr2020b')
# compare_polang(prefix='mfi_recalib4',date='202004b',datestr='apr2020b')
# compare_polang(prefix='mfi_recalib5',date='202004b',datestr='apr2020b')
# compare_polang(prefix='mfi_recalib6',date='202004b',datestr='apr2020b')
# compare_polang(prefix='mfi_recalib7',date='202004b',datestr='apr2020b')
# compare_polang(prefix='mfi_recalib8',date='202004b',datestr='apr2020b')
# compare_polang(prefix='mfi_recalib9',date='202004b',datestr='apr2020b')
# compare_polang(prefix='mfi',date='202011',datestr='nov2020')

# This was the cross-check from 2020 to 2021 changes.
# compare_polang(prefix='mfi',date='202011',datestr='nov2020',indirectory='/Volumes/Toshiba5TB2/maps/quijote_202011/')
# exit()

# Compare 2021 with polang_err
# compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=False,polang_err_threshold=2.0,planckmap='npipe',applyoffsets=False,mapprefix='_period2')
# compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=False,polang_err_threshold=2.0,planckmap='npipe',applyoffsets=False,mapprefix='_period5')
# compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=False,polang_err_threshold=2.0,planckmap='npipe',applyoffsets=False,mapprefix='_period6')

# compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=True,polang_err_threshold=2.0,planckmap='npipe',applyoffsets='wmap',mapprefix='')
# compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=False,polang_err_threshold=2.0,planckmap='npipe',applyoffsets='wmap',mapprefix='_period2')
# compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=True,polang_err_threshold=2.0,planckmap='npipe',applyoffsets='wmap',mapprefix='_period5')
# compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=True,polang_err_threshold=2.0,planckmap='npipe',applyoffsets='wmap',mapprefix='
# _period6')
# compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=False,polang_err_threshold=0.0,planckmap='npipe',applyoffsets='wmap')

# All sky
# compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=False,planckmap='npipe',applyoffsets=False,outextraext='wmap_npipe_nooff')
# compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=False,planckmap='2018',applyoffsets=False,outextraext='wmap_2018_nooff')
# compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=False,planckmap='2015',applyoffsets=False,outextraext='wmap_2015_nooff')
# OLD - with offsets
# compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=False,planckmap='npipe',applyoffsets='planck',outextraext='npipe')
# compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=False,planckmap='2018',applyoffsets='planck',outextraext='2018')
# compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=False,planckmap='2015',applyoffsets='planck',outextraext='2015')

# compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=True,polang_err_threshold=2.0,planckmap='npipe',applyoffsets=False,outextraext='2deg_npipe')
# compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=True,polang_err_threshold=2.0,planckmap='2018',applyoffsets=False,outextraext='2deg_2018')
# compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=True,polang_err_threshold=2.0,planckmap='2015',applyoffsets=False,outextraext='2deg_2015')
#
# compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=True,polang_err_threshold=3.0,planckmap='npipe',applyoffsets=False,outextraext='3deg_npipe')
# compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=True,polang_err_threshold=3.0,planckmap='2018',applyoffsets=False,outextraext='3deg_2018')
# compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=True,polang_err_threshold=3.0,planckmap='2015',applyoffsets=False,outextraext='3deg_2015')
#
#
# compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=True,polang_err_threshold=5.0,planckmap='npipe',applyoffsets=False,outextraext='5deg_npipe')
# compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=True,polang_err_threshold=5.0,planckmap='2018',applyoffsets=False,outextraext='5deg_2018')
# compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=True,polang_err_threshold=5.0,planckmap='2015',applyoffsets=False,outextraext='5deg_2015')

# compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=True,polang_err_threshold=2.0,planckmap='npipe',applyoffsets='wmap',outextraext='wmap')
# compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=True,polang_err_threshold=2.0,planckmap='npipe',applyoffsets='planck',outextraext='npipe')
# compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=True,polang_err_threshold=2.0,planckmap='2018',applyoffsets='planck',outextraext='2018')
# compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=True,polang_err_threshold=2.0,planckmap='2015',applyoffsets='planck',outextraext='2015')

# With static masks
# inputmask='/Users/mpeel/Documents/maps/quijote_masks/masc_polangle_1.fits'
# compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=False,planckmap='npipe',applyoffsets=False,outextraext='wmap_npipe_nooff_mask1_m2',staticmask=True,inputmask=inputmask)
# compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=False,planckmap='2018',applyoffsets=False,outextraext='wmap_2018_nooff_mask1_m2',staticmask=True,inputmask=inputmask)
# compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=False,planckmap='2015',applyoffsets=False,outextraext='wmap_2015_nooff_mask1_m2',staticmask=True,inputmask=inputmask)
# inputmask='/Users/mpeel/Documents/maps/quijote_masks/masc_polangle_2.fits'
# compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=False,planckmap='npipe',applyoffsets=False,outextraext='wmap_npipe_nooff_mask2_w2',staticmask=True,inputmask=inputmask)
# compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=False,planckmap='2018',applyoffsets=False,outextraext='wmap_2018_nooff_mask2_m2',staticmask=True,inputmask=inputmask)
# compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=False,planckmap='2015',applyoffsets=False,outextraext='wmap_2015_nooff_mask2_m2',staticmask=True,inputmask=inputmask)
# inputmask='/Users/mpeel/Documents/maps/quijote_masks/masc_polangle_3.fits'
# compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=False,planckmap='npipe',applyoffsets=False,outextraext='wmap_npipe_nooff_mask3_m2',staticmask=True,inputmask=inputmask)
# compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=False,planckmap='2018',applyoffsets=False,outextraext='wmap_2018_nooff_mask3_m2',staticmask=True,inputmask=inputmask)
# compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=False,planckmap='2015',applyoffsets=False,outextraext='wmap_2015_nooff_mask3_m2',staticmask=True,inputmask=inputmask)

# With static masks and simulations
inputmask='/Users/mpeel/Documents/maps/quijote_masks/masc_polangle_1.fits'
compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=False,planckmap='npipe',applyoffsets=False,outextraext='wmap_npipe_nooff_mask1_sim',staticmask=True,inputmask=inputmask,simnoise=True)
compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=False,planckmap='2018',applyoffsets=False,outextraext='wmap_2018_nooff_mask1_sim',staticmask=True,inputmask=inputmask,simnoise=True)
compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=False,planckmap='2015',applyoffsets=False,outextraext='wmap_2015_nooff_mask1_sim',staticmask=True,inputmask=inputmask,simnoise=True)
inputmask='/Users/mpeel/Documents/maps/quijote_masks/masc_polangle_2.fits'
compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=False,planckmap='npipe',applyoffsets=False,outextraext='wmap_npipe_nooff_mask2_sim',staticmask=True,inputmask=inputmask,simnoise=True)
compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=False,planckmap='2018',applyoffsets=False,outextraext='wmap_2018_nooff_mask2_sim',staticmask=True,inputmask=inputmask,simnoise=True)
compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=False,planckmap='2015',applyoffsets=False,outextraext='wmap_2015_nooff_mask2_sim',staticmask=True,inputmask=inputmask,simnoise=True)
inputmask='/Users/mpeel/Documents/maps/quijote_masks/masc_polangle_3.fits'
compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=False,planckmap='npipe',applyoffsets=False,outextraext='wmap_npipe_nooff_mask3_sim',staticmask=True,inputmask=inputmask,simnoise=True)
compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=False,planckmap='2018',applyoffsets=False,outextraext='wmap_2018_nooff_mask3_sim',staticmask=True,inputmask=inputmask,simnoise=True)
compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=False,planckmap='2015',applyoffsets=False,outextraext='wmap_2015_nooff_mask3_sim',staticmask=True,inputmask=inputmask,simnoise=True)
# inputmask='/Users/mpeel/Documents/maps/quijote_masks/masc_polangle_1.fits'
# compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=False,planckmap='npipe',applyoffsets=False,outextraext='wmap_npipe_nooff_mask1_m2_sim',staticmask=True,inputmask=inputmask,simnoise=True)
# compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=False,planckmap='2018',applyoffsets=False,outextraext='wmap_2018_nooff_mask1_m2_sim',staticmask=True,inputmask=inputmask,simnoise=True)
# compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=False,planckmap='2015',applyoffsets=False,outextraext='wmap_2015_nooff_mask1_m2_sim',staticmask=True,inputmask=inputmask,simnoise=True)
# inputmask='/Users/mpeel/Documents/maps/quijote_masks/masc_polangle_2.fits'
# compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=False,planckmap='npipe',applyoffsets=False,outextraext='wmap_npipe_nooff_mask2_w2_sim',staticmask=True,inputmask=inputmask,simnoise=True)
# compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=False,planckmap='2018',applyoffsets=False,outextraext='wmap_2018_nooff_mask2_m2_sim',staticmask=True,inputmask=inputmask,simnoise=True)
# compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=False,planckmap='2015',applyoffsets=False,outextraext='wmap_2015_nooff_mask2_m2_sim',staticmask=True,inputmask=inputmask,simnoise=True)
# inputmask='/Users/mpeel/Documents/maps/quijote_masks/masc_polangle_3.fits'
# compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=False,planckmap='npipe',applyoffsets=False,outextraext='wmap_npipe_nooff_mask3_m2_sim',staticmask=True,inputmask=inputmask,simnoise=True)
# compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=False,planckmap='2018',applyoffsets=False,outextraext='wmap_2018_nooff_mask3_m2_sim',staticmask=True,inputmask=inputmask,simnoise=True)
# compare_polang(prefix='mfi',date='202103',datestr='apr2021',newformat=True,indirectory='/Users/mpeel/Documents/maps/quijote_202103_tqu_v1.5_noise_v1.0_newwf/',use_polang_err=False,planckmap='2015',applyoffsets=False,outextraext='wmap_2015_nooff_mask3_m2_sim',staticmask=True,inputmask=inputmask,simnoise=True)

exit()


# Run many and do a comparison

prefixes = ['mfi_apr2020b_pwv1', 'mfi_apr2020b_allatonce', 'mfi_apr2020b_altsample1', 'mfi_apr2020b_altsample2', 'mfi_apr2020b_daynight1', 'mfi_apr2020b_daynight2', 'mfi_apr2020b_fivesample1', 'mfi_apr2020b_fivesample2', 'mfi_apr2020b_half1', 'mfi_apr2020b_half2', 'mfi_apr2020b_halfring1', 'mfi_apr2020b_halfring2', 'mfi_apr2020b_period1', 'mfi_apr2020b_period2', 'mfi_apr2020b_period5', 'mfi_apr2020b_period6', 'mfi_apr2020b_pwv2', 'mfi_apr2020b_ring1', 'mfi_apr2020b_ring2', 'mfi_apr2020b_sample1', 'mfi_apr2020b_sample2', 'mfi_apr2020b_tbem1', 'mfi_apr2020b_tbem2', 'mfi_apr2020b_twosample1', 'mfi_apr2020b_twosample2']
valset = []
prefixes.sort()
# for i in range(0,1):
for i in range(0,len(prefixes)):
	newprefix = prefixes[i].replace('_apr2020b','')

	try:
		vals = compare_polang(prefix=newprefix,date='202004b_jk',datestr='apr2020b',indirectory='/Volumes/Toshiba5TB2/mfi/quijote_202004b_jk/')
	except:
		vals = [np.zeros(6), np.zeros(6), np.zeros(6), np.zeros(6)]
	valset.append(vals)

print(valset)

plt.figure(figsize=(7.2*2, 5.4*2), dpi=80)
ax = plt.subplot(111)
lines = {'linestyle': 'None'}
plt.rc('lines', **lines)
index = np.arange(1,7)
# symbols = ['ob','sb','<b']
# symbols2 = ['or','sr','<r']

colours = []
for c in mcolors.BASE_COLORS:
	colours.append(c)
for c in mcolors.TABLEAU_COLORS:
	colours.append(c)
for c in mcolors.CSS4_COLORS:
	colours.append(c)

# for i in range(0,1):
for i in range(0,len(prefixes)):
	newprefix = prefixes[i].replace('_apr2020b','')
	plt.plot(index,np.zeros(len(index)),'-k')
	plt.errorbar(index-0.15,valset[i][0],yerr=valset[i][1],marker='<',color=colours[i])
	plt.errorbar(index-0.08,valset[i][2],yerr=valset[i][3],marker='o',color=colours[i],label=newprefix)

# plt.errorbar(index-0.15,valset[i][0],yerr=diff_311_med_unc,fmt='ob',label='Diff to 311')
# plt.errorbar(index-0.15,diff_311_wgt,yerr=diff_311_wgt_unc,fmt='sb')
# plt.errorbar(index-0.15,diff_311_mask_med,yerr=diff_311_mask_med_unc,fmt='^b')

# plt.errorbar(index-0.08,diff_wmap_wgt,yerr=diff_wmap_wgt_unc,fmt='sr')
# plt.errorbar(index-0.08,diff_wmap_mask_med,yerr=diff_wmap_mask_med_unc,fmt='^r')
# plt.errorbar(index-0.08,diff_wmap_mask_wgt,yerr=diff_wmap_mask_wgt_unc,fmt='<r')
# plt.errorbar(index+0.08,diff_planck_med,yerr=diff_planck_med_unc,fmt='og',label='Diff to Planck')
# plt.errorbar(index+0.08,diff_planck_wgt,yerr=diff_planck_wgt_unc,fmt='sg')
# plt.errorbar(index+0.08,diff_planck_mask_med,yerr=diff_planck_mask_med_unc,fmt='^g')
# plt.errorbar(index+0.08,diff_planck_mask_wgt,yerr=diff_planck_mask_wgt_unc,fmt='<g')

# plt.errorbar(index,diff_taua,yerr=diff_taua_unc,fmt='ok',label='Tau A wide')
# plt.errorbar(index,diff_raster,yerr=diff_raster_unc,fmt='sk',label='Tau A raster')
plt.ylabel('Angle difference [deg]')

box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig('plot_poldifference_jk.png')
