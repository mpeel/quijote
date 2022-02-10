from compare_polang import *
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors

# This is needed if you want to write out lots of plots
mpl.rcParams['agg.path.chunksize'] = 10000


# compare_polang(prefix='mfi',date='201905')
#compare_polang(prefix='mfi302',date='201904',use_variance=False)
#compare_polang(prefix='mfi',date='201907')
# compare_polang(prefix='mfi',date='201911',datestr='nov2019')
# compare_polang(prefix='mfi',date='201911',datestr='nov2019_half2')

# compare_polang(prefix='mfi',date='202003',datestr='mar2020')
# compare_polang(prefix='mfi',date='202004',datestr='apr2020')
# compare_polang(prefix='mfi',date='202004b',datestr='apr2020b')


# Run many and do a comparison

prefixes = ['mfi_apr2020b_pwv1', 'mfi_apr2020b_allatonce', 'mfi_apr2020b_altsample1', 'mfi_apr2020b_altsample2', 'mfi_apr2020b_daynight1', 'mfi_apr2020b_daynight2', 'mfi_apr2020b_fivesample1', 'mfi_apr2020b_fivesample2', 'mfi_apr2020b_half1', 'mfi_apr2020b_half2', 'mfi_apr2020b_halfring1', 'mfi_apr2020b_halfring2', 'mfi_apr2020b_period1', 'mfi_apr2020b_period2', 'mfi_apr2020b_period5', 'mfi_apr2020b_period6', 'mfi_apr2020b_pwv2', 'mfi_apr2020b_ring1', 'mfi_apr2020b_ring2', 'mfi_apr2020b_sample1', 'mfi_apr2020b_sample2', 'mfi_apr2020b_tbem1', 'mfi_apr2020b_tbem2', 'mfi_apr2020b_twosample1', 'mfi_apr2020b_twosample2']
valset = []
prefixes.sort()
# for i in range(0,1):
# for i in range(0,len(prefixes)):
# 	newprefix = prefixes[i].replace('_apr2020b','')

# 	try:
# 		vals = compare_polang(prefix=newprefix,date='202004b_jk',datestr='apr2020b',indirectory='/Volumes/Toshiba5TB2/mfi/quijote_202004b_jk/')
# 	except:
# 		vals = [np.zeros(6), np.zeros(6), np.zeros(6), np.zeros(6)]
# 	valset.append(vals)

print(valset)

valset = [[[-2.777, -3.439,  0.978, -1.425, -1.304,  2.068], [ 0.599,  0.877,  0.299,  0.313,  0.311,  0.538], [-0.566, -7.290, -0.155, -1.957,  0.281,  4.793], [ 0.173,  0.281,  0.041,  0.057,  0.098,  0.178]], [[-2.902, -5.106,  0.913, -1.348, -1.374,  1.331], [ 0.639,  0.898,  0.294,  0.321,  0.377,  0.581], [-0.892, -8.531, -0.197, -2.163,  0.208,  4.469], [ 0.169,  0.262,  0.041,  0.057,  0.095,  0.174]], [[-2.905, -4.193,  0.920, -1.376, -1.451,  2.261], [ 0.632,  0.878,  0.288,  0.326,  0.391,  0.573], [-1.461, -8.274, -0.133, -1.985,  0.308,  4.401], [ 0.166,  0.260,  0.041,  0.056,  0.097,  0.173]], [[-3.456, -6.250,  1.030, -1.612, -1.515,  1.495], [ 0.716,  0.925,  0.311,  0.364,  0.477,  0.646], [-2.184, -11.641, -0.260, -2.265,  0.706,  3.042], [ 0.154,  0.209,  0.041,  0.056,  0.095,  0.170]], [[-2.343, -7.488,  0.681, -1.102, -1.548,  1.224], [ 0.664,  0.894,  0.315,  0.332,  0.386,  0.611], [-1.769, -12.340, -0.402, -1.899, -0.495,  2.665], [ 0.154,  0.228,  0.041,  0.056,  0.094,  0.161]], [[-3.742, -6.454,  0.870, -1.318, -1.258,  1.301], [ 0.668,  0.908,  0.295,  0.334,  0.395,  0.596], [-1.518, -11.056, -0.148, -1.562,  0.033,  4.117], [ 0.163,  0.237,  0.041,  0.056,  0.095,  0.170]], [[-2.431, -4.651,  0.949, -1.378, -1.332,  1.759], [ 0.656,  0.895,  0.310,  0.357,  0.413,  0.595], [-1.386, -8.142, -0.090, -2.097,  0.478,  4.039], [ 0.161,  0.241,  0.041,  0.056,  0.096,  0.171]], [[-3.729, -6.697,  0.870, -1.386, -1.590,  1.653], [ 0.715,  0.920,  0.303,  0.346,  0.464,  0.609], [-1.857, -9.413, -0.070, -1.873,  0.526,  3.636], [ 0.155,  0.218,  0.041,  0.056,  0.097,  0.172]], [[-3.283, -6.256,  0.876, -1.296, -1.302,  1.197], [ 0.638,  0.895,  0.295,  0.320,  0.395,  0.602], [-1.606, -11.773, -0.217, -1.930,  0.097,  3.773], [ 0.159,  0.227,  0.041,  0.056,  0.094,  0.164]], [[-2.798, -3.341,  0.968, -1.355, -1.291,  1.921], [ 0.594,  0.876,  0.299,  0.316,  0.312,  0.546], [-0.586, -7.360, -0.156, -1.940,  0.282,  4.807], [ 0.173,  0.280,  0.041,  0.057,  0.098,  0.178]], [[-26.191, -23.016, -2.375, -109.770, -62.393, -109.316], [ 1.724,  1.713,  6.800,  7.119,  5.193,  9.270], [-20.362, -32.552, -14.853, -19.961, -0.133, -46.740], [ 1.840,  1.726,  0.479,  2.482,  8.464,  2.088]], [[ 0.000,  0.000,  0.000,  0.000,  0.000,  0.000], [ 0.000,  0.000,  0.000,  0.000,  0.000,  0.000], [ 0.000,  0.000,  0.000,  0.000,  0.000,  0.000], [ 0.000,  0.000,  0.000,  0.000,  0.000,  0.000]], [[-1.115, -4.435,  0.089, -1.142,  4.697,  0.167], [ 0.558,  0.863,  0.367,  0.465,  0.636,  0.717], [ 0.138, -9.160, -0.572, -0.261,  23.052,  3.859], [ 0.167,  0.270,  0.040,  0.055,  14.364,  0.106]], [[ 0.000,  0.000,  0.000,  0.000,  0.000,  0.000], [ 0.000,  0.000,  0.000,  0.000,  0.000,  0.000], [ 0.000,  0.000,  0.000,  0.000,  0.000,  0.000], [ 0.000,  0.000,  0.000,  0.000,  0.000,  0.000]], [[-4.358, -5.515,  0.574, -1.739, -1.641,  0.957], [ 0.701,  0.901,  0.305,  0.354,  0.348,  0.577], [-2.486, -10.697, -0.153, -1.875,  0.211,  3.652], [ 0.153,  0.215,  0.039,  0.055,  0.095,  0.175]], [[-2.805, -5.105,  0.868, -1.336, -1.636,  1.672], [ 0.661,  0.904,  0.292,  0.344,  0.445,  0.583], [-0.149, -7.967, -0.115, -0.821,  1.058,  4.028], [ 0.156,  0.235,  0.041,  0.055,  0.095,  0.172]], [[-4.285, -9.111,  0.920, -1.312, -1.335,  1.189], [ 0.708,  0.918,  0.295,  0.342,  0.428,  0.648], [-2.997, -14.993, -0.134, -2.197, -0.421,  3.130], [ 0.152,  0.210,  0.040,  0.057,  0.095,  0.163]], [[-2.707, -6.514,  0.938, -1.176, -1.697,  1.618], [ 0.694,  0.918,  0.302,  0.336,  0.426,  0.612], [-1.502, -9.180, -0.213, -1.815,  0.115,  4.483], [ 0.159,  0.228,  0.041,  0.056,  0.095,  0.171]], [[-3.352, -6.370,  0.725, -1.407, -1.077,  1.742], [ 0.688,  0.913,  0.307,  0.340,  0.409,  0.597], [-1.777, -12.016, -0.172, -2.183,  0.186,  4.156], [ 0.156,  0.226,  0.041,  0.056,  0.096,  0.170]], [[-3.009, -5.280,  1.022, -1.349, -1.321,  1.626], [ 0.649,  0.899,  0.293,  0.346,  0.404,  0.588], [-1.369, -9.447, -0.197, -1.983,  0.111,  4.574], [ 0.167,  0.261,  0.041,  0.057,  0.096,  0.173]], [[-2.484, -4.337,  0.888, -1.327, -1.401,  1.637], [ 0.651,  0.897,  0.300,  0.314,  0.381,  0.586], [-0.266, -9.345, -0.242, -1.927,  0.361,  4.476], [ 0.168,  0.261,  0.041,  0.057,  0.097,  0.174]], [[-2.626, -8.321,  0.794, -1.377, -1.676,  1.443], [ 0.651,  0.918,  0.285,  0.344,  0.414,  0.614], [-0.412, -11.412, -0.130, -1.971, -0.063,  3.739], [ 0.154,  0.228,  0.041,  0.056,  0.095,  0.170]], [[-4.089, -5.341,  0.951, -1.346, -1.601,  1.223], [ 0.715,  0.889,  0.319,  0.357,  0.449,  0.604], [-2.544, -11.211, -0.496, -2.182,  0.273,  3.243], [ 0.157,  0.206,  0.041,  0.056,  0.095,  0.164]], [[-2.801, -3.554,  0.929, -1.349, -1.547,  1.107], [ 0.631,  0.891,  0.300,  0.322,  0.381,  0.586], [-1.506, -8.234, -0.139, -2.017,  0.165,  4.222], [ 0.165,  0.256,  0.041,  0.057,  0.095,  0.173]], [[-2.724, -4.680,  0.805, -1.458, -1.469,  1.618], [ 0.639,  0.888,  0.303,  0.343,  0.406,  0.587], [-1.124, -8.967, -0.330, -2.235,  0.190,  4.231], [ 0.168,  0.257,  0.041,  0.056,  0.097,  0.173]]]

# Std
m217 = []
m219 = []
m311 = []
m313 = []
m417 = []
m419 = []
w217 = []
w219 = []
w311 = []
w313 = []
w417 = []
w419 = []

threshold = 10
for i in range(0,len(valset)):
	if valset[i][0][0] < threshold and valset[i][0][0] > -threshold and valset[i][0][0] != 0:
		m217.append(valset[i][0][0])
	if valset[i][0][1] < threshold and valset[i][0][1] > -20 and valset[i][0][1] != 0:
		m219.append(valset[i][0][1])
	if valset[i][0][2] < threshold and valset[i][0][2] > -threshold and valset[i][0][2] != 0:
		m311.append(valset[i][0][2])
	if valset[i][0][3] < threshold and valset[i][0][3] > -threshold and valset[i][0][3] != 0:
		m313.append(valset[i][0][3])
	if valset[i][0][4] < threshold and valset[i][0][4] > -threshold and valset[i][0][4] != 0:
		m417.append(valset[i][0][4])
	if valset[i][0][5] < threshold and valset[i][0][5] > -threshold and valset[i][0][5] != 0:
		m419.append(valset[i][0][5])

	if valset[i][2][0] < threshold and valset[i][2][0] > -threshold and valset[i][2][0] != 0:
		w217.append(valset[i][2][0])
	if valset[i][2][1] < threshold and valset[i][2][1] > -20 and valset[i][2][1] != 0:
		w219.append(valset[i][2][1])
	if valset[i][2][2] < threshold and valset[i][2][2] > -threshold and valset[i][2][2] != 0:
		w311.append(valset[i][2][2])
	if valset[i][2][3] < threshold and valset[i][2][3] > -threshold and valset[i][2][3] != 0:
		w313.append(valset[i][2][3])
	if valset[i][2][4] < threshold and valset[i][2][4] > -threshold and valset[i][2][4] != 0:
		w417.append(valset[i][2][4])
	if valset[i][2][5] < threshold and valset[i][2][5] > -threshold and valset[i][2][5] != 0:
		w419.append(valset[i][2][5])
print(np.std(m217))
print(np.std(m219))
print(np.std(m311))
print(np.std(m313))
print(np.std(m417))
print(np.std(m419))
print(np.std(w217))
print(np.std(w219))
print(np.std(w311))
print(np.std(w313))
print(np.std(w417))
print(np.std(w419))
exit()


# Plot
plt.figure(figsize=(7.2*1.5, 5.4*1.5), dpi=80)
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
colours.remove('w')

# for i in range(0,1):
for i in range(0,len(prefixes)):
	if np.sum(valset[i]) != 0.0:
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
plt.ylim(-20,10)
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig('plot_poldifference_jk2.pdf')
