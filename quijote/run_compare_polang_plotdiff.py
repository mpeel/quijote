import matplotlib.pyplot as plt
import numpy as np

index = np.arange(1,7)
diff_311_med = [-4.8, -7.7, 0.0, -0.6, -0.7, -3.5]
diff_311_med_unc = [0.4, 0.5, 0.0, 0.3, 0.4, 0.4]
diff_311_wgt = [-1.64, -8.17, 0.0, -0.59, 1.65, 1.79]
diff_311_wgt_unc = [0.13, 0.16, 0.0, 0.06, 0.09, 0.13]
diff_311_mask_med = [-3.1, -6.0, 0.0, -1.5, -0.4, 0.3]
diff_311_mask_med_unc = [0.6, 0.9, 0.0, 0.3, 0.4, 0.6]
diff_311_mask_wgt = [-0.48, -8.5, 0.0, -1.0, 2.63, 3.94]
diff_311_mask_wgt_unc = [0.17, 0.3, 0.0, 0.07, 0.1, 0.17]

diff_wmap_med = [-7.0, -10.8, -1.5, -2.6, -2.2, -4.8]
diff_wmap_med_unc = [0.4, 0.5, 0.3, 0.3, 0.4, 0.4]
diff_wmap_wgt = [-3.23, -10.49, -1.60, -2.34, -0.87, -0.66]
diff_wmap_wgt_unc = [0.11, 0.13, 0.03, 0.05, 0.07, 0.11]
diff_wmap_mask_med = [-2.5, -5.3, 0.4, -1.1, 0.1, 0.8]
diff_wmap_mask_med_unc = [0.6, 0.9, 0.3, 0.3, 0.3, 0.6]
diff_wmap_mask_wgt = [-0.29, -8.8, -0.68, -1.77, 1.74, 3.31]
diff_wmap_mask_wgt_unc = [0.16, 0.3, 0.04, 0.05, 0.09, 0.17]

diff_planck_med = [-6.5, -10.3, -1.1, -2.6, -2.2, -4.8]
diff_planck_med_unc = [0.4, 0.5, 0.3, 0.4, 0.4, 0.4]
diff_planck_wgt = [-3.93, -10.8, -0.85, -1.80, 0.07, 0.61]
diff_planck_wgt_unc = [0.11, 0.13, 0.03, 0.05, 0.07, 0.11]
diff_planck_mask_med = [-3.1, -5.4, -0.2, -1.7, -0.6, 0.2]
diff_planck_mask_med_unc = [0.6, 0.9, 0.3, 0.3, 0.3, 0.6]
diff_planck_mask_wgt = [-0.17, -8.6, -0.05, -1.09, 2.43, 4.53]
diff_planck_mask_wgt_unc = [0.16, 0.2, 0.04, 0.05, 0.09, 0.16]

diff_taua = [0.3, -5.7, -2.1, -0.3, 2.6, -2.2]
diff_taua_unc = [0.7, 1.0, 1.1, 0.9, 1.7, 2.3]
diff_raster = [1.2, -3.3, -0.4, 1.3, 3.2, -0.8]
diff_raster_unc = [0.5, 2.0, 0.8, 0.9, 1.0, 1.1]

lines = {'linestyle': 'None'}
plt.rc('lines', **lines)
plt.plot(index,np.zeros(len(index)),'-k')
plt.errorbar(index-0.15,diff_311_med,yerr=diff_311_med_unc,fmt='ob',label='Diff to 311')
plt.errorbar(index-0.15,diff_311_wgt,yerr=diff_311_wgt_unc,fmt='sb')
plt.errorbar(index-0.15,diff_311_mask_med,yerr=diff_311_mask_med_unc,fmt='^b')
plt.errorbar(index-0.15,diff_311_mask_wgt,yerr=diff_311_mask_wgt_unc,fmt='<b')
plt.errorbar(index-0.08,diff_wmap_med,yerr=diff_wmap_med_unc,fmt='or',label='Diff to WMAP')
plt.errorbar(index-0.08,diff_wmap_wgt,yerr=diff_wmap_wgt_unc,fmt='sr')
plt.errorbar(index-0.08,diff_wmap_mask_med,yerr=diff_wmap_mask_med_unc,fmt='^r')
plt.errorbar(index-0.08,diff_wmap_mask_wgt,yerr=diff_wmap_mask_wgt_unc,fmt='<r')
plt.errorbar(index+0.08,diff_planck_med,yerr=diff_planck_med_unc,fmt='og',label='Diff to Planck')
plt.errorbar(index+0.08,diff_planck_wgt,yerr=diff_planck_wgt_unc,fmt='sg')
plt.errorbar(index+0.08,diff_planck_mask_med,yerr=diff_planck_mask_med_unc,fmt='^g')
plt.errorbar(index+0.08,diff_planck_mask_wgt,yerr=diff_planck_mask_wgt_unc,fmt='<g')

plt.errorbar(index,diff_taua,yerr=diff_taua_unc,fmt='ok',label='Tau A wide')
plt.errorbar(index,diff_raster,yerr=diff_raster_unc,fmt='sk',label='Tau A raster')
plt.ylabel('Angle difference [deg]')
plt.legend()
plt.savefig('plot_poldifference.png')