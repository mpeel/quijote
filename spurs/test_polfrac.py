import healpy as hp
import numpy as np
import matplotlib.pyplot as plt

# Location of data
basedir = '/Users/mpeel/Documents/maps/'
outdir = '/Users/mpeel/Desktop/'
mfi_weighted = 'quijote_202103_tqu_v1.5_noise_v1.0_weighted/256_60.0smoothed_quijotecombwei10_tqu_v1.5_noise_v1.0_-3.1_combine.fits'
commander_sync = 'planck_commander2015/commander_freq_maps/commander_sync_10.0.fits'


# Weighted MFI map
mfiw = hp.read_map(basedir+mfi_weighted,field=None)
mfiw_pol = np.sqrt(mfiw[1]**2+mfiw[2]**2)
mfiw_pol[mfiw[1]==0] = hp.UNSEEN

plotmax_p = 1.5

hp.mollview(mfiw_pol,min=0,max=plotmax_p,cmap='jet',unit='mK CMB',title='')#,title='MFI weighted polarised intensity'
plt.savefig(outdir+'mfi_combine.pdf')
plt.clf()
plt.close()

cmd = hp.read_map(basedir+commander_sync)

polfracmax = 0.7

hp.mollview(cmd,max=plotmax_p/polfracmax,cmap='jet')
plt.savefig(outdir+'commander_10ghz.png')
plt.clf()
plt.close()

hp.mollview(mfiw_pol/cmd,min=0,max=polfracmax,cmap='jet')
plt.savefig(outdir+'polfrac_commander_mfiweighted.png')
plt.clf()
plt.close()
