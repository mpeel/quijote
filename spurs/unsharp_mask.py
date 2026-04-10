import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
from astrocode.astroutils import *
import matplotlib
from astrocode.polfunc import *

mfi_weighted = 'quijote_202103_tqu_v1.5_noise_v1.0_weighted/256_60.0smoothed_quijotecombwei10_tqu_v1.5_noise_v1.0_-3.1_combine.fits'
basedir = '/Users/mpeel/Documents/maps/'
#outdir=basedir+'quijote_202103_tqu_v1.5_noise_v1.0_newwf/spurs/'
outdir=basedir+'2026_loops_spurs/'

fontsize = 18
matplotlib.rcParams.update({'font.size':fontsize})

inputmap, h = hp.read_map(basedir+mfi_weighted,field=None,h=True)
# print(h)
# exit()
polmap = np.sqrt(inputmap[1]**2+inputmap[2]**2)
# unsharpmap_qt = unsharpmap(polmap,15)#np.sqrt(10**2-1**2))
# unsharpmap_qt -= np.min(unsharpmap_qt)
# hp.mollview(unsharpmap_qt,norm='asinh')
# # hp.mollview(unsharpmap_qt,min=-0.5,max=0.5)
# plt.savefig(outdir+'_test_unsharp_15.pdf')
# # hp.write_map(outdir+'_test_unsharp_15.fits',unsharpmap_qt)

test_debias = debias_p_as(inputmap[1],inputmap[2],np.sqrt(inputmap[3]),np.sqrt(inputmap[5]))
test_debias[0][np.where(inputmap[3]>1e10)]=hp.UNSEEN
# hp.mollview(test_debias[1],max=1.0,cmap='jet')
# hp.mollview(test_debias[0],norm='hist')
# plt.savefig(outdir+'_test_deb.pdf')
# hp.mollview(inputmap[1],norm='hist')
# plt.savefig(outdir+'_test_deb_in1.pdf')
# hp.mollview(inputmap[2],norm='hist')
# plt.savefig(outdir+'_test_deb_in2.pdf')
test_plot = np.sqrt(inputmap[3])
test_plot[np.where(inputmap[3]>1e10)]=hp.UNSEEN
# hp.mollview(test_plot,norm='hist')
# plt.savefig(outdir+'_test_deb_in3.pdf')
test_plot = np.sqrt(inputmap[5])
test_plot[np.where(inputmap[5]>1e10)]=hp.UNSEEN
# hp.mollview(np.sqrt(test_plot),norm='hist')
# plt.savefig(outdir+'_test_deb_in5.pdf')
# hp.mollview(polmap,norm='hist')
# plt.savefig(outdir+'_test_deb_in.pdf')

# testmap = bgfmap(polmap, np.sqrt(10.0**2-1.0**2), 10)
testmap = bgfmap(test_debias[0], np.sqrt(10.0**2-1.0**2), 10)
# hp.mollview(testmap,norm='hist',cmap='jet')
# plt.savefig(outdir+'_test_bgfmap_10_deb.pdf')
# hp.mollview(testmap,max=1.0,cmap='jet')
# plt.savefig(outdir+'_test_bgfmap_10_linear_deb.pdf')
testmap[np.where(inputmap[3]>1e10)] = hp.UNSEEN
hp.mollview(testmap,max=1.0,cmap='jet',unit='mK CMB',title='')
plt.savefig(outdir+'fig2_mfi_bgfmap.pdf')


# hp.mollview(polmap,max=1.0,cmap='jet')
# plt.savefig(outdir+'_test_bgfmap_input.pdf')


# This approach doesn't work very well!
# testq = bgfmap(inputmap[1], np.sqrt(10.0**2-1.0**2), 10)
# testu = bgfmap(inputmap[2], np.sqrt(10.0**2-1.0**2), 10)
# hp.mollview(np.sqrt(testq**2+testu**2),max=1.0,cmap='jet')
# plt.savefig(outdir+'_test_bgfmap_qu_10_linear.pdf')
# test_debias = debias_p_mas(testq,testu,inputmap[3],inputmap[5])
# hp.mollview(test_debias[1],max=1.0,cmap='jet')
# plt.savefig(outdir+'_test_bgfmap_qu_10_linear_deb.pdf')
