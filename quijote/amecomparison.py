import healpy as hp
import matplotlib.pyplot as plt

filename = '/Users/mpeel/Documents/Quijote/AME/201903_comparison/combined_gal_hlat_BestPolBonaldi_free_synch_LIN_AMEAmpO.fits'
mapdata = hp.read_map(filename,field=0)
print len(mapdata)
hp.mollview(mapdata,norm='hist',xsize=4000)
plt.savefig(filename+'_plot.pdf')

hp.mollview(mapdata,min=0,max=750, xsize=4000)
plt.savefig(filename+'_cut.pdf')

hp.mollview(mapdata,min=0,max=750, xsize=4000,coord=['G','C'])
plt.savefig(filename+'_cut_celestial.pdf')

hp.gnomview(mapdata,rot=[10.68, 41.27],coord=['G','C'],title='M31')
plt.savefig(filename+'_m31.pdf')

hp.gnomview(mapdata,rot=[160.26,-18.62],title='Perseus')
plt.savefig(filename+'_perseus.pdf')

hp.gnomview(mapdata,rot=[160.26,-12.5],title='California nebula')
plt.savefig(filename+'_california.pdf')

hp.gnomview(mapdata,rot=[92,-37],title='Pegasus plume',reso=5.0,min=0.0)
plt.savefig(filename+'_plume.pdf')

# hp.gnomview(mapdata,rot=[209.01,-19.38],title='Orion')
# plt.savefig(filename+'_orion.pdf')

hp.gnomview(mapdata,rot=[195.2,-12.2],title='Lambda Orionis',reso=5.0)
plt.savefig(filename+'_lorionis.pdf')

