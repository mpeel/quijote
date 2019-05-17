import healpy as hp
import matplotlib.pyplot as plt

directory = '/Users/mpeel/Documents/Quijote/AME/201904_comparison/'
indir = directory+'files/'
outdir = directory+'plots/'
filenames = ['combined_gal_hlat100000_symmetrized_hasreich_updated_quijote2019_ALL_LIN_AMEAmpO.fits','combined_gal_hlat100000_symmetrized_hasreich_updated_quijote2019_ALL_LIN_AMEAmpO1.fits','combined_gal_hlat100000_symmetrized_hasreich_updated_quijote2019_ALL_LIN_AMEAmpO2.fits','combined_gal_hlat100000_symmetrized_hasreich_updated_quijote2019_ALL_LIN_FFAmpO.fits','combined_gal_hlat100000_symmetrized_hasreich_updated_quijote2019_ALL_LIN_SynchAmpO.fits','combined_gal_hlat100000_symmetrized_hasreich_updated_quijote2019_ALL_LIN_ameFreqO.fits','combined_gal_hlat_BestPolBonaldi_free_synch_LIN_AMEAmpO.fits','combined_gal_hlat_BestPolBonaldi_free_synch_LIN_FFAmpO.fits','combined_gal_hlat_BestPolBonaldi_free_synch_LIN_SynchAmpO.fits','combined_gal_hlat_BestPolBonaldi_free_synch_LIN_SynchIndO.fits','combined_gal_hlat_BestPolBonaldi_free_synch_LIN_ameFreqO.fits','combined_gal_hlat_BestPolBonaldi_free_synch_LIN_m60O.fits','combined_gal_hlat_BestPolFixed2Bonaldi_ALL_NOLIN_AMEAmpO.fits','combined_gal_hlat_BestPolFixed2Bonaldi_ALL_NOLIN_FFAmpO.fits','combined_gal_hlat_BestPolFixed2Bonaldi_ALL_NOLIN_SynchAmpO.fits','combined_gal_hlat_BestPolFixed2Bonaldi_ALL_NOLIN_ameFreqO.fits','combined_gal_hlat_BestPolFixed2Bonaldi_ALL_NOLIN_m60O.fits','combined_gal_hlat_galpropBonaldi_symmetrized_hasreich_updated_quijote2019_ALL_LIN_AMEAmpO.fits','combined_gal_hlat_galpropBonaldi_symmetrized_hasreich_updated_quijote2019_ALL_LIN_FFAmpO.fits','combined_gal_hlat_galpropBonaldi_symmetrized_hasreich_updated_quijote2019_ALL_LIN_SynchAmpO.fits','combined_gal_hlat_galpropBonaldi_symmetrized_hasreich_updated_quijote2019_ALL_LIN_ameFreqO.fits','combined_gal_hlat_galpropBonaldi_symmetrized_hasreich_updated_quijote2019_ALL_LIN_m60O.fits','combined_gal_hlat_galpropMIAC_ameFreq2_symmetrized_hasreich_updated_quijote2019_ALL_LIN_AMEAmpO.fits','combined_gal_hlat_galpropMIAC_ameFreq2_symmetrized_hasreich_updated_quijote2019_ALL_LIN_FFAmpO.fits','combined_gal_hlat_galpropMIAC_ameFreq2_symmetrized_hasreich_updated_quijote2019_ALL_LIN_SynchAmpO.fits','combined_gal_hlat_galpropMIAC_ameFreq2_symmetrized_hasreich_updated_quijote2019_ALL_LIN_ameFreqO.fits','combined_gal_hlat_galpropMIAC_ameFreq2_symmetrized_hasreich_updated_quijote2019_ALL_LIN_sigmaO.fits','combined_gal_hlat_galpropMIAC_symmetrized_hasreich_updated_quijote2019_ALL_LIN_AMEAmpO.fits','combined_gal_hlat_galpropMIAC_symmetrized_hasreich_updated_quijote2019_ALL_LIN_FFAmpO.fits','combined_gal_hlat_galpropMIAC_symmetrized_hasreich_updated_quijote2019_ALL_LIN_SynchAmpO.fits','combined_gal_hlat_galpropMIAC_symmetrized_hasreich_updated_quijote2019_ALL_LIN_ameFreqO.fits','combined_gal_hlat_galpropMIAC_symmetrized_hasreich_updated_quijote2019_ALL_LIN_sigmaO.fits','combined_gal_hlat_planckModel_symmetrized_ALL_LIN_AMEAmpO1.fits','combined_gal_hlat_planckModel_symmetrized_ALL_LIN_AMEAmpO2.fits','combined_gal_hlat_planckModel_symmetrized_ALL_LIN_FFAmpO.fits','combined_gal_hlat_planckModel_symmetrized_ALL_LIN_SynchAmpO.fits','combined_gal_hlat_planckModel_symmetrized_ALL_LIN_ameFreqO.fits']
for filename in filenames:
	mapdata = hp.read_map(indir+filename,field=0)

	# Check for bad pixels
	mapdata[mapdata < -1e20] = hp.pixelfunc.UNSEEN
	print(len(mapdata))
	hp.mollview(mapdata,norm='hist',xsize=4000)
	plt.savefig(outdir+filename+'_plot.pdf')

	hp.mollview(mapdata,min=0,max=750, xsize=4000)
	plt.savefig(outdir+filename+'_cut.pdf')

	hp.mollview(mapdata,min=0,max=750, xsize=4000,coord=['G','C'])
	plt.savefig(outdir+filename+'_cut_celestial.pdf')

	hp.gnomview(mapdata,rot=[10.68, 41.27],coord=['G','C'],title='M31')
	plt.savefig(outdir+filename+'_m31.pdf')

	hp.gnomview(mapdata,rot=[160.26,-18.62],title='Perseus')
	plt.savefig(outdir+filename+'_perseus.pdf')

	hp.gnomview(mapdata,rot=[160.26,-12.5],title='California nebula')
	plt.savefig(outdir+filename+'_california.pdf')

	hp.gnomview(mapdata,rot=[92,-37],title='Pegasus plume',reso=5.0,min=0.0)
	plt.savefig(outdir+filename+'_plume.pdf')

	# hp.gnomview(mapdata,rot=[209.01,-19.38],title='Orion')
	# plt.savefig(filename+'_orion.pdf')

	hp.gnomview(mapdata,rot=[195.2,-12.2],title='Lambda Orionis',reso=5.0)
	plt.savefig(outdir+filename+'_lorionis.pdf')

