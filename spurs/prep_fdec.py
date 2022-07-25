# We have fdec-subtracted Planck and WMAP maps, but they lost the variance maps, which we need in one file for the weighted map.
# So, read them in, recombine them, write them back out.
# This won't override any files - to regenerate, delete the old version.
# MP, 25 Jul 2022
import healpy as hp
import astropy.io.fits as fits

# 2015
indir_fdec = '/Users/mpeel/Documents/maps/quijote_fdec/planck2015_tqu_v1.5_noise_v1.0/'
files_fdec = ['512_60.0smoothed_PlanckR2fullbeamNoise_44.1_1024_2015_mKCMBunits_nofdec_ring.fits','512_60.0smoothed_PlanckR2fullbeamNoise_28.4_1024_2015_mKCMBunits_nofdec_ring.fits','512_60.0smoothed_PlanckR2fullbeamNoise_70.4_2048_2015_mKCMBunits_nofdec_ring.fits']
indir_raw = '/Users/mpeel/Documents/maps/planck2015_tqu_v1.5_noise_v1.0/'
files_raw = ['512_60.0smoothed_PlanckR2fullbeamNoise_44.1_1024_2015_mKCMBunits.fits','512_60.0smoothed_PlanckR2fullbeamNoise_28.4_1024_2015_mKCMBunits.fits','512_60.0smoothed_PlanckR2fullbeamNoise_70.4_2048_2015_mKCMBunits.fits']
outdir = '/Users/mpeel/Documents/maps/planck2015_tqu_v1.5_noise_v1.0/'
outfiles = ['512_60.0smoothed_PlanckR2fullbeamFdecNoise_44.1_1024_2015_mKCMBunits.fits','512_60.0smoothed_PlanckR2fullbeamFdecNoise_28.4_1024_2015_mKCMBunits.fits','512_60.0smoothed_PlanckR2fullbeamFdecNoise_70.4_2048_2015_mKCMBunits.fits']
# 2018
indir_fdec = '/Users/mpeel/Documents/maps/quijote_fdec/planck2018_tqu_v1.5_noise_v1.0/'
files_fdec = ['512_60.0smoothed_PlanckR3fullbeamNoise_44.1_1024_2018_mKCMBunits_nofdec_ring.fits','512_60.0smoothed_PlanckR3fullbeamNoise_28.4_1024_2018_mKCMBunits_nofdec_ring.fits','512_60.0smoothed_PlanckR3fullbeamNoise_70.4_1024_2018_mKCMBunits_nofdec_ring.fits']
indir_raw = '/Users/mpeel/Documents/maps/planck2018_tqu_v1.5_noise_v1.0_10k/'
files_raw = ['512_60.0smoothed_PlanckR3fullbeamNoise_44.1_1024_2018_mKCMBunits.fits','512_60.0smoothed_PlanckR3fullbeamNoise_28.4_1024_2018_mKCMBunits.fits','512_60.0smoothed_PlanckR3fullbeamNoise_70.4_1024_2018_mKCMBunits.fits']
outdir = '/Users/mpeel/Documents/maps/planck2018_tqu_v1.5_noise_v1.0_10k/'
outfiles = ['512_60.0smoothed_PlanckR3fullbeamFdecNoise_44.1_1024_2018_mKCMBunits.fits','512_60.0smoothed_PlanckR3fullbeamFdecNoise_28.4_1024_2018_mKCMBunits.fits','512_60.0smoothed_PlanckR3fullbeamFdecNoise_70.4_1024_2018_mKCMBunits.fits']
# 2020
indir_fdec = '/Users/mpeel/Documents/maps/quijote_fdec/planck2020_tqu_v1.5_noise_v1.0/'
files_fdec = ['512_60.0smoothed_PlanckR4fullbeamnodpNoise_28.4_1024_2020_mKCMBunits_nofdec_ring.fits','512_60.0smoothed_PlanckR4fullbeamnodpNoise_44.1_1024_2020_mKCMBunits_nofdec_ring.fits','512_60.0smoothed_PlanckR4fullbeamnodpNoise_70.4_1024_2020_mKCMBunits_nofdec_ring.fits']
indir_raw = '/Users/mpeel/Documents/maps/planck2020_tqu_v1.5_noise_v1.0_10k/'
files_raw = ['512_60.0smoothed_PlanckR4fullbeamnodpNoise_28.4_1024_2020_mKCMBunits.fits','512_60.0smoothed_PlanckR4fullbeamnodpNoise_44.1_1024_2020_mKCMBunits.fits','512_60.0smoothed_PlanckR4fullbeamnodpNoise_70.4_1024_2020_mKCMBunits.fits']
outdir = '/Users/mpeel/Documents/maps/planck2020_tqu_v1.5_noise_v1.0_10k/'
outfiles = ['512_60.0smoothed_PlanckR4fullbeamnodpFdecNoise_28.4_1024_2020_mKCMBunits.fits','512_60.0smoothed_PlanckR4fullbeamnodpFdecNoise_44.1_1024_2020_mKCMBunits.fits','512_60.0smoothed_PlanckR4fullbeamnodpFdecNoise_70.4_1024_2020_mKCMBunits.fits']
# WMAP
# indir_fdec = '/Users/mpeel/Documents/maps/quijote_fdec/wmap9_tqu_v1.5_noise_v1.0/'
# files_fdec = ['512_60.0smoothed_wmap9beamNoise_33.0_512_2013_mKCMBunits_nofdec_ring.fits','512_60.0smoothed_wmap9beamNoise_22.8_512_2013_mKCMBunits_nofdec_ring.fits','512_60.0smoothed_wmap9beamNoise_40.7_512_2013_mKCMBunits_nofdec_ring.fits']
# indir_raw = '/Users/mpeel/Documents/maps/wmap9_tqu_v1.5_noise_v1.0_10k/'
# files_raw = ['512_60.0smoothed_wmap9beamNoise_33.0_512_2013_mKCMBunits.fits','512_60.0smoothed_wmap9beamNoise_22.8_512_2013_mKCMBunits.fits','512_60.0smoothed_wmap9beamNoise_40.7_512_2013_mKCMBunits.fits']
# outdir = '/Users/mpeel/Documents/maps/wmap9_tqu_v1.5_noise_v1.0_10k/'
# outfiles = ['512_60.0smoothed_wmap9beamFdecNoise_33.0_512_2013_mKCMBunits.fits','512_60.0smoothed_wmap9beamFdecNoise_22.8_512_2013_mKCMBunits.fits','512_60.0smoothed_wmap9beamFdecNoise_40.7_512_2013_mKCMBunits.fits']

def get_header_val(hdr,search):
	for i in range(0,len(hdr)):
		if search in hdr[i][0]:
			return hdr[i][1]
	return ''

for i in range(0,len(files_fdec)):
	fdec_maps = hp.read_map(indir_fdec+files_fdec[i],field=None)
	raw_maps,h = hp.read_map(indir_raw+files_raw[i],field=None,h=True)
	print(h)
	column_names = [get_header_val(h,'TTYPE1'), get_header_val(h,'TTYPE2'), get_header_val(h,'TTYPE3'), get_header_val(h,'TTYPE4'), get_header_val(h,'TTYPE5'), get_header_val(h,'TTYPE6'), get_header_val(h,'TTYPE7')]
	column_units = [get_header_val(h,'TUNIT1'), get_header_val(h,'TUNIT2'), get_header_val(h,'TUNIT3'), get_header_val(h,'TUNIT4'), get_header_val(h,'TUNIT5'), get_header_val(h,'TUNIT6'), get_header_val(h,'TUNIT7')]
	extra_header = [get_header_val(h, 'COMMENT')]
	try:
		hp.write_map(outdir+outfiles[i],[fdec_maps[0], fdec_maps[1], fdec_maps[2], raw_maps[3], raw_maps[4], raw_maps[5],raw_maps[6]],column_names=column_names,column_units=column_units,extra_header=extra_header)
	except:
		pass
