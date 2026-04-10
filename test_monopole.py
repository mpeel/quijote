import healpy as hp
import numpy as np

nside = 16
npix = hp.pixelfunc.nside2npix(nside)
testmap = np.ones(npix)
alms = hp.map2alm(testmap)
cls = hp.alm2cl(alms)
print(cls)
cls2 = hp.anafast(testmap)
print(cls2)
new_cls = np.zeros(len(cls))
new_cls[0] = 1.0
new_alms = hp.sphtfunc.synalm(new_cls)
new_map = hp.alm2map(new_alms,nside=nside)
print(new_map)
