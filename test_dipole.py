import healpy as hp
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.coordinates import Galactic, FK5
from astropy import units as u
import matplotlib.pyplot as plt
import test_dipole_coord as coord
from astrocode.astroutils import *

nside = 128
npix = hp.pixelfunc.nside2npix(nside)

dipolemap = dipole_pdf(nside, 1.0, 10.0, 10.0, 10.0)
dipolemap *= 100.0
hp.mollview(dipolemap)

plt.savefig('dipole_input1.pdf')


sub = hp.pixelfunc.remove_dipole(dipolemap, fitval=False)
hp.mollview(sub)
plt.savefig('dipole_output1.pdf')


# map_rot = rotate_map(dipolemap, nside, 100)
angle = np.asarray([1.0, 1.0, 1.0])
print(angle)
print(angle.shape)
map_rot = rotate_map3(dipolemap, 100.0)
# map_rot = rotate_map2(dipolemap, angle, 0.0)
# map_rot = rotate_map(dipolemap, nside, 100)
map_rot *= 0.5
print(len(map_rot))
hp.mollview(map_rot)
plt.savefig('dipole_input2.pdf')
sub = hp.pixelfunc.remove_dipole(map_rot, fitval=False)
hp.mollview(sub)
plt.savefig('dipole_output2.pdf')

comb = dipolemap + map_rot
hp.mollview(comb)
plt.savefig('dipole_combination.pdf')
sub = hp.pixelfunc.remove_dipole(comb, fitval=False)
hp.mollview(sub)
plt.savefig('dipole_combination2.pdf')

# alm = hp.map2alm(declinationmap)
# hp.rotator.rotate_alm(alm, )
