import healpy as hp
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.coordinates import Galactic, FK5
from astropy import units as u
import matplotlib.pyplot as plt
import test_dipole_coord as coord

# Function from Federica, 10 December 2019
def rotate_map(map, nside, dphi):
   rot_map = np.zeros(hp.nside2npix(nside))
   #pix2ang
   phi, theta = hp.pixelfunc.pix2ang(nside, np.arange(hp.nside2npix(nside)), nest=False, lonlat=True)
   #rotphi
   phi = phi+dphi/np.cos(theta)
   #ang2pix
   ipix = hp.pixelfunc.ang2pix(nside, phi, theta, nest=False, lonlat=True)
   #reordeer
   rot_map[np.arange(hp.nside2npix(nside))] = map[ipix]
   return(rot_map)

def rotate_map3(map, dphi):
	R = hp.rotator.Rotator(rot=[dphi,0.0],deg=True)
	return R.rotate_map_alms(map)


# From https://astro.pages.rwth-aachen.de/astrotools/_modules/healpytools.html#dipole_pdf
def dipole_pdf(nside, a, v, y=None, z=None, pdf=True):
    """
    Probability density function of a dipole. Returns 1 + a * cos(theta) for all pixels in hp.nside2npix(nside).

    :param nside: nside of the healpy map
    :param a: amplitude of the dipole, 0 <= a <= 1, automatically clipped
    :param v: either (x, y, z) vector of the pixel center(s) or only x-coordinate
    :param y: y-coordinate(s) of the center
    :param z: z-coordinate(s) of the center
    :param pdf: if false, return unnormalized map (modulation around 1)
    :return: weights
    """
    assert (a >= 0.) and (a <= 1.), "Dipole amplitude must be between 0 and 1"
    a = np.clip(a, 0., 1.)
    if y is not None and z is not None:
        v = np.array([v, y, z], dtype=np.float)

    # normalize to one
    direction = v / np.sqrt(np.sum(v ** 2))
    npix = hp.nside2npix(nside)
    v = np.array(hp.pix2vec(nside, np.arange(npix)))
    cos_angle = np.sum(v.T * direction, axis=1)
    dipole_map = 1 + a * cos_angle

    if pdf is True:
        dipole_map /= np.sum(dipole_map)

    return dipole_map

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
