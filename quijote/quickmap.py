import healpy as hp
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

# Configuration
input_file = '/Users/mpeel/Desktop/NOMINAL60A-130503-2106-000.ctod'
nside = 64

# Set up the map
npix=12*nside*nside
outmap = np.zeros(npix)

# Read in the data
hdul = fits.open(input_file)

# Select what we want to map
hornindex = 0
gl = np.asarray(hdul[1].data.field('GL')[hornindex])
gb = np.asarray(hdul[1].data.field('GB')[hornindex])
data = np.asarray(hdul[1].data.field('DATA')[hornindex])
print(len(gl))
print(len(gb))
print(len(data))

# Convert to Healpix pixels
healpix_pixel = np.asarray(hp.ang2pix(nside, (np.pi/2)-gb*np.pi/180.0, gl*np.pi/180.0))
print(len(healpix_pixel))

# Bin the data into the map
for pixnum in np.unique(healpix_pixel[healpix_pixel > 0]):
	print(pixnum)
	if pixnum < npix:
		outmap[pixnum] = np.sum(data[np.where(healpix_pixel==pixnum)])

# Write out the map and plot it
hp.write_map(input_file.replace('.ctod','_map.fits'),outmap,overwrite=True)
hp.mollview(outmap)
plt.savefig(input_file.replace('.ctod','_map.png'))
