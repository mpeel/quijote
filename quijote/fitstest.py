import numpy as np
from astropy.io import fits
from astropy.table import Table

n = np.arange(1000.0)
hdu = fits.PrimaryHDU(n)
hdu.writeto('test1.fits',overwrite=True)

hdul = fits.open('test1.fits')
hdul.info()
print(hdul[0].data)
print(np.shape(hdul[0].data))
hdul[0].data = hdul[0].data[0:100]
print(np.shape(hdul[0].data))
hdul.writeto('test2.fits',overwrite=True)

hdul2 = fits.open('test2.fits')
hdul2.info()
print(np.shape(hdul[0].data))
print(hdul[0].data)

hdul = fits.open('/Users/mpeel/Documents/QUIJOTE/Jupiter/NOMINAL60A-130503-2106-002.ctod')
hdul.info()
print('-.-')
t = Table(hdul[1].data)
print(t)
new_columns = []

for colname in t.columns:
	col = t[colname]
	print(col.shape)
	if len(col.shape) == 3 or len(col.shape) == 32:
		col = col[:, :100, :]
	if len(col.shape) > 1:
		if col.shape[1] > 1000:
			col = col[:,:100]
	new_columns.append(col)
new_t = Table(new_columns, names=t.columns) 
print('.-.')
print(new_t)
hdul[1].data = np.array(new_t)
hdul.writeto('test_ctod.fits',overwrite=True)

hdul2 = fits.open('test_ctod.fits')
hdul2.info()
print(np.shape(hdul2[0].data))
print(hdul2[1].data)
