import numpy as np
from astropy.io import fits
from astropy.table import Table

# Demo code for fits files
# n = np.arange(1000.0)
# a = ('test','test','test')
# col1 = fits.Column(name="numbers",array=n,format="E")
# col2 = fits.Column(name="strings",array=a,format="10A")
# hdu = fits.BinTableHDU.from_columns([col1,col2])
# # hdu = fits.PrimaryHDU(n)
# hdu.writeto('test1.fits',overwrite=True)

# hdul = fits.open('test1.fits')
# # hdul.info()
# print(hdul[0].data)
# print(hdul[1].data)
# exit()
# print(np.shape(hdul[0].data))
# hdul[0].data = hdul[0].data[0:100]
# print(np.shape(hdul[0].data))
# hdul.writeto('test2.fits',overwrite=True)

# hdul2 = fits.open('test2.fits')
# hdul2.info()
# print(np.shape(hdul[0].data))
# print(hdul[0].data)

# Demo code for ctod changes
hdul = fits.open('/Users/mpeel/Desktop/NOMINAL60A-130503-2106-000.ctod')
hdul.info()
print('-.-')
t = Table(hdul[1].data)
print(t)
print(np.result_type(hdul[1].data))
new_columns = []
mask = np.zeros(len(t[0][0]))
mask[0:100] = 1
for colname in t.columns:
	col = t[colname]
	print(str(col.shape) + ' - ' + str(len(col.shape)))
	if len(col.shape) == 3:
		if col.shape[1] > 1000:
			col = col[:, mask==1, :]
	elif len(col.shape) == 5:
		if col.shape[1] > 1000:
			col = col[:, mask==1, :, :, :]
	elif len(col.shape) > 1:
		if col.shape[1] > 1000:
			col = col[:,mask==1]
	print(' -> ' + str(col.shape))
	new_columns.append(col)
new_t = Table(new_columns, names=t.columns) 
print('.-.')
print(new_t)
hdul[1].data = np.array(new_t)
hdul.writeto('/Users/mpeel/Desktop/test_ctod.fits',overwrite=True)

hdul2 = fits.open('/Users/mpeel/Desktop/test_ctod.fits')
hdul2.info()
print(np.shape(hdul2[0].data))
# print(hdul2[1].data)
