import healpy as hp
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.coordinates import get_body, EarthLocation, SkyCoord
from astropy.time import Time
import astropy.units as u
# Configuration
# input_file = '/Users/mpeel/Desktop/NOMINAL60A-130503-2106-000.ctod'

do_rotated = False
do_rotation = True
use_ecliptic = True
basedir = '/Users/mpeel/Desktop/'
if do_rotated == True:
	input_file = basedir+'JUPITERH3B-131010-0447-JC.ctod'
	input_file = basedir+'JUPITERC-150704-1407-JC.ctod'
	input_file = basedir+'JUPITERH3A-180626-2007-JC.ctod'
	pos = (0,0)
else:
	# input_file = basedir+'JUPITERA-150706-1359.ctod'
	# t = Time('2015-07-06T13:59:00', format='isot', scale='utc')

	input_file = basedir+'JUPITERD-150705-1403.ctod'
	t = Time('2015-07-05T14:03:00', format='isot', scale='utc')



	# input_file = basedir+'JUPITERH3B-131010-0447.ctod'
	# t = Time('2013-10-10T04:47:00', format='isot', scale='utc')

	# input_file = basedir+'JUPITERC-150704-1407.ctod'
	# t = Time('2015-07-01T15:07:00', format='isot', scale='utc')

	# input_file = basedir+'JUPITERH3A-180626-2007.ctod'
	# t = Time('2018-06-26T20:07:00', format='isot', scale='utc')

	# Jupiter position
	loc=EarthLocation(-16.5085,28.3,2390)
	jupiter=get_body("Jupiter",t,loc)
	print(jupiter)
	c_icrs = SkyCoord(ra=jupiter.ra.deg*u.degree, dec=jupiter.dec.deg*u.degree, frame='icrs')
	if use_ecliptic:
		print(c_icrs)
		# ecliptic = c_icrs.GeocentricMeanEcliptic()
		print(c_icrs.geocentrictrueecliptic)
		l_j = c_icrs.geocentrictrueecliptic.lon.rad
		b_j = c_icrs.geocentrictrueecliptic.lat.rad
		pos = (c_icrs.geocentrictrueecliptic.lon.deg,c_icrs.geocentrictrueecliptic.lat.deg)
		# exit()
	else:
		print(c_icrs.galactic)
		pos = (c_icrs.galactic.l.deg,c_icrs.galactic.b.deg)
		l_j = c_icrs.galactic.l.rad
		b_j = c_icrs.galactic.b.rad
	print(l_j)
	print(b_j)

nside = 512

# Set up the map
npix=12*nside*nside
outmap = np.zeros(npix)
outmap[:] = hp.UNSEEN

# Read in the data
hdul = fits.open(input_file)

# Select what we want to map
gl = np.asarray(hdul[1].data.field('GL')[0])
gb = np.asarray(hdul[1].data.field('GB')[0])
data = np.asarray(hdul[1].data.field('DATA')[0])

# Select the horn and channel
hornnum = 3
channum = 23
gl = gl[:,hornnum-1]
gb = gb[:,hornnum-1]
data = data[:,channum-1]
# data = gl.copy()

# Change from Galactic to Ecliptic coordinates
if use_ecliptic:
	print(gl)
	positions = SkyCoord(gl*u.degree, gb*u.degree, frame='galactic')
	gl = positions.geocentrictrueecliptic.lon.deg
	gb = positions.geocentrictrueecliptic.lat.deg
	print(gl)


if do_rotation:
	pos = (0,0)
	if use_ecliptic:
		input_file = input_file.replace('.ctod','_newrot_ecliptic.ctod')
	else:
		input_file = input_file.replace('.ctod','_newrot.ctod')

	#Turn each l and b into a cartesian vector
	phi=gl*np.pi/180.0
	theta=(90.0-gb)*np.pi/180.0
	print(phi)
	print(theta)
	
	#Change to cartesian coordinates.
	x=np.sin(theta)*np.cos(phi)
	y=np.sin(theta)*np.sin(phi)
	z=np.cos(theta)
	data_cart=[[x],[y],[z]]   
	data_cart=np.array(data_cart)        

	#FIRST ROTATION
	
	#create cos and sin
	cos, sin = np.cos(l_j), np.sin(l_j)

	#Create rotation matrix
	Rz=[[cos,sin,0],
		[-sin,cos,0],
		[0,0,1]]
	Rz=np.array(Rz)

	#New cos and sin
	cos2, sin2 = np.cos(b_j), np.sin(b_j)
	Ry=[[cos2,0,sin2],
		[0,1,0],
		[-sin2,0,cos2]]
	Ry=np.array(Ry)
			
	Rz=np.transpose(Rz)
	Ry=np.transpose(Ry)
	answer=Rz.dot(Ry)
	answer=np.transpose(answer)
	
	C_new1=[]  
	for k in range(len(gl)):
		temp4=answer.dot(data_cart[:,0,k])
		C_new1.append(temp4)
	C_new1=np.array(C_new1) 

	#Now convert C_new1 back to galactic and save to a new fits file

	#transform back to spherical coordinates
	phi_new1=np.arctan(C_new1[:,1]/C_new1[:,0])
	# theta_new1=np.arccos(C_new1[:,2]/(np.sqrt(C_new1[:,0]**2+C_new1[:,1]**2+C_new1[:,2]**2)))
	theta_new1 = np.arccos(C_new1[:,2])

	gl = phi_new1*180.0/np.pi
	gb = 90.0-theta_new1*180.0/np.pi

# Convert to Healpix pixels
healpix_pixel = np.asarray(hp.ang2pix(nside, (np.pi/2)-gb*np.pi/180.0, gl*np.pi/180.0))

# Bin the data into the map
for pixnum in np.unique(healpix_pixel[healpix_pixel > 0]):
	# print(pixnum)
	if pixnum < npix:
		outmap[pixnum] = np.median(data[np.where(healpix_pixel==pixnum)])


maxpix = hp.pix2ang(nside, np.argmax(outmap),lonlat=True)
print(maxpix)
print('Diff in lon:' + str(maxpix[0]-pos[0]))
print('Diff in lat:' + str(maxpix[1]-pos[1]))

# Write out the map and plot it
hp.write_map(input_file.replace('.ctod','_map.fits'),outmap,overwrite=True)
hp.mollview(outmap)
plt.savefig(input_file.replace('.ctod','_map.png'))
plt.clf()
hp.gnomview(outmap,reso=4.0,cmap='jet',rot=pos)
plt.savefig(input_file.replace('.ctod','_mapzoom.png'))
