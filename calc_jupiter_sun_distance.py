# Plot the distance from Earth vs. the distance from the Sun for Jupiter
from astropy.time import Time
from astropy.coordinates import solar_system_ephemeris, EarthLocation
from astropy.coordinates import get_body_barycentric, get_body, get_sun
import astropy.units as u
import numpy as np
import astropy as ap
import matplotlib.pyplot as plt

telescope = EarthLocation(lat=28.300467*u.deg, lon=-16.510288*u.deg, height=2390*u.m)

# 1 Jan 2012
mjds = 55927 + np.arange(0,6*365)
times = ap.time.Time(mjds, format='mjd', scale='utc')
jupiter = get_body('jupiter',times, telescope)
sun = get_body('sun',times, telescope)
distance = jupiter.separation(sun)
print(distance.deg)
print(jupiter.distance)
plt.plot(distance.deg, jupiter.distance, '-')
plt.title('Jupiter distance between 1 Jan 2012 and 1 Jan 2018')
plt.xlabel('Distance from Sun (deg)')
plt.ylabel('Distance from Earth (AU)')
plt.savefig('calc_jupiter_sun_distance.png')