import astropy as ap
import numpy as np
import astropy.units as u
from astropy.time import Time, TimezoneInfo
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun, get_moon, get_body
import matplotlib.pyplot as plt
import datetime
from matplotlib.dates import (YEARLY, DateFormatter,
                              rrulewrapper, RRuleLocator, drange)
import datetime

import gspread
# from oauth2client.service_account import ServiceAccountCredentials

gc = gspread.service_account(filename='astrocode-3674ffba7119.json')

sheet = gc.open("TFGI observations").sheet1
# print(sheet.get('A1'))

telescope = EarthLocation(lat=28.300224*u.deg, lon=-16.510113*u.deg, height=2390*u.m)

today = datetime.date.today()
todaystr = today.strftime("%Y-%m-%d 00:00:00")
todayTime = Time(todaystr)
TimeArr = todayTime + np.arange(24) * datetime.timedelta(hours=1)

source_names = ['Tau A', 'Cas A', 'Cyg A']
source_pos = []
for i in range(0,len(source_names)):
	source_pos.append(SkyCoord.from_name(source_names[i]))
	print(source_names[i])
	print(source_pos[i].to_string(style='hmsdms'))

print(TimeArr[0].value[0:10])

vals = np.zeros((24+4,4+len(source_names)))
for i in range(0,24):
	time = TimeArr[i]
	# print(time.value[11:16])

	# Sun
	sun = get_sun(time)
	altaz = sun.transform_to(AltAz(obstime=time,location=telescope))
	# print(altaz.alt.deg)
	vals[i, 0] = np.round(altaz.alt.deg)

	# Moon
	moon = get_moon(time)
	altaz = moon.transform_to(AltAz(obstime=time,location=telescope))
	# print(altaz.alt.deg)
	vals[i, 1] = np.round(altaz.alt.deg)

	# Jupiter
	jupiter = get_body('jupiter',time, telescope)
	altaz = jupiter.transform_to(AltAz(obstime=time,location=telescope))
	# print(altaz.alt.deg)
	vals[i, 2] = np.round(altaz.alt.deg)

	# Venus
	venus = get_body('venus',time, telescope)
	altaz = venus.transform_to(AltAz(obstime=time,location=telescope))
	# print(altaz.alt.deg)
	vals[i, 3] = np.round(altaz.alt.deg)

	for j in range(0,len(source_names)):
		# print(source_names[j])
		altaz = source_pos[j].transform_to(AltAz(obstime=time,location=telescope))
		# print(altaz.alt.deg)
		vals[i, j+4] = np.round(altaz.alt.deg)

sampling = 1000
time_delta = (TimeArr[-1].mjd - TimeArr[0].mjd)/sampling
times_mjd = np.arange(TimeArr[0].mjd,TimeArr[-1].mjd,time_delta)
times = Time(times_mjd,format='mjd')

for j in range(0,len(source_names)):
	el40 = 0
	el75 = 0
	elm75 = 0
	elm40 = 0
	altaz = source_pos[j].transform_to(AltAz(obstime=times,location=telescope))
	for k in range(0,sampling):
		# print(altaz.alt[j].deg)
		if altaz.alt[k].deg > 40 and el40 == 0:
			vals[24,j+4] = float(times[k].datetime.strftime("%H.%M"))
			el40 = 1
			elm40 = 0
		if altaz.alt[k].deg > 75 and el75 == 0:
			vals[25,j+4] = float(times[k].datetime.strftime("%H.%M"))
			el75 = 1
			elm75 = 0
		if altaz.alt[k].deg < 40:
			el40 = 0
		if altaz.alt[k].deg < 75 and elm75 == 0:
			vals[26,j+4] = float(times[k].datetime.strftime("%H.%M"))
			elm75 = 1
		if altaz.alt[k].deg < 40 and elm40 == 0:
			vals[27,j+4] = float(times[k].datetime.strftime("%H.%M"))
			elm40 = 1

print(vals)
sheet.update('B3', vals.tolist())
sheet.update_cell(1, 2, TimeArr[0].value[0:10])
print(sheet.cell(1, 2).value)
