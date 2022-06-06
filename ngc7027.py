# Code to estimate the flux density of NGC 7027 over time.
# Mike Peel    8 Feb 2017    Started
# Mike Peel    3 Mar 2017    Tidy labels

import matplotlib.pyplot as plt
import numpy as np

def calcflux_percent(t, nu, amp, ref_t, ref_nu, change_t, index):
	return (amp[0] * (1.0+change_t[0] * (t - ref_t)) ) * (nu/ref_nu)**index

def calcflux_subtract(t, nu, amp, ref_t, ref_nu, change_t, index):
	return (amp + change_t * (t - ref_t)) * (nu/ref_nu)**index

# Configuration
# Years
startyear = 2000.0
endyear = 2020.0
giveyear = 2017.0

# Hafez et al. (2008)
hafez08_flux = [5.39, 0.04]
hafez08_ref_t = 2003.0
hafez08_ref_nu = 33.0
hafez08_change_t = [-0.17/100.0, 0.03/100.0]
hafez08_index = -0.119

zijlstra08_flux = [5.42971, 5.42971*0.05]
zijlstra08_ref_t = 2000.0
zijlstra08_ref_nu = 30.199517
zijlstra08_change_t = -7.787/1000.0
zijlstra08_index = hafez08_index

perley13_flux = [5.480, 0.013]
perley13_ref_t = 2000.0
perley13_ref_nu = 22.460
perley13_change_t = -7.0/1000.0
perley13_index = hafez08_index

perley13_2_flux = [5.091, 0.042]
perley13_2_ref_t = 2000.0
perley13_2_ref_nu = 43.340
perley13_2_change_t = -5.0/1000.0
perley13_2_index = hafez08_index

fig1 = plt.figure(1)
ax1 = fig1.add_subplot(111)
dates = np.linspace(startyear,endyear, (endyear-startyear)+1)

# Zijlstra
maxval = calcflux_subtract(dates, 30.0, zijlstra08_flux[0]+zijlstra08_flux[1], zijlstra08_ref_t, zijlstra08_ref_nu, zijlstra08_change_t, zijlstra08_index)
minval = calcflux_subtract(dates, 30.0, zijlstra08_flux[0]-zijlstra08_flux[1], zijlstra08_ref_t, zijlstra08_ref_nu, zijlstra08_change_t, zijlstra08_index)
centerval = calcflux_subtract(dates, 30.0, zijlstra08_flux[0], zijlstra08_ref_t, zijlstra08_ref_nu, zijlstra08_change_t, zijlstra08_index)
ax1.fill_between(dates, maxval, minval, color='m',alpha=0.5)
plt.plot(dates, centerval,'m',label='Zijlstra et al. (2008)')
print 'Zijlstra model:'
print str(centerval[dates == giveyear]) + ' + ' + str(maxval[dates == giveyear] - centerval[dates == giveyear]) + ' - ' + str(abs(minval[dates == giveyear] - centerval[dates == giveyear]))

# Hafez
maxval = calcflux_percent(dates, 30.0, [hafez08_flux[0]+hafez08_flux[1]], hafez08_ref_t, hafez08_ref_nu, hafez08_change_t, hafez08_index)
minval = calcflux_percent(dates, 30.0, [hafez08_flux[0]-hafez08_flux[1]], hafez08_ref_t, hafez08_ref_nu, hafez08_change_t, hafez08_index)
centerval = calcflux_percent(dates, 30.0, hafez08_flux, hafez08_ref_t, hafez08_ref_nu, hafez08_change_t, hafez08_index)
ax1.fill_between(dates, maxval, minval, color='blue',alpha=0.5)
plt.plot(dates, centerval,'black', label='Hafez et al. (2008)')
print 'Hafez model:'
print str(centerval[dates == giveyear]) + ' + ' + str(maxval[dates == giveyear] - centerval[dates == giveyear]) + ' - ' + str(abs(minval[dates == giveyear] - centerval[dates == giveyear]))

# Perley
maxval = calcflux_subtract(dates, 30.0, perley13_flux[0]+perley13_flux[1], perley13_ref_t, perley13_ref_nu, perley13_change_t, perley13_index)
minval = calcflux_subtract(dates, 30.0, perley13_flux[0]-perley13_flux[1], perley13_ref_t, perley13_ref_nu, perley13_change_t, perley13_index)
centerval = calcflux_subtract(dates, 30.0, perley13_flux[0], perley13_ref_t, perley13_ref_nu, perley13_change_t, perley13_index)
ax1.fill_between(dates, maxval, minval, color='g',alpha=0.5, label='Perley et al. (2013) [22GHz]')
plt.plot(dates, centerval,'g')
print 'Perley model (extrapolate from 22.46GHz):'
print str(centerval[dates == giveyear]) + ' + ' + str(maxval[dates == giveyear] - centerval[dates == giveyear]) + ' - ' + str(abs(minval[dates == giveyear] - centerval[dates == giveyear]))

# Perley
maxval = calcflux_subtract(dates, 30.0, perley13_2_flux[0]+perley13_2_flux[1], perley13_2_ref_t, perley13_2_ref_nu, perley13_2_change_t, perley13_2_index)
minval = calcflux_subtract(dates, 30.0, perley13_2_flux[0]-perley13_2_flux[1], perley13_2_ref_t, perley13_2_ref_nu, perley13_2_change_t, perley13_2_index)
centerval = calcflux_subtract(dates, 30.0, perley13_2_flux[0], perley13_2_ref_t, perley13_2_ref_nu, perley13_2_change_t, perley13_2_index)
ax1.fill_between(dates, maxval, minval, color='lightgreen',alpha=0.5, label='Perley et al. (2013) [43GHz]')
plt.plot(dates, centerval,'lightgreen')
print 'Perley model (extrapolate from 43.34GHz):'
print str(centerval[dates == giveyear]) + ' + ' + str(maxval[dates == giveyear] - centerval[dates == giveyear]) + ' - ' + str(abs(minval[dates == giveyear] - centerval[dates == giveyear]))

plt.title('Flux density of NGC 7027 over time')
plt.xlabel('Time [year]')
plt.ylabel('Flux density [Jy]')
l = plt.legend(prop={'size':8})
l.set_zorder(20)
plt.savefig('ngc7027.png')