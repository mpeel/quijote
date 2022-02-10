import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#CALCULATING FIT TO THE DATA
def power(x,a,b):
	print(a)
	print(b)
	print(x)
	return a*x**(b)
#CURVE FIT TO THE DATA WITH POWER LAW
flux_normal = [59.17924124, 62.55201676, 30.77045018, 33.28932108, 67.63237391, 77.90297961]
frequency = [17., 19., 11., 13., 17., 19.]
print(flux_normal)
print(frequency)
#bounds=((0,0),(0.36,1.78))
k,c=curve_fit(power,frequency,flux_normal)
print(k[0],k[1])
x=np.arange(10,20,0.1)
y=k[0]*x**k[1]
plt.ion()
plt.plot(x,y)
plt.plot(frequency,flux_normal,'.')
plt.show()
input('test')