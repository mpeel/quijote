import numpy as np
import matplotlib.pyplot as plt

tatm,count,sigma=np.loadtxt('/Users/mpeel/Desktop/Datos30GHz0.txt',unpack=True)
weightsum = np.zeros(len(sigma))
weightsum[0] = 1.0/sigma[0]**2
for i in range(1,len(sigma)):
	weightsum[i] = weightsum[i-1] + np.power(sigma[i], -2)

plt.plot(sigma,label='sigma')
plt.plot(np.sqrt(1.0/weightsum),label='weighted uncertainty')
plt.legend()
plt.yscale('log')
plt.savefig('atmcheck.pdf')
