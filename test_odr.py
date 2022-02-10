import numpy as np
from scipy import odr

def linfit3(param, x):
	return param[0]*x+param[1]


x = np.arange(0,10)+np.random.normal(scale=0.1, size=10)
y = np.arange(0,10)
sigma = np.ones(len(x))
cov = np.zeros((len(x),len(x)))
for i in range(0,len(x)):
	cov[i][i] = 1.0

params = [0,0]
odr_model = odr.Model(linfit3)
dataset = odr.RealData(x, y, sx=sigma, covy=cov)
odr_run = odr.ODR(dataset, odr_model, beta0=params)
out = odr_run.run()
param_est = out.beta
sigma_param_est = out.sd_beta
print(param_est)
print(sigma_param_est)