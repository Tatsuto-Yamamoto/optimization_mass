import scipy.optimize
import numpy as np
# 質量の実験値	
mex_top = 173000
mex_up = 2.2


def objective_function(a,b,c,d,f,h,g):
	return ((a-b) * np.exp((a-b)/3) * np.sinh(c/3 - c*d) + (a-b)*np.sinh(c*d) - \
		c * np.exp((a-b)/3) * np.cosh(c/3 - c*d) + c*np.cosh(c*d)) / ((a-b) * np.exp((a-b)/3) * np.sinh(c - c*d)+\
		(-(a-b))* np.sinh(2*c/3 - c*d) - c*np.exp((a-b)/3) * np.cosh(c - c*d)+\
		c*np.cosh(2*c/3 - c*d)) / (mex_up/mex_top)

def gradient(theta):
    '''勾配'''
    return 2 * (theta - 2)

theta_opt = scipy.optimize.fmin_bfgs(f=objective_function, x0=[5.0,-14.48,8.6,0.28,-280000.0,34.9,17.0])
print(theta_opt)