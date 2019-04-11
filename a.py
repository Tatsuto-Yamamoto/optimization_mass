import numpy as np
import sympy as sym

# 質量の実験値	
mex_top = 173000
mex_up = 2.2

def ups(x): #理論式と実験値の比を返す		
	a = x[0]
	b = x[1]
	c = x[2]
	d = x[3]
	f = x[4]
	h = x[5]
	g = x[6]
	#print(c/3 - c*d)
	p1 = (a-b) * np.exp((a-b)/3) * np.sinh(c/3 - c*d) + (a-b)*np.sinh(c*d)
	print("p1=", p1)
	p2 = c * np.exp((a-b)/3) * np.cosh(c/3 - c*d) - c*np.cosh(c*d)
	p3 = (a-b) * np.exp((a-b)/3) * np.sinh(c - c*d)
	p4 = (-(a-b))* np.sinh(2*c/3 - c*d) - c*np.exp((a-b)/3) * np.cosh(c - c*d)
	p5 = c*np.cosh(2*c/3 - c*d)
	p6 = (mex_up/mex_top)
	r1 = (p1-p2)/(p3+p4+p5)/p6
	r2 = ((a-b) * np.exp((a-b)/3) * np.sinh(c/3 - c*d) + (a-b)*np.sinh(c*d) - \
		c * np.exp((a-b)/3) * np.cosh(c/3 - c*d) + c*np.cosh(c*d)) / ((a-b) * np.exp((a-b)/3) * np.sinh(c - c*d)+\
		(-(a-b))* np.sinh(2*c/3 - c*d) - c*np.exp((a-b)/3) * np.cosh(c - c*d)+\
		c*np.cosh(2*c/3 - c*d)) / (mex_up/mex_top)
	print(r1,r2)
	return r2
x1 = np.array([5,-14.48,8.6,0.28,-280000,34.9,17])
ups(x1)