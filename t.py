import numpy as np
from scipy.optimize import minimize
import sympy as sym

init = np.array([5.0,-14.48,8.6,0.28,-280000.0,34.9,17.0],dtype=np.float128)


def ups(x): #理論式と実験値の比を返す		
	#sym.setter(over="Ignore")
	a = x[0]
	b = x[1]
	c = x[2]
	d = x[3]
	f = x[4]
	h = x[5]
	g = x[6]
	r = ((a-b) * sym.exp((a-b)/3) * sym.sinh(c/3 - c*d) + (a-b)*sym.sinh(c*d) - \
		c * sym.exp((a-b)/3) * sym.cosh(c/3 - c*d) + c*sym.cosh(c*d)) / ((a-b) * sym.exp((a-b)/3) * sym.sinh(c - c*d)+\
		(-(a-b))* sym.sinh(2*c/3 - c*d) - c*sym.exp((a-b)/3) * sym.cosh(c - c*d)+\
		c*sym.cosh(2*c/3 - c*d)) / (mex_up/mex_top)
	return r
def charms(x):
	a = x[0]
	b = x[1]
	c = x[2]
	d = x[3]
	f = x[4]
	h = x[5]
	g = x[6]

	r  = ((a-b)*sym.exp((a-b)/3)*sym.sinh(2*c/3 - c*d) + (-(a-b))*sym.sinh(c/3 - c*d) -
    c*sym.exp((a-b)/3)*sym.cosh(2* c/3 - c*d) +
    c*sym.cosh(c/3 - c*d))/((a-b)*sym.exp((a-b)/3) *sym.sinh(c - c*d) + (-(a-b))*
     sym.sinh(2*c/3 - c*d) - c*sym.exp((a-b)/3)* sym.cosh(c - c*d) +
    c*sym.cosh(2*c/3 - c*d))/(1275/173000)
	return r
def downs():
	# LaTeXで数式を表示
	sym.init_printing()
	sym.var('a')
	sym.var('b')
	sym.var('c')
	sym.var('d')
	sym.var('f')
	sym.var('h')
	sym.var('g')
	r = ((a - b)**2 - c**2)/((a - f)**2 -
    c**2) * sym.sqrt((f*(-sym.exp(-2*b/3) + 1))/(b*(-sym.exp(-2*f/3) +
       1))) * ((a - f)*sym.exp((a - f)/3)*sym.sinh(c/3 - c*d) + (a - f)*
      sym.sinh(c*d) - c*sym.exp((a - f)/3)*sym.cosh(c/3 - c*d) +
     c*sym.cosh(-c*d))/((a - b)*sym.exp((a - b)/3)* sym.sinh(c - c*d) + (-a + b)*
      sym.sinh(2*c/3 - c*d) - c*sym.exp((a - b)/3)* sym.cosh(c - c*d) +
     c*sym.cosh(2*c/3 - c*d))/(4.7/173000)
	
	#display(r)
	return r
downs()
def stranges(x):
	a = x[0]
	b = x[1]
	c = x[2]
	d = x[3]
	f = x[4]
	h = x[5]
	g = x[6]

	r = ((a - b)**2 - c**2)/((a - f)**2 -
    c**2)*sym.sqrt((f*(-sym.exp(-2*b/3) + 1))/(b*(-sym.exp(-2*f/3) +
       1))) * ((a - f)*sym.exp((a - f)/3)*sym.sinh(2 * c/3 - c*d) + (-a + f)*
      sym.sinh(c/3 - c*d) - c*sym.exp((a - f)/3)*sym.cosh(2*c/3 - c*d) +
     c*sym.cosh(c/3 - c*d))/((a - b)*
      sym.exp((a - b)/3)* sym.sinh(c - c*d) + (-a + b)*sym.sinh(2*c/3 - c*d) -
     c*sym.exp((a - b)/3)* sym.cosh(c - c*d) + c*sym.cosh(2*c/3 - c*d))/(95/
    173000)
	return r
def botoms(x):
	a = x[0]
	b = x[1]
	c = x[2]
	d = x[3]
	f = x[4]
	h = x[5]
	g = x[6]

	r = ((a - b)**2 - c**2)/((a - f)**2 -
    c**2) * sym.sqrt((f*(-sym.exp(-2*b/3) + 1))/(b*(-sym.exp(-2*f/3) +
       1))) * ((a - f)*sym.exp((a - f)/3)*sym.sinh(c - c*d) + (-a + f)*
      sym.sinh(2 * c/3 - c*d) - c*sym.exp((a - f)/3)*sym.cosh(c - c*d) +
     c*sym.cosh(2 * c/3 - c*d))/((a - b)*
      sym.exp((a - b)/3) * sym.sinh(c - c*d) + (-a + b)*sym.sinh(2*c/3 - c*d) -
     c*sym.exp((a - b)/3) * sym.cosh(c - c*d) + c*sym.cosh(2*c/3 - c*d))/(41.8/
    1730)

	return r
def electrons(x):
	a = x[0]
	b = x[1]
	c = x[2]
	d = x[3]
	f = x[4]
	h = x[5]
	g = x[6]

	r = ((a - b)**2 - c**2)/((g - h)**2 - c**2) * \
		sym.sqrt((g*(sym.exp(2*a/3) - 1))/(a*(sym.exp(2*g/3) - 1)))*sym.sqrt((h*(-sym.exp(-2*b/3) + 1))/(b*(-sym.exp(-2*h/3) + 1))) * ((g - h)*
		sym.exp((g - h)/3)*sym.sinh(c/3 - c*d) + (g - h)*sym.sinh(c*d) -
		c*sym.exp((g - h)/3)*sym.cosh(c/3 - c*d) +
		c*sym.cosh(-c*d))/((a - b)*sym.exp((a - b)/3)* sym.sinh(c - c*d) + (-a + b)*
		sym.sinh(2*c/3 - c*d) - c*sym.exp((a - b)/3)* sym.cosh(c - c*d) +
	c*sym.cosh(2*c/3 - c*d))/(0.510999/173000)

	return r
def muons(x):
	a = x[0]
	b = x[1]
	c = x[2]
	d = x[3]
	f = x[4]
	h = x[5]
	g = x[6]

	r = ((a - b)**2 - c**2)/((g - h)**2 - c**2)*\
 		sym.sqrt((g*(sym.exp(2*a/3) - 1))/(a*(sym.exp(2*g/3) - 1)))*\
 		sym.sqrt((h*(-sym.exp(-2*b/3) + 1))/(b*(-sym.exp(-2*h/3) + 1))) * ((g - h)*
      	sym.exp((g - h)/3)*sym.sinh(2*c/3 - c*d) + (-g + h)*sym.sinh(c/3 - c*d) -
     	c*sym.exp((g - h)/3)*sym.cosh(2*c/3 - c*d) +
     	c*sym.cosh(c/3 - c*d))/((a - b)*
      	sym.exp((a - b)/3) * sym.sinh(c - c*d) + (-a + b)*sym.sinh(2*c/3 - c*d) - \
     	c*sym.exp((a - b)/3) * sym.cosh(c - c*d) + c*sym.cosh(2*c/3 - c*d))/(105.65/
    	173000)

	return r
def tauons(x):
	a = x[0]
	b = x[1]
	c = x[2]
	d = x[3]
	f = x[4]
	h = x[5]
	g = x[6]

	r = ((a - b)**2 - c**2)/((g - h)**2 - c**2)*\
	 	sym.sqrt((g*(sym.exp(2*a/3) - 1))/(a*(sym.exp(2*g/3) - 1)))*\
 		sym.sqrt((h*(-sym.exp(-2*b/3) + 1))/(b*(-sym.exp(-2*h/3) + 1))) * ((g - h)*
	 	sym.exp((g - h)/3)*sym.sinh(c - c*d) + (-g + h)*sym.sinh(2*c/3 - c*d) -
	 	c*sym.exp((g - h)/3)*sym.cosh(c - c*d) +
	 	c*sym.cosh(2*c/3 - c*d))/((a - b)*
	 	sym.exp((a - b)/3) *sym.sinh(c - c*d) + (-a + b)*sym.sinh(2*c/3 - c*d) -
	 	c*sym.exp((a - b)/3) *sym.cosh(c - c*d) + c*sym.cosh(2*c/3 - c*d))/(1776.86/
	 	173000)
