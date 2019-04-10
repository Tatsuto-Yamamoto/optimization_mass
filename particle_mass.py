# 与えたパラメータに対して、質量とその勾配を返す
import numpy as np
import sympy as sym

# 以下の値が、"いい感じ" な値らしい
# x = np.array([5,-14.48,8.6,0.28,-280000,34.9,17])

# 質量の実験値	
mex_top = 173000
mex_up = 2.2



def ups(x): #理論式と実験値の比を返す		
	#np.seterr(over="ignore")
	#x = np.array(x, dtype = np.float128)
	a = x[0]
	b = x[1]
	c = x[2]
	d = x[3]
	f = x[4]
	h = x[5]
	g = x[6]
	r = ((a-b) * np.exp((a-b)/3) * np.sinh(c/3 - c*d) + (a-b)*np.sinh(c*d) - \
		c * np.exp((a-b)/3) * np.cosh(c/3 - c*d) + c*np.cosh(c*d)) / ((a-b) * np.exp((a-b)/3) * np.sinh(c - c*d)+\
		(-(a-b))* np.sinh(2*c/3 - c*d) - c*np.exp((a-b)/3) * np.cosh(c - c*d)+\
		c*np.cosh(2*c/3 - c*d)) / (mex_up/mex_top)
	return r
def charms(x):
	a = x[0]
	b = x[1]
	c = x[2]
	d = x[3]
	f = x[4]
	h = x[5]
	g = x[6]
	r  = ((a-b)*np.exp((a-b)/3)*np.sinh(2*c/3 - c*d) + (-(a-b))*np.sinh(c/3 - c*d) -
    c*np.exp((a-b)/3)*np.cosh(2* c/3 - c*d) +
    c*np.cosh(c/3 - c*d))/((a-b)*np.exp((a-b)/3) *np.sinh(c - c*d) + (-(a-b))*
     np.sinh(2*c/3 - c*d) - c*np.exp((a-b)/3)* np.cosh(c - c*d) +
    c*np.cosh(2*c/3 - c*d))/(1275/173000)
	return r
def downs(x):

	a = x[0]
	b = x[1]
	c = x[2]
	d = x[3]
	f = x[4]
	h = x[5]
	g = x[6]
	# p1 = ((a - b)**2.0 - c**2.0)/((a - f)**2.0 -\
 #    c**2.0) * np.sqrt((f*(-np.exp(-2.0*b/3.0) + 1.0))/(b*(-np.exp(-2.0*f/3.0) +\
 #       1.0)))
	# p2 = ((a - f)*np.exp((a - f)/3.0)*np.sinh(c/3.0 - c*d) + (a - f) * 
	# 		np.sinh(c*d) - c*np.exp((a - f)/3.0)*np.cosh(c/3.0 - c*d) +\
	# 		c*np.cosh(-c*d))
	# p3 = ((a - b)*np.exp((a - b)/3.0)* np.sinh(c - c*d) + (-a + b)*\
	# 	np.sinh(2.0*c/3.0 - c*d) - c*np.exp((a - b)/3.0)* np.cosh(c - c*d) +\
	# 		c*np.cosh(2.0*c/3.0 - c*d))/(4.7/173000.0)
	# o = p1*p2/p3
	# print("p1 =", p1)
	# print("p2 = ", p2)
	# print("p3 =", p3)
	r = ((a - b)**2.0 - c**2.0)/((a - f)**2.0 -\
		c**2.0) * np.sqrt((f*(-np.exp(-2.0*b/3.0) + 1.0))/(b*(-np.exp(-2.0*f/3.0) +\
		1.0))) * ((a - f)*np.exp((a - f)/3.0)*np.sinh(c/3.0 - c*d) + (a - f)*\
		np.sinh(c*d) - c*np.exp((a - f)/3.0)*np.cosh(c/3.0 - c*d) +\
		c*np.cosh(-c*d))/((a - b)*np.exp((a - b)/3.0)* np.sinh(c - c*d) + (-a + b)*\
		np.sinh(2.0*c/3.0 - c*d) - c*np.exp((a - b)/3.0)* np.cosh(c - c*d) +\
		c*np.cosh(2.0*c/3.0 - c*d))/(4.7/173000.0)
	return r
def stranges(x):
	a = x[0]
	b = x[1]
	c = x[2]
	d = x[3]
	f = x[4]
	h = x[5]
	g = x[6]

	r = ((a - b)**2 - c**2)/((a - f)**2 -
    c**2)*np.sqrt((f*(-np.exp(-2*b/3) + 1))/(b*(-np.exp(-2*f/3) +
       1))) * ((a - f)*np.exp((a - f)/3)*np.sinh(2 * c/3 - c*d) + (-a + f)*
      np.sinh(c/3 - c*d) - c*np.exp((a - f)/3)*np.cosh(2*c/3 - c*d) +
     c*np.cosh(c/3 - c*d))/((a - b)*
      np.exp((a - b)/3)* np.sinh(c - c*d) + (-a + b)*np.sinh(2*c/3 - c*d) -
     c*np.exp((a - b)/3)* np.cosh(c - c*d) + c*np.cosh(2*c/3 - c*d))/(95/
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
    c**2) * np.sqrt((f*(-np.exp(-2*b/3) + 1))/(b*(-np.exp(-2*f/3) +
       1))) * ((a - f)*np.exp((a - f)/3)*np.sinh(c - c*d) + (-a + f)*
      np.sinh(2 * c/3 - c*d) - c*np.exp((a - f)/3)*np.cosh(c - c*d) +
     c*np.cosh(2 * c/3 - c*d))/((a - b)*
      np.exp((a - b)/3) * np.sinh(c - c*d) + (-a + b)*np.sinh(2*c/3 - c*d) -
     c*np.exp((a - b)/3) * np.cosh(c - c*d) + c*np.cosh(2*c/3 - c*d))/(41.8/
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
		np.sqrt((g*(np.exp(2*a/3) - 1))/(a*(np.exp(2*g/3) - 1)))*np.sqrt((h*(-np.exp(-2*b/3) + 1))/(b*(-np.exp(-2*h/3) + 1))) * ((g - h)*
		np.exp((g - h)/3)*np.sinh(c/3 - c*d) + (g - h)*np.sinh(c*d) -
		c*np.exp((g - h)/3)*np.cosh(c/3 - c*d) +
		c*np.cosh(-c*d))/((a - b)*np.exp((a - b)/3)* np.sinh(c - c*d) + (-a + b)*
		np.sinh(2*c/3 - c*d) - c*np.exp((a - b)/3)* np.cosh(c - c*d) +
	c*np.cosh(2*c/3 - c*d))/(0.510999/173000)

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
 		np.sqrt((g*(np.exp(2*a/3) - 1))/(a*(np.exp(2*g/3) - 1)))*\
 		np.sqrt((h*(-np.exp(-2*b/3) + 1))/(b*(-np.exp(-2*h/3) + 1))) * ((g - h)*
      	np.exp((g - h)/3)*np.sinh(2*c/3 - c*d) + (-g + h)*np.sinh(c/3 - c*d) -
     	c*np.exp((g - h)/3)*np.cosh(2*c/3 - c*d) +
     	c*np.cosh(c/3 - c*d))/((a - b)*
      	np.exp((a - b)/3) * np.sinh(c - c*d) + (-a + b)*np.sinh(2*c/3 - c*d) - \
     	c*np.exp((a - b)/3) * np.cosh(c - c*d) + c*np.cosh(2*c/3 - c*d))/(105.65/
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
	 	np.sqrt((g*(np.exp(2*a/3) - 1))/(a*(np.exp(2*g/3) - 1)))*\
 		np.sqrt((h*(-np.exp(-2*b/3) + 1))/(b*(-np.exp(-2*h/3) + 1))) * ((g - h)*
	 	np.exp((g - h)/3)*np.sinh(c - c*d) + (-g + h)*np.sinh(2*c/3 - c*d) -
	 	c*np.exp((g - h)/3)*np.cosh(c - c*d) +
	 	c*np.cosh(2*c/3 - c*d))/((a - b)*
	 	np.exp((a - b)/3) *np.sinh(c - c*d) + (-a + b)*np.sinh(2*c/3 - c*d) -
	 	c*np.exp((a - b)/3) *np.cosh(c - c*d) + c*np.cosh(2*c/3 - c*d))/(1776.86/
	 	173000)

	return r


def dups(x): #upsの勾配
	a = sym.Symbol('a')
	b = sym.Symbol('b')
	c = sym.Symbol('c')
	d = sym.Symbol('d')
	f = sym.Symbol('f')
	h = sym.Symbol('h')
	g = sym.Symbol('g')

	a_now = x[0]
	b_now = x[1]
	c_now = x[2]
	d_now = x[3]
	f_now = x[4]
	h_now = x[5]
	g_now = x[6]

	# sympyを用いて、文字式として保持する
	mass = ((a-b) * sym.exp((a-b)/3) * sym.sinh(c/3 - c*d) + (a-b)*sym.sinh(c*d) - 
      c * sym.exp((a-b)/3) * sym.cosh(c/3 - c*d) + c*sym.cosh(c*d)) / ((a-b) * sym.exp((a-b)/3) * sym.sinh(c - c*d) 
		+ (-(a-b))* sym.sinh(2*c/3 - c*d) - c*sym.exp((a-b)/3) * sym.cosh(c - c*d)
		 + c*sym.cosh(2*c/3 - c*d)) / (mex_up/mex_top)

    # df* は * で偏微分した結果を値で返す
    # ここは、for で効率化できる・・・	
	dfa = sym.diff(mass,a).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])
	dfb = sym.diff(mass,b).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])
	dfc = sym.diff(mass,c).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])
	dfd = sym.diff(mass,d).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])
	dff = sym.diff(mass,f).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])
	dfh = sym.diff(mass,h).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])
	dfg = sym.diff(mass,g).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])

	grad = np.array([dfa,dfb,dfc,dfd,dff,dfh,dfg],dtype = "float64") #ここで、型を浮動小数点で統一する	
	return grad

def dcharms(x): #upsの勾配
	a = sym.Symbol('a')
	b = sym.Symbol('b')
	c = sym.Symbol('c')
	d = sym.Symbol('d')
	f = sym.Symbol('f')
	h = sym.Symbol('h')
	g = sym.Symbol('g')

	a_now = x[0]
	b_now = x[1]
	c_now = x[2]
	d_now = x[3]
	f_now = x[4]
	h_now = x[5]
	g_now = x[6]

	mass = ((a-b)*sym.exp((a-b)/3)*sym.sinh(2*c/3 - c*d) + (-(a-b))*sym.sinh(c/3 - c*d) -
    c*sym.exp((a-b)/3)*sym.cosh(2* c/3 - c*d) +
    c*sym.cosh(c/3 - c*d))/((a-b)*sym.exp((a-b)/3) *sym.sinh(c - c*d) + (-(a-b))*
     sym.sinh(2*c/3 - c*d) - c*sym.exp((a-b)/3)* sym.cosh(c - c*d) +
    c*sym.cosh(2*c/3 - c*d))/(1275/173000)

	dfa = sym.diff(mass,a).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])
	dfb = sym.diff(mass,b).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])
	dfc = sym.diff(mass,c).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])
	dfd = sym.diff(mass,d).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])
	dff = sym.diff(mass,f).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])
	dfh = sym.diff(mass,h).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])
	dfg = sym.diff(mass,g).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])

	grad = np.array([dfa,dfb,dfc,dfd,dff,dfh,dfg],dtype = "float64") #ここで、型を浮動小数点で統一する	
	return grad

def ddowns(x): #upsの勾配
	a = sym.Symbol('a')
	b = sym.Symbol('b')
	c = sym.Symbol('c')
	d = sym.Symbol('d')
	f = sym.Symbol('f')
	h = sym.Symbol('h')
	g = sym.Symbol('g')

	a_now = x[0]
	b_now = x[1]
	c_now = x[2]
	d_now = x[3]
	f_now = x[4]
	h_now = x[5]
	g_now = x[6]

	mass = ((a - b)**2 - c**2)/((a - f)**2 -
    c**2) * sym.sqrt((f*(-sym.exp(-2*b/3) + 1))/(b*(-sym.exp(-2*f/3) +
       1))) * ((a - f)*sym.exp((a - f)/3)*sym.sinh(c/3 - c*d) + (a - f)*
      sym.sinh(c*d) - c*sym.exp((a - f)/3)*sym.cosh(c/3 - c*d) +
     c*sym.cosh(-c*d))/((a - b)*sym.exp((a - b)/3)* sym.sinh(c - c*d) + (-a + b)*
      sym.sinh(2*c/3 - c*d) - c*sym.exp((a - b)/3)* sym.cosh(c - c*d) +
     c*sym.cosh(2*c/3 - c*d))/(4.7/173000)

	dfa = sym.diff(mass,a).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])
	dfb = sym.diff(mass,b).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])
	dfc = sym.diff(mass,c).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])
	dfd = sym.diff(mass,d).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])
	dff = sym.diff(mass,f).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])
	dfh = sym.diff(mass,h).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])
	dfg = sym.diff(mass,g).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])

	grad = np.array([dfa,dfb,dfc,dfd,dff,dfh,dfg],dtype = "float64") #ここで、型を浮動小数点で統一する	
	return grad
def dstranges(x): #upsの勾配
	a = sym.Symbol('a')
	b = sym.Symbol('b')
	c = sym.Symbol('c')
	d = sym.Symbol('d')
	f = sym.Symbol('f')
	h = sym.Symbol('h')
	g = sym.Symbol('g')

	a_now = x[0]
	b_now = x[1]
	c_now = x[2]
	d_now = x[3]
	f_now = x[4]
	h_now = x[5]
	g_now = x[6]

	mass = ((a - b)**2 - c**2)/((a - f)**2 -
    c**2)*sym.sqrt((f*(-sym.exp(-2*b/3) + 1))/(b*(-sym.exp(-2*f/3) +
       1))) * ((a - f)*sym.exp((a - f)/3)*sym.sinh(2 * c/3 - c*d) + (-a + f)*
      sym.sinh(c/3 - c*d) - c*sym.exp((a - f)/3)*sym.cosh(2*c/3 - c*d) +
     c*sym.cosh(c/3 - c*d))/((a - b)*
      sym.exp((a - b)/3)* sym.sinh(c - c*d) + (-a + b)*sym.sinh(2*c/3 - c*d) -
     c*sym.exp((a - b)/3)* sym.cosh(c - c*d) + c*sym.cosh(2*c/3 - c*d))/(95/
    173000)

	dfa = sym.diff(mass,a).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])
	dfb = sym.diff(mass,b).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])
	dfc = sym.diff(mass,c).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])
	dfd = sym.diff(mass,d).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])
	dff = sym.diff(mass,f).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])
	dfh = sym.diff(mass,h).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])
	dfg = sym.diff(mass,g).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])

	grad = np.array([dfa,dfb,dfc,dfd,dff,dfh,dfg],dtype = "float64") #ここで、型を浮動小数点で統一する	
	return grad
def dbotoms(x): #upsの勾配
	a = sym.Symbol('a')
	b = sym.Symbol('b')
	c = sym.Symbol('c')
	d = sym.Symbol('d')
	f = sym.Symbol('f')
	h = sym.Symbol('h')
	g = sym.Symbol('g')

	a_now = x[0]
	b_now = x[1]
	c_now = x[2]
	d_now = x[3]
	f_now = x[4]
	h_now = x[5]
	g_now = x[6]

	mass = ((a - b)**2 - c**2)/((a - f)**2 -
    c**2) * sym.sqrt((f*(-sym.exp(-2*b/3) + 1))/(b*(-sym.exp(-2*f/3) +
       1))) * ((a - f)*sym.exp((a - f)/3)*sym.sinh(c - c*d) + (-a + f)*
      sym.sinh(2 * c/3 - c*d) - c*sym.exp((a - f)/3)*sym.cosh(c - c*d) +
     c*sym.cosh(2 * c/3 - c*d))/((a - b)*
      sym.exp((a - b)/3) * sym.sinh(c - c*d) + (-a + b)*sym.sinh(2*c/3 - c*d) -
     c*sym.exp((a - b)/3) * sym.cosh(c - c*d) + c*sym.cosh(2*c/3 - c*d))/(41.8/
    1730)

	dfa = sym.diff(mass,a).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])
	dfb = sym.diff(mass,b).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])
	dfc = sym.diff(mass,c).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])
	dfd = sym.diff(mass,d).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])
	dff = sym.diff(mass,f).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])
	dfh = sym.diff(mass,h).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])
	dfg = sym.diff(mass,g).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])

	grad = np.array([dfa,dfb,dfc,dfd,dff,dfh,dfg],dtype = "float64") #ここで、型を浮動小数点で統一する	
	return grad
def delectrons(x): #upsの勾配
	a = sym.Symbol('a')
	b = sym.Symbol('b')
	c = sym.Symbol('c')
	d = sym.Symbol('d')
	f = sym.Symbol('f')
	h = sym.Symbol('h')
	g = sym.Symbol('g')

	a_now = x[0]
	b_now = x[1]
	c_now = x[2]
	d_now = x[3]
	f_now = x[4]
	h_now = x[5]
	g_now = x[6]

	mass = ((a - b)**2 - c**2)/((g - h)**2 - c**2) * \
		sym.sqrt((g*(sym.exp(2*a/3) - 1))/(a*(sym.exp(2*g/3) - 1)))*sym.sqrt((h*(-sym.exp(-2*b/3) + 1))/(b*(-sym.exp(-2*h/3) + 1))) * ((g - h)*
		sym.exp((g - h)/3)*sym.sinh(c/3 - c*d) + (g - h)*sym.sinh(c*d) -
		c*sym.exp((g - h)/3)*sym.cosh(c/3 - c*d) +
		c*sym.cosh(-c*d))/((a - b)*sym.exp((a - b)/3)* sym.sinh(c - c*d) + (-a + b)*
		sym.sinh(2*c/3 - c*d) - c*sym.exp((a - b)/3)* sym.cosh(c - c*d) +
	c*sym.cosh(2*c/3 - c*d))/(0.510999/173000)


	dfa = sym.diff(mass,a).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])
	dfb = sym.diff(mass,b).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])
	dfc = sym.diff(mass,c).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])
	dfd = sym.diff(mass,d).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])
	dff = sym.diff(mass,f).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])
	dfh = sym.diff(mass,h).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])
	dfg = sym.diff(mass,g).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])

	grad = np.array([dfa,dfb,dfc,dfd,dff,dfh,dfg],dtype = "float64") #ここで、型を浮動小数点で統一する	
	return grad
def dmuons(x): #upsの勾配
	a = sym.Symbol('a')
	b = sym.Symbol('b')
	c = sym.Symbol('c')
	d = sym.Symbol('d')
	f = sym.Symbol('f')
	h = sym.Symbol('h')
	g = sym.Symbol('g')

	a_now = x[0]
	b_now = x[1]
	c_now = x[2]
	d_now = x[3]
	f_now = x[4]
	h_now = x[5]
	g_now = x[6]

	mass = ((a - b)**2 - c**2)/((g - h)**2 - c**2)*\
 		sym.sqrt((g*(sym.exp(2*a/3) - 1))/(a*(sym.exp(2*g/3) - 1)))*\
 		sym.sqrt((h*(-sym.exp(-2*b/3) + 1))/(b*(-sym.exp(-2*h/3) + 1))) * ((g - h)*
      	sym.exp((g - h)/3)*sym.sinh(2*c/3 - c*d) + (-g + h)*sym.sinh(c/3 - c*d) -
     	c*sym.exp((g - h)/3)*sym.cosh(2*c/3 - c*d) +
     	c*sym.cosh(c/3 - c*d))/((a - b)*
      	sym.exp((a - b)/3) * sym.sinh(c - c*d) + (-a + b)*sym.sinh(2*c/3 - c*d) - \
     	c*sym.exp((a - b)/3) * sym.cosh(c - c*d) + c*sym.cosh(2*c/3 - c*d))/(105.65/
    	173000)

	dfa = sym.diff(mass,a).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])
	dfb = sym.diff(mass,b).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])
	dfc = sym.diff(mass,c).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])
	dfd = sym.diff(mass,d).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])
	dff = sym.diff(mass,f).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])
	dfh = sym.diff(mass,h).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])
	dfg = sym.diff(mass,g).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])

	grad = np.array([dfa,dfb,dfc,dfd,dff,dfh,dfg],dtype = "float64") #ここで、型を浮動小数点で統一する	
	return grad
def dtauons(x): #upsの勾配
	a = sym.Symbol('a')
	b = sym.Symbol('b')
	c = sym.Symbol('c')
	d = sym.Symbol('d')
	f = sym.Symbol('f')
	h = sym.Symbol('h')
	g = sym.Symbol('g')

	a_now = x[0]
	b_now = x[1]
	c_now = x[2]
	d_now = x[3]
	f_now = x[4]
	h_now = x[5]
	g_now = x[6]

	mass = ((a - b)**2 - c**2)/((g - h)**2 - c**2)*\
	 	sym.sqrt((g*(sym.exp(2*a/3) - 1))/(a*(sym.exp(2*g/3) - 1)))*\
 		sym.sqrt((h*(-sym.exp(-2*b/3) + 1))/(b*(-sym.exp(-2*h/3) + 1))) * ((g - h)*
	 	sym.exp((g - h)/3)*sym.sinh(c - c*d) + (-g + h)*sym.sinh(2*c/3 - c*d) -
	 	c*sym.exp((g - h)/3)*sym.cosh(c - c*d) +
	 	c*sym.cosh(2*c/3 - c*d))/((a - b)*
	 	sym.exp((a - b)/3) *sym.sinh(c - c*d) + (-a + b)*sym.sinh(2*c/3 - c*d) -
	 	c*sym.exp((a - b)/3) *sym.cosh(c - c*d) + c*sym.cosh(2*c/3 - c*d))/(1776.86/
	 	173000)

	dfa = sym.diff(mass,a).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])
	dfb = sym.diff(mass,b).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])
	dfc = sym.diff(mass,c).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])
	dfd = sym.diff(mass,d).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])
	dff = sym.diff(mass,f).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])
	dfh = sym.diff(mass,h).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])
	dfg = sym.diff(mass,g).subs([(a,a_now), (b, b_now),(c,c_now),(d,d_now),(f,f_now),(h,h_now),(g,g_now)])

	grad = np.array([dfa,dfb,dfc,dfd,dff,dfh,dfg], dtype=np.float128) #ここで、型を浮動小数点で統一する	
	return grad

# x = np.array([5,-14.48,8.6,0.28,-280000,34.9,17])
# print(x)
# print(ups(x))
# print(dups(x))
