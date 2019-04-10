# 与えたパラメータに対して、質量とその勾配を返す
import numpy as np
import sympy as sym

# 以下の値が、"いい感じ" な値らしい
# x = np.array([5,-14.48,8.6,0.28,-280000,34.9,17])

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
	print(c/3 - c*d)
	p1 = (a-b) * np.exp((a-b)/3) * np.sinh(c/3 - c*d) + (a-b)*np.sinh(c*d)
	print("p1=", p1)
	p2 = c * np.exp((a-b)/3) * np.cosh(c/3 - c*d) + c*np.cosh(c*d)
	p3 = (a-b) * np.exp((a-b)/3) * np.sinh(c - c*d)
	p4 = (-(a-b))* np.sinh(2*c/3 - c*d) - c*np.exp((a-b)/3) * np.cosh(c - c*d)
	p5 = c*np.cosh(2*c/3 - c*d)
	p6 = (mex_up/mex_top)
	r1 = (p1-p2)/(p3-p4+p5)/p6
	r2 = ((a-b) * np.exp((a-b)/3) * np.sinh(c/3 - c*d) + (a-b)*np.sinh(c*d) - \
		c * np.exp((a-b)/3) * np.cosh(c/3 - c*d) + c*np.cosh(c*d)) / ((a-b) * np.exp((a-b)/3) * np.sinh(c - c*d)+\
		(-(a-b))* np.sinh(2*c/3 - c*d) - c*np.exp((a-b)/3) * np.cosh(c - c*d)+\
		c*np.cosh(2*c/3 - c*d)) / (mex_up/mex_top)
	print(r1,r2)
	return r2

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



# x = np.array([5,-14.48,8.6,0.28,-280000,34.9,17])
# print(x)
# print(ups(x))
# print(dups(x))




