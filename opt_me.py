import numpy as np
from scipy.optimize import minimize
import particle_mass as pmas

def func(x):
	# f = (1 - pmas.ups(x))**2+(1 - pmas.charms(x))**2+(1 - pmas.downs(x))**2+(1 - pmas.stranges(x))**2+(1 - pmas.botoms(x))**2+(1 - pmas.electrons(x))**2+(1 - pmas.muons(x))**2+(1 - pmas.tauons(x))**2
	f = np.log(pmas.ups(x))**2+np.log(pmas.charms(x))**2+np.log(pmas.downs(x))**2+np.log(pmas.stranges(x))**2+np.log(pmas.botoms(x))**2+np.log(pmas.electrons(x))**2+np.log(pmas.muons(x))**2+np.log(pmas.tauons(x))**2
	return f
def gradient(x):
	# df = -2 * (1 - pmas.ups(x)) * pmas.dups(x)-2 * (1 - pmas.charms(x)) * pmas.dcharms(x)-2 * (1 - pmas.downs(x)) * pmas.ddowns(x)-2 * (1 - pmas.stranges(x)) * pmas.dstranges(x)-2 * (1 - pmas.botoms(x)) * pmas.dbotoms(x)-2 * (1 - pmas.electrons(x)) * pmas.delectrons(x)-2 * (1 - pmas.muons(x)) * pmas.dmuons(x)-2 * (1 - pmas.tauons(x)) * pmas.dtauons(x)
	df = 2*np.log(pmas.ups(x))*pmas.dups(x)/pmas.ups(x)+2*np.log(pmas.charms(x))*pmas.dcharms(x)/pmas.charms(x)+2*np.log(pmas.downs(x))*pmas.ddowns(x)/pmas.downs(x)+2*np.log(pmas.stranges(x))*pmas.dstranges(x)/pmas.stranges(x)+2*np.log(pmas.botoms(x))*pmas.dbotoms(x)/pmas.botoms(x)+2*np.log(pmas.electrons(x))*pmas.delectrons(x)/pmas.electrons(x)+2*np.log(pmas.muons(x))*pmas.dmuons(x)/pmas.muons(x)+2*np.log(pmas.tauons(x))*pmas.dtauons(x)/pmas.tauons(x)
	# df = np.array(df)
	# print(df)
	return df
init = np.array([5.0,-14.48,8.6,0.28,-28.0,34.9,17.0],dtype=np.float128)
#init = np.array([0,0,0,0,0,0,0],dtype=np.float128)

res1 = minimize(func, init, jac=gradient,
                  method='SLSQP', options={"maxiter":200})
print(res1)
x = res1.x
print("befor=",func(init)) 
print("after=",func(x))
print('up =',pmas.ups(x))
print('charm =',pmas.charms(x))
print('down =',pmas.downs(x))
print('strange =',pmas.stranges(x))
print('botom =',pmas.botoms(x))
print('electron =',pmas.electrons(x))
print('muon =',pmas.muons(x))
print('tauon =',pmas.tauons(x))