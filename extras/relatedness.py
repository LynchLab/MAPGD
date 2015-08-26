import sys
import math
from scipy import stats
import scipy
import numpy as np
from scipy.optimize import minimize
import rml
import random

def getdll (model, params):
	ret=[]
	maxll=model(params)
	z=0
	for x in range(0, 8):
	#	ret.append(getci(model, params, x-z, 'l', maxll) )
	#	ret.append(getci(model, params, x-z, 'u', maxll) )
		ret.append(getci(model, params, x-z, 'dll', maxll) )
	return ret


class oneparam:
	def __init__ (self, model, params, i, maxll, b ):
		self.p=np.copy(params)
		self.i=i
		self.maxll=maxll
		self.mode=(b=='u')
		self.init=np.copy(params[i])
		self.model=model
	def get (self, x):
		if self.mode:
			self.p[self.i]=self.init+abs(x)
		else:
			self.p[self.i]=self.init-abs(x)
		return abs(self.model(self.p)-self.maxll-3.84)

class eightparam:
        def __init__ (self, model, params, i, maxll, b ):
                self.p=np.copy(params)
                self.i=i
                self.maxll=maxll
                self.free=np.copy(params).tolist()
		self.free=self.free[0:i]+self.free[i+1:]
                self.model=model
        def get (self, x):
                return self.model(np.insert(x, self.i, 0) )

		
def getci (model, params, x, b, maxll):
	OneModel=oneparam(model, params, x, maxll, b)
	if b=='l':
		q=minimize(OneModel.get, [0.01], method='SLSQP')
		return params[x]-abs(q.x[0])
	if b=='u':
		q=minimize(OneModel.get, [0.01], method='SLSQP')
		return params[x]+abs(q.x[0])
	if b=='dll':
		EightModel=eightparam(model, params, x, maxll, b)
		q=minimize(EightModel.get, EightModel.free, method='SLSQP')
		return maxll-q.fun #EightModel.get(EightModel.free)

def Model(params):
	fA, fC, r, sA, sC, z1, z2=params
#	e, z1, z2=[0, 0, 0]
	rml.fullModel(0, fA, fC, r, sA, sC, z1, z2)
	return (-rml.get_ll())

def FullModel(params):
	e, fA, fC, r, sA, sC, z1, z2=params
	rml.fullModel(0, fA, fC, r, sA, sC, z1, z2)
	return (-rml.get_ll()+abs(e) )

size=rml.getsize(sys.argv[1])

R=[]
print "total sites, (line A:min-max, line C:min-max), success, e, dll, fA, dll, fC, dll, r, dll, sAC, dll, sCA, dll, z1, dll, z2, dll, null_ll, fit_ll"

model=FullModel

for x in range(0, size):
	for y in range(x+1, size):
		rml.read(sys.argv[1], x, y)
#		params=rml.estimate()
		params=[0.,0.,0.,0.25,0.,0.,0.,0.25]
		SLSQP=[]
		#for x in range (0, 10):
		#	params=[0.,0.,0.,random.random(),0.,0.,0., random.random()]

		SLSQP.append( minimize(model, params, method='SLSQP', options={'ftol': 1e-11}) )
		#	BFGS = minimize(model, params, method='BFGS', options={'gtol': 1e-11} )
		#	L_BFGS_B = minimize(model, params, method='L-BFGS-B', options={'gtol': 1e-11} )
		#	COBYLA = minimize(model, params, method='COBYLA')
#			NEWTON-CG = minimize(model, params, method='DOGLEG')


		#	print "SL", SLSQP.fun, 		SLSQP.x[1:4] 
		#	print "BF", BFGS.fun,		BFGS.x[1:4]
		#	print "L_", L_BFGS_B.fun,	L_BFGS_B.x[1:4]
		#	print "CO", COBYLA.fun,		COBYLA.x[1:4]
		sorted(SLSQP, key=lambda this : this.fun)
		results=SLSQP[0]
		#quit()
		fit=np.copy(results.x).tolist() 
		dlls=getdll(model, np.copy(results.x) )
		out=[]
		null_ll=FullModel( [0.,0.,0.,0.,0.,0.,0.,0.] ) 
		out.append(str(results.success))
		fit_ll=model( results.x ) 
		best_ll=FullModel( [0.,0.,0.,0.25,0.,0.,0.,0.25] )
		for z in range(0, 8):
			fit[z]='{:.4f}'.format(fit[z])
			dlls[z]='{:.4f}'.format(dlls[z])
		for z in range(0, 8):
			out.append(fit[z]+", "+dlls[z])
		out.append('{:.4f}'.format(null_ll) )
		out.append('{:.4f}'.format(fit_ll) )
		print ", ".join(out)
		print (best_ll-fit_ll)
		print (null_ll-fit_ll)
		quit()
