import sys
import math
from scipy import stats
import scipy
import numpy as np
from scipy.optimize import minimize
from scipy.optimize import basinhopping
from scipy.optimize import fmin_ncg
import rml
import random

class my_results:
	def __init__(self, func):
		self.func=func

def my_minimize (model, params):
#	return scipy.optimize.minimize(model, params, method='SLSQP').x
	nR=np.copy(map(float, params) )
	rml.set_max_P(1.0)
	ret=my_results(0)
	ret.niter=0
	while True:
#		print params
		if(rml.get_max_P()<0.75):	
			rml.set_max_P(0.75)
		ret.jac=Model_jac(params)	
		ret.hess=Model_hess(params)	
		if (np.linalg.det(ret.hess)==0):
			for x in range(0, len(params) ):
				params[x]=params[x]/2.
			continue
		R=np.linalg.inv( np.matrix(ret.hess ) )*np.matrix( ret.jac ).transpose()
		R=np.array(R.transpose() )[0]
		for x in range(0, len(R) ):
			nR[x]=params[x]-R[x]
		params=np.copy(nR)
		if (max ( abs(R) ) < 10**-6):
			break
		ret.niter+=1
		if (ret.niter>99):
			ret.success=False
			break
	ret.x=params
	ret.func=model(params)
	return params

def my_dll (model, params, x):
#	return scipy.optimize.minimize(model, params, method='SLSQP').fun
	ret=my_results(model(params) )
	ret.x=params
	params=[0.,0.,0.,0.,0.,0.,0.]
	nR=np.copy(map(float, params) )
	ret.niter=0
	P=rml.get_max_P()
	while True:
		if(rml.get_max_P()<0.75):	
			rml.set_max_P(0.75)

		ret.hess=Model_hess(params)
		ret.jac=Model_jac(params)

		ret.hess=np.delete(ret.hess, (x), axis=0)
		ret.hess=np.delete(ret.hess, (x), axis=1)
		ret.jac=np.delete(ret.jac, (x), axis=0)

		if (np.linalg.det(ret.hess)==0):
			for x in range(0, len(params) ):
				params[x]=params[x]/2.
			continue
		try:
			R=np.linalg.inv( np.matrix(ret.hess ) )*np.matrix( ret.jac ).transpose()
		except:
			print params
			print ret.hess
			print np.linalg.inv(ret.hess)
			print ret.jac
			quit()
		R=np.array(R.transpose() )[0]

		for y in range(0, len(params) ):
			if y<x:
				nR[y]=params[y]-R[y]
			elif y>x:
				nR[y]=params[y]-R[y-1]
		params=np.copy(nR)
		if (max ( abs(R) ) < 10**-4):
			break
		ret.niter+=1
		if (ret.niter>99):
			ret.success=False
			print "DANM!"
			break
	ret.func=model(params)-model(ret.x)
	ret.x=params
	rml.set_max_P(P)
	return ret.func


def getdll (model, params):
	ret=[]
	maxll=model(params)
	z=0
	for x in range(0, len(params) ):
	#	ret.append(getci(model, params, x-z, 'l', maxll) )
	#	ret.append(getci(model, params, x-z, 'u', maxll) )
		ret.append(my_dll(model, params, x) )
	return ret

def Model(params):
	fA, fC, r, sA, sC, z1, z2=params
	rml.fullModel(0, fA, fC, r, sA, sC, z1, z2)
	return (-rml.get_ll() )

def Model_jac(params):
	fA, fC, r, sA, sC, z1, z2=params
	rml.fullModel(0, fA, fC, r, sA, sC, z1, z2)
	return (np.multiply( rml.get_jac(), 1) )

def Model_hess(params):
	fA, fC, r, sA, sC, z1, z2=params
	rml.fullModel(0, fA, fC, r, sA, sC, z1, z2)
	return ( np.multiply( np.array(rml.get_hess()).reshape(7,7), 1) )

def FullModel(params):
	e, fA, fC, r, sA, sC, z1, z2=params
	rml.fullModel(e, fA, fC, r, sA, sC, z1, z2)
	return (-rml.get_ll() )

size=rml.getsize(sys.argv[1])

R=[]
print "total sites, line A, min, max, line C, min, max, success, e, dll, fA, dll, fC, dll, r, dll, sAC, dll, sCA, dll, z1, dll, z2, dll, null_ll, fit_ll, max_P"

model=Model

for x in range(0, size):
	for y in range(x+1, size):
		cov=rml.read(sys.argv[1], x, y)
		out=map(str, cov) 
		if cov[0]!=0:
			params=[0.,0.,0.,0.,0.,0.,0.,0.]

			maximum=my_minimize(model, params[1:] )
			fit=np.copy(maximum).tolist() 
			dlls=getdll( model, np.copy(maximum) )
			fit=[0]+fit
			dlls=[0]+dlls
			null_ll=Model( [0.,0.,0.,0.,0.,0.,0.] ) 
			fit_ll=model( maximum ) 
			out.append("true")
			for z in range(0, 8):
				fit[z]='{:.4f}'.format(fit[z])
				dlls[z]='{:.4f}'.format(dlls[z])
			for z in range(0, 8):
				out.append(fit[z]+", "+dlls[z])
			out.append('{:.4f}'.format(null_ll) )
			out.append('{:.4f}'.format(fit_ll) )
			out.append( str(rml.get_max_P()) )
			print ", ".join(out)
