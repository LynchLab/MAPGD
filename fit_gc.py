import sys
import scipy
import numpy
from scipy.optimize import minimize

def read_data (name):
    data=[]
    File=open(name)
    for line in File:
        line=line.strip('\n').split(',')
        if len(line)==4:
            data.append([int(line[2]), float(line[3])])
    return data

lfac=[0,0]

def log_fac(n):
    if (n>0):
        return scipy.special.gammaln(n+1)
    else:
        return 0
    global lfac
    N=int(n)
    try :
        return lfac[N]
    except:
        for i in xrange(2,n+1): 
            lfac.append(lfac[i-1]+numpy.log(i) )
        return lfac[n]
              
def log_binomial(n,k):
    if n<0 or k<0 or n-k<0:
        print n, k, n-k
    return (log_fac(n)-log_fac(k)-log_fac(n-k))

def init (data):
    global k_min
    k_min=200
    global k_max
    k_max=600
    return [0.97,0.0,-1.18,21.23]
    return [0.1,0.0,5000,0.0] 
 
def objective (param):
    global data
    global s
    global k_min
    global k_max
    ll=0
    x0, b0, x1, b1 =param
    C=-1e300
    for pair in data:
        k, g = pair
        p = x0+b0*g
        r = x1+b1*g
        if (g > 0.3 and g < 0.5 and k > k_min and k < k_max ):
            if(p > 0 and r > 1 and p < 1 ):
                ll+=log_binomial(k+r-1, k) +r*numpy.log( 1.-p )+k*numpy.log( p )
            else:
                ll+=C
#            print p, r, g
    print "Results:", x0, b0, x1, b1, ll
    return -ll
 
#making the data
data=read_data(sys.argv[1])
 
#guessing that parameters from the data (using the simulation parameters would
#be cheeting).
param=init(data)
 
#finding the maximum.
res = minimize(objective, param, method='nelder-mead', options={'xtol': 1e-8, 'disp': True, 'maxiter':10000, 'maxfev':10000} )
 
#and were done. The squaring done here is just to undo a transfermation that is
#done to the parameters in the objective function.
print res
print "Two parameter model: x0 = ", res.x[0], "x1 = ", res.x[1], " log likelihood = ", -res.fun
