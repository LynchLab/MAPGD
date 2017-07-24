import sympy
#import sympy.Q as Query 
import sys
from sympy.utilities.codegen import codegen
from sympy.printing import print_ccode

#The Data
Q=sympy.Symbol('pair.m', positive=True)

MM2=sympy.Symbol('pair.Y_MM', positive=True)
Mm2=sympy.Symbol('pair.Y_Mm', positive=True)
mm2=sympy.Symbol('pair.Y_mm', positive=True)

MM1=sympy.Symbol('pair.X_MM', positive=True)
Mm1=sympy.Symbol('pair.X_Mm', positive=True)
mm1=sympy.Symbol('pair.X_mm', positive=True)

data=[MM2, Mm2, mm2, MM1, Mm1, mm1]
consts=[]

#The Parameters
g_XY=sympy.Symbol('rel.gamma_XY_')
g_YX=sympy.Symbol('rel.gamma_YX_')
T=sympy.Symbol('rel.theta_XY_')
D=sympy.Symbol('rel.Delta_XY_')
d=sympy.Symbol('rel.delta_XY_')
F_Y=sympy.Symbol('rel.f_Y_')
F_X=sympy.Symbol('rel.f_X_')

#An array of the parameters.
params=[F_X, F_Y, T, g_XY, g_YX, d, D]

constants={}
keys=[]

testpow=F_X**2
testadd=F_X+1

class Node:
	def __init__(self, name):
		self.to=[]
		self.fr=[]
		self.lv=-1
		self.name=name
class Dag:
	def __init__(self, size):
		self.size=size
		self.node=[]
		for x in range(0, size):
			self.node.append(Node(x) )
		self.root=[]
	def add_edge(self, x, y):
		self.node[x].to.append(self.node[y])
		self.node[y].fr.append(self.node[x])
	def set_roots(self):
		for x in range(0, self.size):
			if len(self.node[x].fr)==0:
				self.root.append(self.node[x])
	def set_lvs(self):
		for r in self.root:
			set_lv(r, 0)
	def sort(self):
		self.node.sort(key=lambda x: x.name, reverse=True)

def set_lv (b, lv):
	if lv>b.lv:
		b.lv=lv
	for d in b.to:
		set_lv(d, lv+1)

def dag_sort(items):
	print "Haldo!"
	size=len(items)
	dag=Dag(size)
	for x in range(0, size):
		for y in range(0, size):
			if y!=x:
				if expr_in(constants[items[x]], constants[items[y]]):
					dag.add_edge(x, y)
	dag.set_roots()
	dag.set_lvs()
	dag.sort()
	rtn=[]
	for x in range(0, size):
		print dag.node[x].lv
		rtn.append(items[dag.node[x].name])	
	return rtn
	

def expr_in(lhs, rhs):
	for arg in lhs.args:
        	if arg==rhs:
                	return True
               	else:
               		if(expr_in(arg, rhs) ):
				return True
        return False
	
	
def check_subexp():
        global data
        global constants
	global keys
	for key1 in keys:
		exp=constants[key1]
		for key2 in keys:
			print constants
			subed=exact_sub(exp, constants[key2], key2)
			pre(subed)
		

def exact_sub(exp, sub, sym):
        oldargs=exp.args
        newargs=[]
        print exp, sub
	if len(oldargs)>0:
 	       for arg in oldargs:
	                if arg==sub:
	                        print "a sub!!"
	                        newargs.append(sym)
	                else:
	                        newargs.append(exact_sub(arg, sub, sym) )
	       exp=exp.func(*newargs)
        return exp


def pre(expr):
        global data
        global constants
	global keys
	print "checking ", sympy.srepr(expr)
        if list(set(data).intersection(expr.atoms(sympy.Symbol)))==[]:
		print "nul intersection"
		notthere=True
		for key in keys:
			if expr==constants[key]:
				notthere=False
				break
		if(notthere):
			print "not in constants"
       	        	if len(expr.atoms(sympy.Symbol))>1:
				print "not atomic"
                        	keys.append(sympy.Symbol('con.c['+str(len(keys))+']') )
				constants[keys[-1]]=expr
			elif type(expr) is type(testpow):
				print "is pow"
	                        keys.append(sympy.Symbol('con.c['+str(len(keys))+']') )
				constants[keys[-1]]=expr
			elif type(expr) is type(testadd):
				print "is add"
	                        keys.append(sympy.Symbol('con.c['+str(len(keys))+']') )
				constants[keys[-1]]=expr
        else :
                for arg in expr.args:
                        pre(arg)

testme=((MM1*F_X*MM2*F_X)*(4*F_Y-1))
pre(testme)
#check_subexp()
keys=dag_sort(keys)
print keys
for key in keys:
	testme=exact_sub(testme, constants[key], key)
print testme
