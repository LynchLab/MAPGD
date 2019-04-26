import sympy
import sys
from sympy.utilities.codegen import codegen
from sympy.printing import print_ccode

#The Data
def set_data(new_data):
	global data
	data=new_data
def set_params(new_params):
	global params
	params=new_params
	global testpow
	testpow=params[0]**2
	global testadd
	testadd=params[0]+1

data=[]
params=[]

constants={}
keys=[]

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
                self.node.sort(key=lambda x: x.lv)

def set_lv (b, lv):
        if lv>b.lv:
                b.lv=lv
        for d in b.to:
                set_lv(d, lv+1)

def dag_sort(items):
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
                rtn.append(items[dag.node[x].name])
        return rtn

def dag_replace(items):
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
	#we need to start at root here
        for x in range(0, size):
        	for y in range(x+1, size):
			constants[items[x]]=exact_sub(constants[items[x]], constants[items[y]], items[y])
        rtn=[]
        for x in reversed(range(0, size)):
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
        if len(oldargs)>0:
               for arg in oldargs:
                        if arg==sub:
                                newargs.append(sym)
                        else:
                                newargs.append(exact_sub(arg, sub, sym) )
               exp=exp.func(*newargs)
        return exp

def pre(expr):
        global data
        global constants
        global keys
        if list(set(data).intersection(expr.atoms(sympy.Symbol)))==[]:
                notthere=True
                for key in keys:
                        if expr==constants[key]:
                                notthere=False
                                break
                if(notthere):
                        if len(expr.atoms(sympy.Symbol))>1:
                                keys.append(sympy.Symbol('con.c['+str(len(keys))+']') )
                                constants[keys[-1]]=expr
                        elif type(expr) is type(testpow):
                                keys.append(sympy.Symbol('con.c['+str(len(keys))+']') )
                                constants[keys[-1]]=expr
                        elif type(expr) is type(testadd):
                                keys.append(sympy.Symbol('con.c['+str(len(keys))+']') )
                                constants[keys[-1]]=expr

        else :
                for arg in expr.args:
			pre(arg)

