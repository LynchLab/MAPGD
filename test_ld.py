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

        else :
                for arg in expr.args:
			pre(arg)
#The three componenents of the likelihood function. H00 is homozygous for the major allele, H01 heterozygous, and H11 homozygous minor.

P=1-Q
#sympy.Q.positive(P)
#var=P*(1.-P)
#std=var**(0.5);
#skew=(1.-2.*P)/std;
#kurtosis=1./(1.-P)+1./P-3.;

#mmmmi=(Q**4+(F_X+F_Y+4*T)*var*Q**2-2*(g_XY+g_YX)*skew*Q+d*kurtosis+D*var**2);
#MMMMi=(P**4+(F_X+F_Y+4*T)*var*P**2-2*(g_XY+g_YX)*skew*P+d*kurtosis+D*var**2);
#Mmmmi=(2*P*Q**3+2*(F_Y*P*Q-F_X*Q**2-2*T*(Q**2-P*Q))*var-4*g_XY*(P-Q)*skew+2*g_YX*(Q-P)*skew-2*d*kurtosis-2*D*var**2);
#mmMmi=(2*P*Q**3+2*(F_X*P*Q-F_Y*Q**2-2*T*(Q**2-P*Q))*var-4*g_YX*(P-Q)*skew+2*g_XY*(Q-P)*skew-2*d*kurtosis-2*D*var**2);
#MmMMi=(2*Q*P**3+2*(F_Y*P*Q-F_X*P**2-2*T*(P**2-P*Q))*var-4*g_XY*(Q-P)*skew+2*g_YX*(P-Q)*skew-2*d*kurtosis-2*D*var**2);
#MMMmi=(2*Q*P**3+2*(F_X*P*Q-F_Y*P**2-2*T*(P**2-P*Q))*var-4*g_YX*(Q-P)*skew+2*g_XY*(P-Q)*skew-2*d*kurtosis-2*D*var**2);
#MmMmi=(4*P**2*Q**2+4*(T*(P**2+Q**2-2*P*Q)-F_X*P*Q-F_Y*P*Q)*var+4*g_XY*(P-Q)*skew+4*g_YX*(Q-P)*skew+4*d*kurtosis+4*D*var**2);
#MMmmi=(P**2*Q**2+(F_Y*P**2+F_X*Q**2-4*T*P*Q)*var+2*g_XY*P*skew-2*g_YX*Q*skew+d*kurtosis+D*var**2);
#mmMMi=(P**2*Q**2+(F_X*P**2+F_Y*Q**2-4*T*P*Q)*var+2*g_YX*P*skew-2*g_XY*Q*skew+d*kurtosis+D*var**2);

A=P
C=P

Va=A*(1.-A)#;    //variance of the two haploid genomes of A 
Vc=C*(1.-C)#;    // "   "        "       "       "     of B
Sa=sympy.sqrt(Va)#;    //standard deviation of haploid genomes A
Sc=sympy.sqrt(Vc)#;    // and "        "       "       "       C.

E_A2  =(F_X*Va+2.*(A**2)+Va)/2.#; //Expectation of A^2
E_C2  =(F_Y*Vc+2.*(C**2)+Vc)/2.#; //
E_AC  =T*Sa*Sc+A*C;
ga=(1.-2.*A)/Sa;
gc=(1.-2.*C)/Sc;
E_A2C =(g_XY*Va*Sc*ga+A*A*C+Va*(T*(1.+2.*A)+F_X*C+C/(1-C) ) )/2.;
E_AC2 =(g_YX*Vc*Sa*gc+C*C*A+Vc*(T*(1.+2.*C)+F_Y*A+A/(1-A) ) )/2.;
ka=1./(1.-A)+1./A-3.;
kc=1./(1.-C)+1./C-3.;
E_A2C2=(d*sympy.sqrt(ka*kc)+D)*Va*Vc+A*A*C*C+F_X*Va*C*C+F_Y*Vc*A*A+4.*T*Sa*Sc*A*A+C*2.*g_XY*Va*Sc*ga+2.*A*g_YX*Vc*Sa*gc;

e=0

mmmm=1-6*P+0.0*e+2*E_A2+2*E_C2+8.0*E_AC-4*E_A2C-4*E_AC2+1*E_A2C2;
Mmmm=0+4*P+2.0*e-4*E_A2+0*E_C2-10.*E_AC+8*E_A2C+4*E_AC2-2*E_A2C2;
MMmm=0-1*P-0.5*e+2*E_A2+0*E_C2+2.0*E_AC-4*E_A2C-0*E_AC2+1*E_A2C2;
mmMm=0+4*P-2.0*e+0*E_A2-4*E_C2-10.*E_AC+4*E_A2C+8*E_AC2-2*E_A2C2;
MmMm=0+0*P+0.0*e+0*E_A2+0*E_C2+12.*E_AC-8*E_A2C-8*E_AC2+4*E_A2C2;
MMMm=0+0*P+0.0*e+0*E_A2+0*E_C2-2.0*E_AC+4*E_A2C+0*E_AC2-2*E_A2C2;
mmMM=0-1*P+0.5*e+0*E_A2+2*E_C2+2.0*E_AC+0*E_A2C-4*E_AC2+1*E_A2C2;
MmMM=0+0*P+0.0*e+0*E_A2+0*E_C2-2.0*E_AC+0*E_A2C+4*E_AC2-2*E_A2C2;
MMMM=0+0*P+0.0*e+0*E_A2+0*E_C2+0.0*E_AC+0*E_A2C+0*E_AC2+1*E_A2C2;

#S=(1-(mmmmk+MMMMk+Mmmmk+mmMmk+MmMMk+MMMmk+MmMmk+MMmmk+mmMMk))**2
S=1

#The log likelihood equation
lnL=sympy.log( mmmm*mm1*mm2+mmMm*mm1*Mm2+mmMM*mm1*MM2+Mmmm*Mm1*mm2+MmMm*Mm1*Mm2+MmMM*Mm1*MM2+MMmm*MM1*mm2+MMMm*MM1*Mm2+MMMM*MM1*MM2 )

system_eq=[]

#We first need the three equations we are going to try and set to zero, i.e. the first partial derivitives wrt e h and F.
print "/*This code was automatically generated by "+str(sys.argv[0])+"*/\n"
print "#include \"allele.h\""
print "#include \"quartet.h\""
print "#include \"typedef.h\""
print "#include \"constants.h\""
print 
#numpy.set_printoptions(precission=18)


for x in range(0, 7):
	system_eq.append(sympy.diff(lnL, params[x]) )
	print "inline float_t H"+str(x)+" (const Genotype_pair &pair, const Constants <float_t, const std::pair<const Genotype_pair &, const Relatedness &> > &con) {"
	sys.stdout.write("\treturn ")
	out=(system_eq[-1])
	out=sympy.simplify(out)
	pre(out)
	keys=dag_sort(keys)
	for key in keys:
		out=exact_sub(out, constants[key], key)
	string=sympy.printing.ccode(out)
	print string, ";\n}\n"

#Then we need to make the Jacobian, which is a matrix with ...
for x in range(0, 7):
	for y in range(0, 7):
		print "inline float_t J"+str(x)+str(y)+" (const Genotype_pair &pair, const Constants <float_t, const std::pair <const Genotype_pair &, const Relatedness &> > &con) {"
		sys.stdout.write("\treturn ")
                out=(sympy.diff(system_eq[x], params[y]))
		#out=sympy.simplify(out)
		pre(out)
		keys=dag_sort(keys)
		for key in keys:
			out=exact_sub(out, constants[key], key)
	        string=sympy.printing.ccode(out)
		print string, ";\n}\n"

print "inline float_t lnL_NR (const Genotype_pair &pair, const Relatedness &rel) {"
sys.stdout.write("\treturn ")
print_ccode(sympy.simplify(lnL))
print ";\n}\n"

keys=dag_replace(keys)

for key in keys:
	print "inline void set_"+str(key).split('.')[1].replace('[','').replace(']','') +" (Constants <float_t, const std::pair <const Genotype_pair &, const Relatedness &> > &con, const std::pair <const Genotype_pair &, const Relatedness &> &d ) {"
	print "\tconst Genotype_pair *pair=&d.first;"
	print "\tconst Relatedness *rel=&d.second;"
	print "\tcon.c["+str(key).split('.')[1].replace('[','').replace(']','')+"]=",
        string=sympy.printing.ccode(constants[key])
        string=string.replace("rel.", "rel->").replace("pair.","pair->")
	print string, ";\n}\n"

print "static void (**cfn)(Constants <float_t, const std::pair <const Genotype_pair &, const Relatedness &> > &, const std::pair <const Genotype_pair &, const Relatedness &> &)["+str(len(constants))+"]={",
for key in keys:
	print ", "+"&set_c"+str(key).split('.')[1].replace('[','').replace(']','')
print "};"

print "#define REL_CNTS\t"+str(len(keys) )
print "#define REL_ARRAY\tcfn"

quit()
