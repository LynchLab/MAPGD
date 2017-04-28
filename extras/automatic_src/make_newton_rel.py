import sympy 
import sys
from sympy.utilities.codegen import codegen
from sympy.printing import print_ccode

#The Data
P=sympy.Symbol('P')
Q=1-P

MM1=sympy.Symbol('pair.X_MM')
Mm1=sympy.Symbol('Mm1')
mm1=sympy.Symbol('mm1')

MM2=sympy.Symbol('MM2')
Mm2=sympy.Symbol('Mm2')
mm2=sympy.Symbol('mm2')

#The Parameters
F_x=sympy.Symbol('rel.Delta_XY_')
F_y=sympy.Symbol('F_y')
T=sympy.Symbol('T')
g_xy=sympy.Symbol('g_xy')
g_yx=sympy.Symbol('g_yx')
D=sympy.Symbol('D')
d=sympy.Symbol('d')

params=[F_x, F_y, T, g_xy, g_yx, D, d]

#The three componenents of the likelihood function. H00 is homozygous for the major allele, H01 heterozygous, and H11 homozygous minor.

M2=
M3=
M4=

H00=Q**4+(F_x+F_y+4*T)*M2*Q**2-2*(g_xy+g_yx)*M3*Q+d*M4+D*M2**2
H01=
H02=
H10=2*P*Q**3+2*(F_y*P*Q-F_x*Q**2-2*T*(Q**2-P*Q) )*M2+4*g_xy*M3*Q+2*g_yx*(Q-P)*M3-2*d*M4-2*d*M2**2
H11=4*P**2*Q**2+4*(T*(P**2+Q**2-2*P*Q)-(F_x+F_y)*(P*Q) )*M2+4*(g_xy+g_yx)*(P-Q)*M3+4*d*M4+4*d*M2**2
H12=
H20=
H21=
H22=

#The log likelihood equation
lnL=sympy.log( H00*sympy.exp(-lmm1-lmm2)+H01*(sympy.exp(-lmm1-lMm2))+H02*sympy.exp(-lmm1-lMM2)+H10*sympy.exp(-lMm1-lmm2)+H11*sympy.exp(-lMm1-lMm2)+H12*sympy.exp(-lMm1-lMM2)+H20*sympy.exp(-lMM1-lmm2)+H21*sympy.exp(-lMM1-lMm2)+H22*sympy.exp(-lMM1-lMM2) )
#lnL=-e**2-h**2-F**2

system_eq=[]

#We first need the three equations we are going to try and set to zero, i.e. the first partial derivitives wrt e h and F.
print "/*This code was automatically generated by "+str(sys.argv[0])+"*/\n"
print "#include \"allele.h\""
print "#include \"quartet.h\""
print "#include \"typedef.h\""
print 
#numpy.set_printoptions(precission=18)
for x in range(0, 7):
	system_eq.append(sympy.diff(lnL, params[x]) )
	print "inline float_t H"+str(x)+" (const Genotype_pair &pair, const Relatedness &rel) {"
	sys.stdout.write("\treturn ")
	print_ccode(  system_eq[-1] )
	print ";\n}\n"

#Then we need to make the Jacobian, which is a matrix with ...
for x in range(0, 7):
	for y in range(0, 7):
		print "inline float_t J"+str(x)+str(y)+" (const Genotype_pair &pair, const Relatedness &rel) {"
		sys.stdout.write("\treturn ")
		print_ccode(sympy.diff(system_eq[x], params[y]))
		print ";\n}\n"

print "inline float_t lnL_NR (const quartet_t &q, const Allele &a) {"
sys.stdout.write("\treturn ")
print_ccode(sympy.simplify(lnL) )
print ";\n}\n"

quit()
