import random
import math
import numpy

e=0.001
P=0.1
Q=0.1
N=200

for x in range(0, 100000):
	[P1, Q1, E1]=numpy.random.multinomial(N,[(1-e)*P+(1-P)*e/3, (1-e)*(1-P)+P*e/3, 2*e/3] )
	[P2, Q2, E2]=numpy.random.multinomial(N,[(1-e)*Q+(1-Q)*e/3, (1-e)*(1-Q)+Q*e/3, 2*e/3] )
	print '\t'.join(["scaffold_1", str(x), "N", str(P1+Q1+E1), 'T'*P1+'A'*Q1+'G'*E1, 'F'*(P1+Q1+E1), str(P2+Q2+E2), 'T'*P2+'A'*Q2+'G'*E2, 'F'*(P2+Q2+E2)])
