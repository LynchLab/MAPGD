#I honestly have no idea what this script is intended to show. I guess that proview doesn't crash?
#!/bin/python
import random
import math
import numpy

e=0.001
P=0.1
Q=0.1
N=200
LEN=10000

for x in range(0, LEN):
	[A1, T1, C1, G1]=numpy.random.multinomial(N,[(1-e)*P+(1-P)*e/3, (1-e)*(1-P)+P*e/3, 2*e/3] )
	[a1, t1, c1, g1]=numpy.random.multinomial(N,[(1-e)*P+(1-P)*e/3, (1-e)*(1-P)+P*e/3, 2*e/3] )
	print '\t'.join(["scaffold_1", str(x), "N", str(P1+Q1+E1), random.shuffle('T'*P1+'A'*Q1+'G'*E1), 'F'*(P1+Q1+E1+p1+q1+), str(P2+Q2+E2), random.shuffle('T'*P2+'A'*Q2+'G'*E2), 'F'*(P2+Q2+E2)])
