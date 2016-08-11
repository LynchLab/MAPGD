import random

SUM={}
NSUM={}

for x in range (0, 10000):
	A=random.randint(0, 10)
	B=random.randint(0, 10)
	C=random.randint(0, 10)
	D=random.randint(0, 10)
	Z=(A+B)
	c=int(C<A)
	d=int(D<A)
	try:
		SUM[Z]+=c
		SUM[Z]+=d
		NSUM[Z]+=2
	except:
		SUM[Z]=c+d
		NSUM[Z]=2
	print c, d, A, Z, float(SUM[Z])/float(NSUM[Z])

#for x in range (0, 20):
	
