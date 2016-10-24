#!/usr/bin/python
import random
import math
import numpy
import sys
from scipy.stats import triang

import argparse

parser = argparse.ArgumentParser(description='Makes the files test.txt and test2.txt for IBD simulations.')
parser.add_argument('-N', metavar='--number', type=int, default=5000,
                   help='number of SNPs')
parser.add_argument('-c', metavar='--cov', type=float, default=10,
                   help='simulated depth of coverage')
parser.add_argument('-e', metavar='--error', type=float, default=0.001,
                   help='mean error rate')
parser.add_argument('-s', metavar='--std', type=float, default=0.001,
                   help='stdev error rate')
args = parser.parse_args()

SUM=sum(S)

File=open("test.txt", 'w')
File2=open("test2.txt", 'w')

emean=args.e
estd=args.s
MAX=args.N
N1=args.c
N2=args.c

skipped=0

MINP=0.005
MAXP=0.5

File2.write( "@NAME:SCAFFOLDS	VERSION:0.4.1	FORMAT:TEXT	CONCATENATED\n")
File2.write( "@NAME       	LENGTH\n")
File2.write( "scaffold_1	100000\n")
File2.write( "@END_TABLE\n")
File2.write( "@NAME:GENOTYPES	VERSION:0.4.1	FORMAT:TEXT	CONCATENATED	INDEXED\n")
File2.write( "@SCFNAME	POS	MN_FREQ	Sample_1	Sample_2\n")

A=0
C=0
skipped=0
x=0

while x<MAX+skipped:
	x+=1	
	lastA=A
	lastC=C

	P=0

	while( (P<=MINP) or (P>MAXP) ):
		P = (numpy.random.pareto(0.1, 1)[0]+1 )/1000
#		print P
#		P = round(triang.rvs(0.1), 2)
#		P = round(random.random(), 2)

	Q=1-P

	td=random.random()

#	PA1, PA2=get_conditional(P, lastA, rsq)
	
#	PC1, PC2=get_conditional(P, lastC, rsq)

	mmmm, Mmmm, MMmm, mmMm, MmMm, MMMm, mmMM, MmMM, MMMM=getMs(P, P, S)

	if( min([mmmm, Mmmm, MMmm, mmMm, MmMm, MMMm, mmMM, MmMM, MMMM])<-0.0001):
		skipped+=1
		continue

	Mmmm+=mmmm
	MMmm+=Mmmm

	mmMm+=MMmm
	MmMm+=mmMm
	MMMm+=MmMm

	mmMM+=MMMm
	MmMM+=mmMM

	if (td<mmmm):
		A=0
		C=0

	elif (td<Mmmm):
		A=1
		C=0	

	elif (td<MMmm):
		A=2
		C=0
	
	elif (td<mmMm):
		A=0
		C=1	
	elif (td<MmMm):
		A=1
		C=1	
	elif (td<MMMm):
		A=2
		C=1	
	elif (td<mmMM):
		A=0
		C=2	
	elif (td<MmMM):
		A=1
		C=2	
	else:
		A=2
		C=2	

	M1=0
	m1=0
	e1=0
	M2=0
	m2=0
	e2=0

	if A==2:
		[M1, m1, e1]=numpy.random.multinomial(N1,[(1-e), e/3, 2*e/3])
	elif A==1:
		[M1, m1, e1]=numpy.random.multinomial(N1,[(1-e)/2.0+e/6, (1-e)/2.0+e/6, 2*e/3])
	else:
		[M1, m1, e1]=numpy.random.multinomial(N1,[e/3.0, (1-e), 2*e/3])
	if C==2:
		[M2, m2, e2]=numpy.random.multinomial(N2,[(1-e), e/3, 2*e/3])
	elif C==1:
		[M2, m2, e2]=numpy.random.multinomial(N2,[(1-e)/2.0+e/6, (1-e)/2.0+e/6, 2*e/3])
	else:
		[M2, m2, e2]=numpy.random.multinomial(N2,[e/3.0, (1-e), 2*e/3])

	if (sigma>0):
		den = numpy.random.binomial(sigma, P)
		sampP = round(float(den)/float(sigma), 2)
		if (l2i):
			P=(float(A+C)+float(den))/float(sigma+4)
			Q=1-P
		else:
			P=sampP
			Q=1-P
	if P==0:
		skipped+=1
		continue

	MM1=lnc*M1+lne*(m1+e1)#-math.log(P**2)
	Mm1=lnch*(M1+m1)+lne*e1#-math.log(2*P*Q)
	mm1=lnc*m1+lne*(M1+e1)#-math.log(Q**2)

	MM2=lnc*M2+lne*(m2+e2)#-math.log(P**2)
	Mm2=lnch*(M2+m2)+lne*e2#-math.log(2*P*Q)
	mm2=lnc*m2+lne*(M2+e2)#-math.log(Q**2)

	norm1=math.log(math.exp(MM1)+math.exp(Mm1)+math.exp(mm1) )
	norm2=math.log(math.exp(MM2)+math.exp(Mm2)+math.exp(mm2) )

	File2.write( '\t'.join(map(str, ["scaffold_1", x-skipped, 1-P] ) )+'\t')
	File.write( '\t'.join(map(str, [-MM1+norm1, -Mm1+norm1, -mm1+norm1, N1, -MM2+norm2, -Mm2+norm2, -mm2+norm2, N2, P] ) )+'\n')
	File2.write( '/'.join(map(str, [-MM1+norm1, -Mm1+norm1, -mm1+norm1, int(N1)] ) )+'\t')
	File2.write( '/'.join(map(str, [-MM2+norm2, -Mm2+norm2, -mm2+norm2, int(N2)] ) )+'\n')
File.close()
dig=4
print  round(FA, dig), round(FC, dig), round(r, dig), round(sA, dig), round(sC, dig), round(z1, dig), round(z2, dig)
File2.write( "@END_TABLE\n")
