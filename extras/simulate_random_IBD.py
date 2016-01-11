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
                   help='simulated deapth of coverage')
parser.add_argument('-s', metavar='--sigma', type=float, default=0.0,
                   help='simulated sampling error')
args = parser.parse_args()



def getMs (P, S):
	FA, FC, r, sA, sC, z1, z2=S
	Q=1-P

        V=P*Q 
        sigma=math.sqrt(V)

        E_A2  =(FA*V+2.*P**2+V)/2.     
        E_C2  =(FC*V+2.*P**2+V)/2.     
        E_AC  =r*sigma*sigma+P**2               
        g=(1.-2.*P)/sigma   
       	E_A2C =(sA*V*sigma*g+P**3+V*(r*(1+2.*P)+FA*P+P/(1-P) ) )/2.
        E_AC2 =(sC*V*sigma*g+P**3+V*(r*(1+2.*P)+FC*P+P/(1-P) ) )/2.
        k=1./(1.-P)+1./P-3.
        E_A2C2=(z1*k+z2)*V*V+P*P*P*P+FA*V*P*P+FC*V*P*P+4*r*sigma*sigma*P*P+P*2*sA*V*sigma*g+2*P*sC*V*sigma*g

        mmmm=1-6*P+2*E_A2+2*E_C2+8.0*E_AC-4*E_A2C-4*E_AC2+1*E_A2C2
        Mmmm=0+4*P-4*E_A2+0*E_C2-10.*E_AC+8*E_A2C+4*E_AC2-2*E_A2C2
        MMmm=0-1*P+2*E_A2+0*E_C2+2.0*E_AC-4*E_A2C-0*E_AC2+1*E_A2C2
        mmMm=0+4*P+0*E_A2-4*E_C2-10.*E_AC+4*E_A2C+8*E_AC2-2*E_A2C2
        MmMm=0+0*P+0*E_A2+0*E_C2+12.*E_AC-8*E_A2C-8*E_AC2+4*E_A2C2
        MMMm=0+0*P+0*E_A2+0*E_C2-2.0*E_AC+4*E_A2C+0*E_AC2-2*E_A2C2
        mmMM=0-1*P+0*E_A2+2*E_C2+2.0*E_AC+0*E_A2C-4*E_AC2+1*E_A2C2
        MmMM=0+0*P+0*E_A2+0*E_C2-2.0*E_AC+0*E_A2C+4*E_AC2-2*E_A2C2
        MMMM=0+0*P+0*E_A2+0*E_C2+0.0*E_AC+0*E_A2C+0*E_AC2+1*E_A2C2

	return mmmm, Mmmm, MMmm, mmMm, MmMm, MMMm, mmMM, MmMM, MMMM

def test (S):
	for x in range (1, 10):
		mmmm, Mmmm, MMmm, mmMm, MmMm, MMMm, mmMM, MmMM, MMMM=getMs(float(x)/10., S)
		if ( min([mmmm, Mmmm, MMmm, mmMm, MmMm, MMMm, mmMM, MmMM, MMMM])<-0.0001):
			return False
	return True

Fa=0
Fc=1
r=2
sac=3
Sca=4
z1=5
z2=6

Z= [0,1,2,3,4,5,6]

S=[random.random()-0.5, random.random()-0.5, random.random()-0.5, random.random()-0.5, random.random()-0.5, random.random()-0.5, random.random()-0.5 ]

random.shuffle(Z)

for x in range(0, random.randint(0, 6) ):
	S[Z[x]]=0      

SUM=sum(S)

File=open("test.txt", 'w')
File2=open("test2.txt", 'w')

e=0.001
lnc=math.log(1.0-e)
lnch=math.log((1.0-e)/2.0+e/6)
lne=math.log(e)

count=[]
for z in range (0, 11):
	count.append([])
	for x in range(0, 3):
		count[z].append([])
		for y in range(0, 3):
			count[z][x].append(0)

skipped=0

MAX=args.N
N1=args.c
N2=args.c
sigma=args.s

MINP=0.1
MAXP=0.4

while ( not ( test(S) ) ):
	S=[random.random()-0.5, random.random()-0.5, random.random()-0.5, random.random()-0.5, random.random()-0.5, random.random()-0.5, random.random()-0.5 ]
	for x in range(0, random.randint(0, 6) ):
		S[Z[x]]=0
FC, FA, r, sC, sA, z1, z2=S
for x in range(0, MAX):
	P=0
	while(P<MINP or P>MAXP):
		P = round(triang.rvs(0.1), 2)
	Q=1-P
	td=random.random()

	mmmm, Mmmm, MMmm, mmMm, MmMm, MMMm, mmMM, MmMM, MMMM=getMs(P, S)

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
		OP=P
		den = numpy.random.binomial(sigma, P) 
		if den==0:
			den=1
		if den==sigma:
			den=sigma-1
		P = round(float(den)/float(sigma), 2)
		Q=1-P
		
	#	print den, sigma, OP, P

	MM1=lnc*M1+lne*(m1+e1)#-math.log(P**2)
	Mm1=lnch*(M1+m1)+lne*e1#-math.log(2*P*Q)
	mm1=lnc*m1+lne*(M1+e1)#-math.log(Q**2)

	MM2=lnc*M2+lne*(m2+e2)#-math.log(P**2)
	Mm2=lnch*(M2+m2)+lne*e2#-math.log(2*P*Q)
	mm2=lnc*m2+lne*(M2+e2)#-math.log(Q**2)

	norm1=math.log(math.exp(MM1)+math.exp(Mm1)+math.exp(mm1) )
	norm2=math.log(math.exp(MM2)+math.exp(Mm2)+math.exp(mm2) )

	File2.write( '\t'.join(map(str, [1, x, "C", "A", P, e] ) )+'\t')
	File.write( '\t'.join(map(str, [-MM1+norm1, -Mm1+norm1, -mm1+norm1, N1, -MM2+norm2, -Mm2+norm2, -mm2+norm2, N2, P] ) )+'\n')
	File2.write( '\t'.join(map(str, [-MM1+norm1, -Mm1+norm1, -mm1+norm1, N1, -MM2+norm2, -Mm2+norm2, -mm2+norm2, N2] ) )+'\n')
File.close()
dig=4
print  round(FA, dig), round(FC, dig), round(r, dig), round(sA, dig), round(sC, dig), round(z1, dig), round(z2, dig)
