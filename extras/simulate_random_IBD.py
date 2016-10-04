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
parser.add_argument('-s', metavar='--sigma', type=float, default=0.0,
                   help='simulated sampling error')
parser.add_argument('-r', metavar='--rsq', type=float, default=0.0,
                   help='ld between adjacent SNPs')
parser.add_argument('-l', metavar='--l2i', type=bool, default=False,
                   help='include the observed individual in allele frequency estimates')
parser.add_argument('-I', metavar='--inbred', type=bool, default=False,
                   help='include the observed individual in allele frequency estimates')
parser.add_argument('-M', metavar='--cM', type=float, default=400,
                   help='genome size in centiMorgans.')
args = parser.parse_args()



def getMs (PA, PC, S):
	P=PA
	P=PC

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
		mmmm, Mmmm, MMmm, mmMm, MmMm, MMMm, mmMM, MmMM, MMMM=getMs(float(x)/10., float(x)/10., S)
		if ( min([mmmm, Mmmm, MMmm, mmMm, MmMm, MMMm, mmMM, MmMM, MMMM])<-0.0001):
			return False
	return True

#def get_conditional(X, rsq):
	

Fa=0
Fc=1
r=2
sac=3
Sca=4
z1=5
z2=6

Z= [0,1,2,3,4,5,6]

S=[random.random()-0.5, random.random()-0.5, random.random()-0.5, random.random()-0.5, random.random()-0.5, random.random()-0.5, random.random()-0.5 ]

#S=[0,0, 0.10, 0, 0, 0, 0]

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
l2i=False

MAX=args.N
N1=args.c
N2=args.c
sigma=args.s
ld=args.r
l2i=args.l
inbred=args.I
genome_size=args.M

MINP=0.005
MAXP=0.5

while ( not ( test(S) ) ):
	S=[random.random()-0.5, random.random()-0.5, random.random()-0.5, random.random()-0.5, random.random()-0.5, random.random()-0.5, random.random()-0.5 ]
	for x in range(0, random.randint(0, 6) ):
		S[Z[x]]=0

if (inbred):
	S=[0.5, 0.5, 0.5, 0.25, 0.25, 0.125, 0.375]
else:
	S=[0.0, 0.0, 0.25, 0.0, 0.0, 0.0, 0.25]
FC, FA, r, sC, sA, z1, z2=S



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

AMOM=random.randint(0,1)
ADAD=random.randint(0,1)
CMOM=random.randint(0,1)
CDAD=random.randint(0,1)

chromosomes=14
cMperSNP=min(0.5, float(genome_size)/float(MAX)/100.*chromosomes)
#print cMperSNP
markers_per_chromosome=MAX/float(chromosomes)
NEXT_CHROMOSOME=0
count=0

Theta=0
Delta=0
GammaAC=0
GammaCA=0
delta=0
fA=0
fC=0

while x<MAX+skipped:

		
	x+=1	

	P=0

	while( (P<=MINP) or (P>MAXP) ):
#		P = (numpy.random.pareto(0.1, 1)[0]+1 )/1000
#		print P
#		P = round(triang.rvs(0.1), 2)
		P = round(random.random(), 2)/2.

	Q=1-P

	MOM=[int(random.random()<P),int(random.random()<P)]
	if not(inbred):
		DAD=[int(random.random()<P),int(random.random()<P)]
	else:
		DAD=MOM


	if (x>NEXT_CHROMOSOME):
		AMOM=random.randint(0,1)
		ADAD=random.randint(0,1)
		CMOM=random.randint(0,1)
		CDAD=random.randint(0,1)
#		print AMOM, ADAD, CMOM, CDAD
		NEXT_CHROMOSOME+=markers_per_chromosome
	else:
		if(random.random()<cMperSNP):
			AMOM=int(AMOM==0)
		if(random.random()<cMperSNP):
			ADAD=int(ADAD==0)
		if(random.random()<cMperSNP):
			CMOM=int(CMOM==0)
		if(random.random()<cMperSNP):
			CDAD=int(CDAD==0)
	if (inbred):
		delta+=(float(AMOM==CMOM and ADAD==CDAD and AMOM==ADAD) )

		Theta+=(float(AMOM==CMOM) )/4.
		Theta+=(float(AMOM==CDAD) )/4.
		Theta+=(float(ADAD==CMOM) )/4.
		Theta+=(float(ADAD==CDAD) )/4.

		Delta+=float(AMOM==CMOM and ADAD==CDAD and AMOM!=ADAD)
		Delta+=float(AMOM==CDAD and ADAD==CMOM and AMOM!=ADAD)
		Delta+=float(AMOM==ADAD and CDAD==CMOM and AMOM!=CDAD)
		fA+=float(AMOM==ADAD)
		fC+=float(CMOM==CDAD)
		GammaAC+=(float(AMOM==ADAD and (ADAD==CDAD) ) )/2.
		GammaAC+=(float(AMOM==ADAD and (ADAD==CMOM) ) )/2.
		GammaCA+=(float(CMOM==CDAD and (ADAD==CDAD) ) )/2.
		GammaCA+=(float(CMOM==CDAD and (CDAD==AMOM) ) )/2.
	else:
		Theta+=(float(AMOM==CMOM)+float(ADAD==CDAD) )/4.
		Delta+=(float(AMOM==CMOM and ADAD==CDAD) )

	A=MOM[AMOM]+DAD[ADAD]
	C=MOM[CMOM]+DAD[CDAD]

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
print fC/float(MAX), fA/float(MAX), Theta/float(MAX), GammaCA/float(MAX), GammaAC/float(MAX), delta/float(MAX), Delta/float(MAX)
