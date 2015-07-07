#!/usr/bin/python
import sys
import argparse
import math
import argparse

parser = argparse.ArgumentParser(description='Utilities to help with the analysis of the output of mapgd')
parser.add_argument('proFile', metavar='proFile', type=str, nargs=1,
                   help='the name of a proFile')
parser.add_argument('mapFile', metavar='mapFile', type=str, nargs=1,
                   help='the name of a mapFile')
parser.add_argument('-p', metavar='--maf', type=float, default=0.0,
                   help='minimum allele frequency for analysis')
parser.add_argument('-P', metavar='--max-maf', type=float, default=0.5, 
                   help='maximum allele frequency for analysis')
parser.add_argument('-z', metavar='--zmin', type=float, default=0.0,
                   help='minimum z-score for a site')
parser.add_argument('-l', metavar='--minll', type=float, default=0.0,
                   help='minimum ll-score for a site (analogous to minQ)')
parser.add_argument('-c', metavar='--cmin', type=int, default=0,
                   help='minimum population coverage for analysis')
parser.add_argument('-C', metavar='--cmax', type=int, default=sys.maxint,
                   help='maximum population coverage for analysis')
parser.add_argument('-e', metavar='--error', type=float, default=0.1, 
                   help='maximum estimated sequencing error')
parser.add_argument('-H', metavar='--hwe', type=float, default=0.0, 
                   help='maximum ll-score for hw diseqaulibrium')
args = parser.parse_args()

proFile=open(args.proFile[0])
mapFile=open(args.mapFile[0])

#0	1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	
#id1	id2	ref	major	minor	cov	M	m	error	null_e	f	MM	Mm	mm	h	polyll	HWEll	gof	eff_chrom	N	N_excluded	model_ll

scaffold=0
site=1	
ref_nuc=2
major_allele=3
minor_allele=4
pop_coverage=5
best_p=6
best_q=7
best_error=8
null_e=9
best_f=10
best_MM=11
best_Mm=12
best_mm=13
best_h=14
pol_llstat=15
HWE_llstat=16
Z_score=17

atoi={'A':0, 'C':1, 'G':2, 'T':3}

def likelihoods_uniform(calls, major, minor, error, p):
        error=max(0.001, error)
        M=calls[major]
        n=sum(calls)
        p2=math.log(1.0-error)
        notp2=math.log(error)
        p1=math.log( (1.0-error)/2.0+error/6.0)
        notp1=math.log(1-( (1.0-error)/2.0+error/6.0) )
        p0=math.log(error/3.)
        notp0=math.log(1-error/3.)

        MM=M*p2+notp2*(n-M)+math.log(p**2)
        Mm=M*p1+notp1*(n-M)+math.log(2*p*(1-p) )
        mm=M*p0+notp0*(n-M)+math.log( (1-p)**2 )
        [E1, E2, E3]=sorted([MM, Mm, mm])
        N=math.log(math.exp(E1)+math.exp(E2)+math.exp(E3) )
        if n>=1:
                return [-MM+N, -Mm+N, -mm+N, n, '|' ]
        else:
                return [0, 0, 0, sum(calls), '|' ]

def likelihoods_emperical(calls, major, minor, error, p, pMM, pMm, pmm):
#        error=max(0.001, error)
        M=calls[major]
        n=sum(calls)
        p2=math.log(1.0-error)
        notp2=math.log(error)
        p1=math.log( (1.0-error)/2.0+error/6.0)
        notp1=math.log(1-( (1.0-error)/2.0+error/6.0) )
        p0=math.log(error/3.)
        notp0=math.log(1-error/3.)

        MM=M*p2+notp2*(n-M)+math.log(pMM)
        Mm=M*p1+notp1*(n-M)+math.log(pMm)+math.log(0.5)
        mm=M*p0+notp0*(n-M)+math.log(pmm)
        [E1, E2, E3]=sorted([MM, Mm, mm])
        N=math.log(math.exp(E1)+math.exp(E2)+math.exp(E3) )
	return [-MM+N, -Mm+N, -mm+N, n]

	
poly={}
for line in mapFile:
	line=line.split()
	if line[0]=="id1":
		continue
	if line[pol_llstat]=="*":
		continue
	if float(line[best_error])<=args.e:
		if float(line[pol_llstat])>=args.l:
			if float(line[pop_coverage])>=args.c and float(line[pop_coverage])<=args.C:
				if float(line[best_q])<=args.P and float(line[best_q])>=args.p:
					try:
						poly[line[scaffold]][line[site]]=[line[major_allele], line[minor_allele], line[best_p], line[best_error], line[best_MM], line[best_Mm], line[best_mm]]
					except:
						poly[line[scaffold]]={}
						poly[line[scaffold]][line[site]]=[line[major_allele], line[minor_allele], line[best_p], line[best_error], line[best_MM], line[best_Mm], line[best_mm]]

mapFile.close()

for line in proFile:
	line=line.split()
	try:
		this=poly[line[scaffold]][line[site]]
		out=[line[scaffold], line[site]]+this[0:4]
		for x in range(3, len(line) ):
			calls=map(int, line[x].split('/') )
			out+=(likelihoods_emperical(calls, atoi[this[0]], atoi[this[1]], float(this[3]), float(this[2]), float(this[4]), float(this[5]), float(this[6]) ))
		print '\t'.join(map(str, out) )

	except:
		C=0
proFile.close()
