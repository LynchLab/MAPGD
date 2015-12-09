#!/usr/bin/python
import sys
import argparse
import math
import argparse

parser = argparse.ArgumentParser(description='Utilities to help with the analysis of the output of mapgd')
parser.add_argument('-p', metavar='--maf', type=float, default=0.0,
                   help='minimum allele frequency for analysis')
parser.add_argument('-P', metavar='--max-maf', type=float, default=0.5, 
                   help='maximum allele frequency for analysis')
parser.add_argument('-z', metavar='--zmin', type=float, default=-sys.float_info.max,
                   help='minimum z-score for a site')
parser.add_argument('-N', metavar='--nmax', type=float, default=0.0,
                   help='maximum number of excluded sites')
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
parser.add_argument('--pro', metavar='proFile', type=str, nargs=1,
                   help='the name of a proFile',  required=True)
parser.add_argument('--map', metavar='mapFile', type=str, nargs=1,
                   help='the name of a mapFile',  required=True)
parser.add_argument('--mode', metavar='mode', type=str, nargs=1,
                   help=' <F|G|H>',  required=True)
args = parser.parse_args()

proFile=open(args.pro[0])
mapFile=open(args.map[0])

F=False
HET=False
GENOTYPE=False
COVERAGE=False

if args.mode[0]=='F':
	F=True
if args.mode[0]=='H':
	HET=True
if args.mode[0]=='G':
	GENOTYPE=True
if args.mode[0]=='C':
	COVERAGE=True

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
N_included=19
N_excluded=20


atoi={'A':0, 'C':1, 'G':2, 'T':3, '.':5}

def likelihoods_uniform(calls, major, minor, error, p, pMM, pMm, pmm):
        error=max(0.001, error)
        M=calls[major]
	m=calls[minor]
        n=sum(calls)
	E=n-M-m

	lnc=math.log(1-error)
	lnch=math.log( (1-error)/2.+error/6. )
	notc=math.log(error/3.)

	MM=M*lnc+notc*(m+E)-math.log(p**2)
	Mm=(M+m)*lnch+notc*(E)-math.log(2*p*(1-p ) )
	mm=m*lnc+notc*(M+E)-math.log( (1-p)**2 )

        [E1, E2, E3]=sorted([MM, Mm, mm])
        N=math.log(math.exp(E1)+math.exp(E2)+math.exp(E3) )
        if n>=1:
                return [-MM+N, -Mm+N, -mm+N, n]
        else:
                return [0, 0, 0, sum(calls)]

def likelihoods_emperical(calls, major, minor, error, p, pMM, pMm, pmm):
#        error=0.005#max(prior, error)
        error=max(0.001, error)
        M=calls[major]
	m=calls[minor]
	n=sum(calls)
	E=n-M-m

	lnc=math.log(1-error)
	lnch=math.log( (1-error)/2.+error/6. )
	notc=math.log(error/3.)

	MM=M*lnc+notc*(m+E)
	Mm=(M+m)*lnch+notc*(E)-math.log(1.01)
	mm=m*lnc+notc*(M+E)

        [E1, E2, E3]=sorted([MM, Mm, mm])
	try:
	        N=math.log(math.exp(E1)+math.exp(E2)+math.exp(E3) )
		return [-MM+N, -Mm+N, -mm+N, n]
	except:
		return [1/3,1/3,1/3, 0]

	
name={}

#E=0
#N=0
#for line in mapFile:
#	line=line.split()
#        try:
#		if line[0]=="scaffold" or "id1":
#                        continue
#                if line[0][0]=='@':
#                        continue
#		if line[pol_llstat]=="." or line[pol_llstat]=="*" or line[pol_llstat]=="NA":
#                        continue
#		if line[best_error]=="." or line[best_error]=="*" or line[best_error]=="NA":
#                        continue
#        except:
#		print "could not parse line ", line
#		exit(0)
#	try:
#		E+=float(line[best_error])
#	except:
#		print line
#		exit(0)
#	N+=1
#
#prior=(E/N)

COV_SUM=0
COV_N=0
SUM_N=0

HEADER=True

#for mapline in mapFile:
while(True):
	while(True):
		line=mapFile.readline().split()
		try:
			if line[0]=="scaffold" or line[0]=="id1":
				continue
			if line[0][0]=='@':
				continue
			if line[pol_llstat]=="." or line[pol_llstat]=="*" or line[pol_llstat]=="NA":
				continue
			if line[best_error]=="." or line[best_error]=="*" or line[best_error]=="NA":
				continue
		except:
			print "could not parse line ", line
			exit(0)
		if float(line[best_error])<=args.e:
			if float(line[pol_llstat])>=args.l:
				if float(line[best_q])<=args.P and float(line[best_q])>=args.p:
					if not (COVERAGE):
						if float(line[N_excluded])<=args.N:
							if float(line[Z_score])>=args.z:
								if float(line[pop_coverage])>=args.c and float(line[pop_coverage])<=args.C:
									this=[line[major_allele], line[minor_allele], line[best_p], line[best_error], line[best_MM], line[best_Mm], line[best_mm]]
									next_scaffold=line[scaffold]
									next_site=line[site]
								break
					else:
						COV_SUM+=float(line[pop_coverage])
						COV_N+=1
						SUM_N+=float(line[N_included])
						if COV_N>0 and SUM_N>0:
							AVG_N=float(SUM_N)/float(COV_N) 
							print round(float(COV_SUM)/float(COV_N), 2), round(float(COV_SUM)/float(COV_N)/AVG_N, 2)
	while (True):
		line=proFile.readline().split()
		if line[0][0]=="@":
			if line[0]=="@ID":
				for x in range(4, len(line) ):
					name[str(x-4)]=line[x]
			continue
		if (HEADER):
			if(GENOTYPE):
				out=[]
				for x in range(3, len(line) ):
					out.append(name[str(x-3)]+"_MM")
					out.append(name[str(x-3)]+"_Mm")
					out.append(name[str(x-3)]+"_mm")
					out.append(name[str(x-3)]+"_N")
				print "@scaffold\tbp\tmajor\tminor\tM\terror\t"+'\t'.join(map(str, out) )
			elif(HET or F):
				out=[]
				for x in range(3, len(line) ):
					out.append(x-3)
				print "@scaffold\tbp\tmajor\tminor\tM\terror\t"+'\t'.join(map(str, out) )
			HEADER=False
		if next_scaffold!=line[scaffold] or next_site!=line[site]:
			continue
		out=[next_scaffold, next_site]+this[0:4]
		if this[0]!='.' and this[1]!='.' and this[0]!='*' and this[1]!='*':
			for x in range(3, len(line) ):
				calls=map(int, line[x].split('/') )
				if (GENOTYPE):
					out+=(likelihoods_emperical(calls, atoi[this[0]], atoi[this[1]], float(this[3]), float(this[2]), float(this[4]), float(this[5]), float(this[6]) ))
				elif (HET):
					M=atoi[this[0]]
					m=atoi[this[1]]
					if (calls[M]+calls[m]>0):
						out.append(float(calls[M])/(float(calls[M])+float(calls[m]) ) )
					else:
						out.append('.')
				elif (F):
					Mm=float(this[5])
					mm=float(this[6])
					MM=float(this[4])
					p=float(this[2])
					if (p!=0 and p!=1):
						f=1-Mm/(2*(1-p)*p )
						out.append(f)
					else:
						out.append('.')	
					#f=1-(f-1)
			print '\t'.join(map(str, out) )
		break

