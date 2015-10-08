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
parser.add_argument('--pro', metavar='proFile', type=str, nargs=1,
                   help='the name of a proFile',  required=True)
parser.add_argument('--map', metavar='mapFile', type=str, nargs=1,
                   help='the name of a mapFile',  required=True)
parser.add_argument('--mode', metavar='mode', type=str, nargs=1,
                   help=' <F|G|H|R|C>',  required=True)
args = parser.parse_args()

proFile=open(args.pro[0])
mapFile=open(args.map[0])

F=False
HET=False
GENOTYPE=False
COVERAGE=False
RELATEDNESS=False

if args.mode[0]=='F':
	F=True
if args.mode[0]=='H':
	HET=True
if args.mode[0]=='G':
	GENOTYPE=True
if args.mode[0]=='C':
	COVERAGE=True
if args.mode[0]=='R':
	RELATEDNESS=True
	z_n=0
	z_d=0
	r_n=0
	r_d=0
	FX_n=0
	FX_d=0
	FY_n=0
	FY_d=0
	X=0
	Y=6

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

atoi={'A':0, 'C':1, 'G':2, 'T':3, '.':5}

def likelihoods_uniform(calls, major, minor, error, p, pMM, pMm, pmm):
        error=max(0.001, error)
        M=calls[major]
	m=calls[minor]
        n=sum(calls)
	E=n-M-m

	lnc=math.log(1-error)
	lnch=math.log( (1-error)/2+error/6. )
	notc=math.log(error/3)

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
	lnch=math.log( (1-error)/2+error/6. )
	notc=math.log(error/3)

	MM=M*lnc+notc*(m+E)
	Mm=(M+m)*lnch+notc*(E)-math.log(1.01)
	mm=m*lnc+notc*(M+E)

        [E1, E2, E3]=sorted([MM, Mm, mm])
        N=math.log(math.exp(E1)+math.exp(E2)+math.exp(E3) )
	return [-MM+N, -Mm+N, -mm+N, n]

	
poly={}

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

mapFile.seek(0)

COV_SUM=0
COV_N=0

HEADER=True

for line in mapFile:
	line=line.split()
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
#		print line[pol_llstat]
		if float(line[pol_llstat])>=args.l:
			if float(line[best_q])<=args.P and float(line[best_q])>=args.p:
				if not (COVERAGE):
					if float(line[pop_coverage])>=args.c and float(line[pop_coverage])<=args.C:
						try:
							poly[line[scaffold]][line[site]]=[line[major_allele], line[minor_allele], line[best_p], line[best_error], line[best_MM], line[best_Mm], line[best_mm]]
						except:
							poly[line[scaffold]]={}
							poly[line[scaffold]][line[site]]=[line[major_allele], line[minor_allele], line[best_p], line[best_error], line[best_MM], line[best_Mm], line[best_mm]]
				else:
					COV_SUM+=int(line[pop_coverage])
					COV_N+=1
					print float(COV_SUM)/float(COV_N)
mapFile.close()


for line in proFile:
	line=line.split()
	if line[0][0]=="@":
		continue
	if (HEADER):
		if(GENOTYPE):
			out=[]
			for x in range(3, len(line) ):
				out.append(str(x-3)+"MM")
				out.append(str(x-3)+"Mm")
				out.append(str(x-3)+"mm")
			print "@scaffold\tbp\tmajor\tminor\tM\terror\t"+'\t'.join(map(str, out) )
		elif(HET or F):
			out=[]
			for x in range(3, len(line) ):
				out.append(x-3)
			print "@scaffold\tbp\tmajor\tminor\tM\terror\t"+'\t'.join(map(str, out) )
		HEADER=False
	try:
		this=poly[line[scaffold]][line[site]]
	except:
		continue
	out=[line[scaffold], line[site]]+this[0:4]
	if this[0]!='.' and this[1]!='.':
		for x in range(3, len(line) ):
			calls=map(int, line[x].split('/') )
			if (GENOTYPE):
				out+=(likelihoods_emperical(calls, atoi[this[0]], atoi[this[1]], float(this[3]), float(this[2]), float(this[4]), float(this[5]), float(this[6]) ))
			elif (RELATEDNESS):
				if x-3==X:
#					print line, this
					try:
						Lx=(likelihoods_emperical(calls, atoi[this[0]], atoi[this[1]], float(this[3]), float(this[2]), float(this[4]), float(this[5]), float(this[6]) ))
					except:
						print this, line
						quit()				
				elif x-3==Y:
					Mm=float(this[5])
					mm=float(this[6])
					MM=float(this[4])
					Ly=(likelihoods_emperical(calls, atoi[this[0]], atoi[this[1]], float(this[3]), float(this[2]), float(this[4]), float(this[5]), float(this[6]) ))
					SSxy=math.exp(-Lx[0])*math.exp(-Ly[0])+math.exp(-Lx[1])*math.exp(-Ly[1])+math.exp(-Lx[2])*math.exp(-Ly[2])
					ZZxy=math.exp(-Lx[0])*math.exp(-Ly[0])+math.exp(-Lx[1])*math.exp(-Ly[1])+math.exp(-Lx[2])*math.exp(-Ly[2])+math.exp(-Lx[0])*math.exp(-Lx[2])+math.exp(-Lx[2])*math.exp(-Ly[0])
					Mx=1./3.
					My=1./3.
					SSx=math.exp(-Lx[0])**2+math.exp(-Lx[1])**2+math.exp(-Lx[2])**2
					SSy=math.exp(-Ly[0])**2+math.exp(-Ly[1])**2+math.exp(-Ly[2])**2
					Vx=SSx-Mx**2
					Vy=SSy-My**2
					#print SSxy, ( ( SSxy-(Mm**2)-(mm**2)-(MM**2) ) )
					#FX_n+=?
					#FX_d+=?
					r_n+=( ( (1.+(SSxy-Mx*My)/math.sqrt(Vx*Vy) )/2-(Mm**2)-(mm**2)-(MM**2) )/(1-Mm**2-mm**2-MM**2 ) )*math.sqrt(Vx*Vy)
					r_d+=math.sqrt(Vx*Vy)*3.
					if Mm!=0:
						z_n+=( ( (1.+(ZZxy-Mx*My)/math.sqrt(Vx*Vy) )/2.-(Mm**2)-(mm**2)-(MM**2)-2*MM*mm )/(1-Mm**2-mm**2-MM**2-MM*mm*2 ) )*math.sqrt(Vx*Vy)
						z_d+=math.sqrt(Vx*Vy)*3.
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
		if (not (RELATEDNESS) ):
			print '\t'.join(map(str, out) )
if (RELATEDNESS):
	print r_n/r_d
	print z_n/z_d
	print FY_n/FY_d
proFile.close()

