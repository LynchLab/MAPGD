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
#ref_nuc=2
major_allele=2
minor_allele=3
pop_coverage=4
best_p=5
best_q=6
best_error=7
null_e=8
best_f=9
best_MM=10
best_Mm=11
best_mm=12
best_h=13
pol_llstat=14
HWE_llstat=15
Z_score=16

atoi={'A':0, 'C':1, 'G':2, 'T':3}

def likelihoods(calls, major, minor, error, p):
	error=0.01
	M=calls[major]
	m=calls[minor]
	E=sum(calls)-M-m
	lnc=math.log(1.0-error)
	lnch=math.log( (1.0-error)/2.0+error/6.0) 
	lne=math.log(error)
	MM=lnc*M+lne*(m+E)
	Mm=lnch*(M+m)+lne*E
	mm=lnc*m+lne*(M+E)
	N=math.log(0.5)*(M+m-1)
	if M+m>1:
		return [-MM, -Mm, -mm, sum(calls)]
	else:
		return [0, 0, 0, sum(calls)]
	
poly={}
for line in mapFile:
	line=line.split()
	if line[0]=="id1":
		continue
	if line[pol_llstat]=="*":
		continue
#	print line
#	print line[best_error], line[pol_llstat], line[pop_coverage]
	if float(line[best_error])<args.e:
		if float(line[pol_llstat])>args.l:
			if float(line[pop_coverage])>args.c and float(line[pop_coverage])<args.C:
				if float(line[best_q])<args.P and float(line[best_q])>args.p:
					print line[best_q]
					try:
						poly[line[scaffold]][line[site]]=[line[major_allele], line[minor_allele], line[best_p], line[best_error]]
					except:
						poly[line[scaffold]]={}
						poly[line[scaffold]][line[site]]=[line[major_allele], line[minor_allele], line[best_p], line[best_error]]

mapFile.close()

for line in proFile:
	line=line.split()
	try:
		this=poly[line[scaffold]][line[site]]
		out=[line[scaffold], line[site]]+this
		for x in range(3, len(line) ):
			calls=map(int, line[x].split('/') )
			out+=(likelihoods(calls, atoi[this[0]], atoi[this[1]], float(this[3]), float(this[2])))
		print '\t'.join(map(str, out) )

	except:
		C=0
proFile.close()
