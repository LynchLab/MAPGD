#!/usr/bin/python
import sys
import argparse

parser = argparse.ArgumentParser(description='scans a gcf file for markers')
parser.add_argument('genotypic.tsv', metavar='genotypic.tsv', type=str, nargs=1,
                   help='the name of a tab seperated genotypic liklelihood file')
args = parser.parse_args()

proFile=open(args.proFile[0])
mapFile=open(args.mapFile[0])


sites={}
def comp (c, a):
	ret={'W':{'A':'T', 'T':'A'}, 'S':{'G':'C', 'C':'G'}, 'K':{'G':'T', 'T':'G'}, 'M':{'A':'C', 'C':'A'}, 'Y':{'C':'T', 'T':'C'}, 'R':{'A':'G', 'G':'A'} }
	try:
		return ret[a][c]
	except:
		print c, a
		return ret[a][c]

def toi (c):
	ret={'A':0, 'C':1, 'G':2, 'T':3}
	return ret[c]

markerFile=open(tsvFile)

for line in markerFile:
	line=line.strip('\n').split()
	scaffold=line[0].split('_')[1]
	bp=line[1]
	try:
		sites[scaffold][bp]=comp(line[2], line[3])
	except:
		sites[scaffold]={}
		try:
			sites[scaffold][bp]=comp(line[2], line[3])
		except:
			print line	

print sites.keys()
processed=0
File=open(gfcFile)
num=[]
denom=[]

for x in range(0, 96):
	num.append(0)
	denom.append(0)

for line in File:
	line=line.split()
	scaffold=line[0]
	bp=line[1]
	try:
#		print scaffold, bp
		asex=sites[scaffold][bp]
		MAJOR=line[2]
		MINOR=line[3]
		print scaffold, bp, asex, MAJOR, MINOR
		calls=map(float, line[6:] )
		for y in range(0, 96):
			lcal=calls[y*4:(y+1)*4]
			if lcal[3]>4:
				if lcal[1]<lcal[0] and lcal[1]<lcal[2]:
					num[y]+=1
				denom[y]+=1
	except:
		processed+=1
for x in range(0, 96):
	if denom[x]>0:
		print x, float(num[x])/float(denom[x]), denom[x]
		
