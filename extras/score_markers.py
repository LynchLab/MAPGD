#!/usr/bin/python
import sys
import argparse

parser = argparse.ArgumentParser(description='scans a gcf file for markers')
parser.add_argument('gcfFile', metavar='gcfFile', type=str, nargs=1,
                   help='the name of a genome file')
parser.add_argument('tsvFile', metavar='tsvFile', type=str, nargs=1,
                   help='the name of a tab seperated file containing markers')
args = parser.parse_args()

sites={}
def comp (c, a):
	ret={'W':{'A':'T', 'T':'A'}, 'S':{'G':'C', 'C':'G'}, 'K':{'G':'T', 'T':'G'}, 'M':{'A':'C', 'C':'A'}, 'Y':{'C':'T', 'T':'C'}, 'R':{'A':'G', 'G':'A'} }
	try:
		return ret[a][c]
	except:
		if a not in ret.keys():
			print a, " is WSKY or R"
			return None
		else:
			print c, "is not a base in ", a
			return None

def toi (c):
	ret={'A':0, 'C':1, 'G':2, 'T':3}
	return ret[c]

markerFile=open(args.tsvFile[0])

G=0
E=0

for line in markerFile:
	if line[0]=='@':
		continue
	line=line.strip('\n').split()
	scaffold=line[0].split('_')[1]
	bp=line[1]
	K=comp(line[2], line[3])
	if K!=None:
		G+=1
	else:
		E+=1
	try:
		sites[scaffold][bp]=K
	except:
		sites[scaffold]={}
		sites[scaffold][bp]=K
print G

processed=0
File=open(args.gcfFile[0])
num=[]
denom=[]
Pulex=[]
Pulicaria=[]

for x in range(0, 96):
	num.append(0)
	denom.append(0)
	Pulex.append(0)
	Pulicaria.append(0)

print sites.keys()

for line in File:
	line=line.split()
	scaffold=str(line[0])
	bp=str(line[1])
	try:
		asex=sites[scaffold][bp]
		MAJOR=line[2]
		MINOR=line[3]
		calls=map(float, line[6:] )
		for y in range(0, 96):
			lcal=calls[y*4:(y+1)*4]
			if lcal[3]>4:
				if lcal[1]<lcal[0] and lcal[1]<lcal[2]:
					num[y]+=1
				if lcal[0]<lcal[1] and lcal[0]<lcal[2]:
					if MAJOR!=asex:
						Pulex[y]+=1
					else:
						Pulicaria[y]+=1
				if lcal[2]<lcal[0] and lcal[2]<lcal[1]:
					if MAJOR==asex:
						Pulex[y]+=1
					else:
						Pulicaria[y]+=1
				
				denom[y]+=1
	except:
		processed+=1
print processed
for x in range(0, 96):
	if denom[x]>0:
		print x, float(num[x])/float(denom[x]), float(Pulex[x])/float(denom[x]), float(Pulicaria[x])/float(denom[x]), denom[x]
