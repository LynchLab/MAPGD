import random
import sys


def randbase (di):
	return di[random.randint(0, 1)]

File=open(sys.argv[1])
header=File.readline().split()

Major=2
Minor=3
Start=6
N=(len(header)-6)/3

seq=[]
for x in range(0, N):
	seq.append([])

for line in File:
	line=line.strip('\n').split()
	di=[line[Major], line[Minor] ]
	for x in range(0, N):
		G=map(float, line[x*4+Start:x*4+3+Start])
		if min(G)<0.1:
			I=G.index(min(G))
			if I==0:
				seq[x].append(di[0])
			elif I==1:
				seq[x].append(randbase(di) )
			else:
				seq[x].append(di[1])
		else:
			seq[x].append('N')
for x in range(0, N):
	print '>', x
	print ''.join(seq[x])
