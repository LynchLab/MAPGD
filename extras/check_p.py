#This compares the allele frequency reported in the map file to the allele frequency obtained from running mapgdutils.py. They should be nearly identical. Ideally they should
#be exactally identical, but since I haven't implemented a real numeric estimation method yet the allele frequencies aren't as good as they should be.

#!/bin/python

import math
import sys
File=open(sys.argv[1])
for line in File:
	line=line.split()
	P=float(line[4])
	ll=map(float, line[6:])
	S=0
	for x in range(0, len(ll), 4):
		S+=math.exp(-ll[x])*2+math.exp(-ll[x+1])
	print S, P*98*2
