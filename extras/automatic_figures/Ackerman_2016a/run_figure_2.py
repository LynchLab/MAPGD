#!/usr/bin/python
import os
import math

B=[[],[],[],[],[],[]]
A=[[],[],[],[],[],[]]
N=0

n_loci=[5000, 100000]
cov=[3, 10]

cov=3
n_loci=5000	
for x in range (0, 6):	
	A[x]=[0,0,0,0,0,0,0,0,0]
	B[x]=[0,0,0,0,0,0,0,0,0]

for y in range(0, 1000):
	os.system("python simulate_random_IBD.py -c "+str(cov)+" -s 0.00 -N "+str(n_loci)+" > sim/temp.txt")
	File=open("sim/temp.txt")
	E=map(float, File.read().split()[0:11])
	File.close()

	N+=1
#	os.system("cat test2.txt | ~/src/mapgd/bin/mapgd filter -q "+str(0.01)+"> sim/test2.txt")
	os.system("cat test2.txt | ~/src/mapgd/bin/mapgd relatedness > sim/temp2.txt")
	File=open("sim/temp2.txt")
	File.readline()
	O=[]
	line=File.readline().split('\t')
	line=File.readline().strip('\t\n').split('\t')
#	print line
	O.append(float(0))
	for z in range(2, len(line), 2):
		O.append(float(line[z]))
	if len(O)<7:
		continue
	File.close()

	if( E[1]!=0):
		print 1, E[1], O[1]
	if (E[0]!=0):
		print 2, E[0], O[2]
	if (E[2]!=0):
		print 3, E[2], O[3]
	if (E[4]!=0):
		print 4, E[4], O[4]
	if (E[3]!=0):
		print 5, E[3], O[5]
	if (E[5]!=0):
		print 6, E[5], O[6]
	if (E[6]!=0):
		print 7, E[6], O[7]
