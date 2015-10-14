import os
import math

A=[0,0,0,0,0,0,0,0,0]
B=[0,0,0,0,0,0,0,0,0]
N=0

for y in range(0, 10000):
	os.system("python simulate_random_IBD.py 3 10.00 10000 > sim/temp.txt")
	File=open("sim/temp.txt")
	E=map(float, File.read().split()[0:11])
	File.close()
	os.system("../bin/convert test2.txt sim/test2.dat > /dev/null")
	os.system("python relatedness.py --bcf sim/test2.dat -a 0 -A 0 -b 1 -B 1 > sim/temp2.txt")
	File=open("sim/temp2.txt")
	File.readline()
	O=[]
	line=File.readline().split(',')
	success=line[7]
	print success
	for x in range(8, len(line), 2):
		O.append(float(line[x]))
	if len(O)<7:
		continue
	File.close()
#	if float(E[7])>0.25:
#		continue
	print y, "all", 0-O[0], -E[1]+O[1], -E[0]+O[2], -E[2]+O[3], -E[4]+O[4], -E[3]+O[5], -E[5]+O[6], -E[6]+O[7]
	N+=1
	A[0]=A[0]*(N-1)/N-(0-O[0])/N
	A[1]=A[1]*(N-1)/N-(E[1]-O[1])/N
	A[2]=A[2]*(N-1)/N-(E[0]-O[2])/N
	A[3]=A[3]*(N-1)/N-(E[2]-O[3])/N
	A[4]=A[4]*(N-1)/N-(E[4]-O[4])/N
	A[5]=A[5]*(N-1)/N-(E[3]-O[5])/N
	A[6]=A[6]*(N-1)/N-(E[5]-O[6])/N
	A[7]=A[7]*(N-1)/N-(E[6]-O[7])/N
	A[8]=A[8]*(N-1)/N-(E[5]+E[6]-O[6]-O[7])/N

	B[0]=B[0]*(N-1)/N+(0-O[0])**2/N
	B[1]=B[1]*(N-1)/N+(E[1]-O[1])**2/N
	B[2]=B[2]*(N-1)/N+(E[0]-O[2])**2/N
	B[3]=B[3]*(N-1)/N+(E[2]-O[3])**2/N
	B[4]=B[4]*(N-1)/N+(E[4]-O[4])**2/N
	B[5]=B[5]*(N-1)/N+(E[3]-O[5])**2/N
	B[6]=B[6]*(N-1)/N+(E[5]-O[6])**2/N
	B[7]=B[7]*(N-1)/N+(E[6]-O[7])**2/N
	B[8]=B[8]*(N-1)/N+(E[6]+E[5]-O[6]-O[7])**2/N


	print y, "e", 0, O[0]
	print y, "fA", E[1], O[1]
	print y, "fC", E[0], O[2]
	print y, "r", E[2], O[3]
	print y, "sA", E[4], O[4]
	print y, "sC", E[3], O[5]
	print y, "z1", E[5], O[6]
	print y, "z2", E[6], O[7]

	print y, "running e, fA, fC, r, sAC, sCA, dAC, DAC, ZAC"
	print y, "running", ', '.join(map(str, (round(x*10**3, 2) for x in A ) ) )
	print y, "running", ', '.join(map(str, (round( x*10**4, 2) for x in B ) ) )
#	print y, "all-2", B
