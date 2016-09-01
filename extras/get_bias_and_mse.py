import os
import math

B=[[],[],[],[],[],[]]
A=[[],[],[],[],[],[]]
N=0
		
for x in range (0, 6):	
	A[x]=[0,0,0,0,0,0,0,0,0]
	B[x]=[0,0,0,0,0,0,0,0,0]

for y in range(0, 1000):
	os.system("python simulate_random_IBD.py -c 3 -s 0.00 -N 100000 > sim/temp.txt")
	File=open("sim/temp.txt")
	E=map(float, File.read().split()[0:11])
	File.close()

	N+=1
	for x in range (0, 6):	
		if (x==0):
			os.system("cat test2.txt | python filter.py "+str(0.99)+"> sim/test2.txt")
		else:
			os.system("cat test2.txt | python filter.py "+str(1.-0.05*x)+"> sim/test2.txt")
		os.system("cat sim/test2.txt | ../bin/mapgd relatedness > sim/temp2.txt")
		File=open("sim/temp2.txt")
		File.readline()
		O=[]
		line=File.readline().split('\t')
		line=File.readline().strip('\t\n').split('\t')

		O.append(float(0))
		for z in range(2, len(line), 2):
			O.append(float(line[z]))
		if len(O)<7:
			continue
		File.close()

		A[x][0]=A[x][0]*(N-1)/N-(0-O[0])/N
		A[x][1]=A[x][1]*(N-1)/N-(E[1]-O[1])/N
		A[x][2]=A[x][2]*(N-1)/N-(E[0]-O[2])/N
		A[x][3]=A[x][3]*(N-1)/N-(E[2]-O[3])/N
		A[x][4]=A[x][4]*(N-1)/N-(E[4]-O[4])/N
		A[x][5]=A[x][5]*(N-1)/N-(E[3]-O[5])/N
		A[x][6]=A[x][6]*(N-1)/N-(E[5]-O[6])/N
		A[x][7]=A[x][7]*(N-1)/N-(E[6]-O[7])/N
		A[x][8]=A[x][8]*(N-1)/N-(E[5]+E[6]-O[6]-O[7])/N

		B[x][0]=B[x][0]*(N-1)/N+(0-O[0])**2/N
		B[x][1]=B[x][1]*(N-1)/N+(E[1]-O[1])**2/N
		B[x][2]=B[x][2]*(N-1)/N+(E[0]-O[2])**2/N
		B[x][3]=B[x][3]*(N-1)/N+(E[2]-O[3])**2/N
		B[x][4]=B[x][4]*(N-1)/N+(E[4]-O[4])**2/N
		B[x][5]=B[x][5]*(N-1)/N+(E[3]-O[5])**2/N
		B[x][6]=B[x][6]*(N-1)/N+(E[5]-O[6])**2/N
		B[x][7]=B[x][7]*(N-1)/N+(E[6]-O[7])**2/N
		B[x][8]=B[x][8]*(N-1)/N+(E[6]+E[5]-O[6]-O[7])**2/N
	

	print y, "running e, fA, fC, r, sAC, sCA, dAC, DAC, ZAC"

	print y, "MSE  (x10^4) P>0.01", ', '.join(map(str, (round(x*10**4, 2) for x in B[0] ) ) )
	print y, "MSE  (x10^4) P>0.05", ', '.join(map(str, (round(x*10**4, 2) for x in B[1] ) ) )
	print y, "MSE  (x10^4) P>0.10:", ', '.join(map(str, (round(x*10**4, 2) for x in B[2] ) ) )
	print y, "MSE  (x10^4) P>0.15:", ', '.join(map(str, (round(x*10**4, 2) for x in B[3] ) ) )
	print y, "MSE  (x10^4) P>0.20:", ', '.join(map(str, (round(x*10**4, 2) for x in B[4] ) ) )
	print y, "MSE  (x10^4) P>0.25:", ', '.join(map(str, (round(x*10**4, 2) for x in B[5] ) ) )

	print y, "Bias (x10^3) P>0.01", ', '.join(map(str, (round(x*10**3, 2) for x in A[0] ) ) )
	print y, "Bias (x10^3) P>0.05", ', '.join(map(str, (round(x*10**3, 2) for x in A[1] ) ) )
	print y, "Bias (x10^3) P>0.10:", ', '.join(map(str, (round(x*10**3, 2) for x in A[2] ) ) )
	print y, "Bias (x10^3) P>0.15:", ', '.join(map(str, (round(x*10**3, 2) for x in A[3] ) ) )
	print y, "Bias (x10^3) P>0.20:", ', '.join(map(str, (round(x*10**3, 2) for x in A[4] ) ) )
	print y, "Bias (x10^3) P>0.25:", ', '.join(map(str, (round(x*10**3, 2) for x in A[5] ) ) )
#	print y, "all-2", B
