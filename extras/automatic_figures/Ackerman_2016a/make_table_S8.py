import os
import math

B=[[],[],[]]
A=[[],[],[]]
B2=[[],[],[]]
A2=[[],[],[]]
BI=[[],[],[]]
AI=[[],[],[]]
BI2=[[],[],[]]
AI2=[[],[],[]]
N=0
		
for i in range (0, 3):
	A[i]=[[],[],[]]
	B[i]=[[],[],[]]
	B2[i]=[[],[],[]]
	A2[i]=[[],[],[]]
	AI[i]=[[],[],[]]
	BI[i]=[[],[],[]]
	BI2[i]=[[],[],[]]
	AI2[i]=[[],[],[]]
	for j in range (0, 3):	
		A[i][j]=[0,0,0,0,0,0,0,0,0]
		B[i][j]=[0,0,0,0,0,0,0,0,0]
		B2[i][j]=[0,0,0,0,0,0,0,0,0]
		A2[i][j]=[0,0,0,0,0,0,0,0,0]
		AI[i][j]=[0,0,0,0,0,0,0,0,0]
		BI[i][j]=[0,0,0,0,0,0,0,0,0]
		BI2[i][j]=[0,0,0,0,0,0,0,0,0]
		AI2[i][j]=[0,0,0,0,0,0,0,0,0]

CM=[1000000]

NUMBER_REP=500
NUMBER_LOCI=10000

COV=[3,3,3]
sigma=["-s 0.0", "-s 188.0 -l True", "-s 94.0 -l True"]
print1=["0.00", "0.05", "0.10"]
print2=["$infty$", "96", "48"]

for i in range (0, 3):
	for j in range (0, 3):
		N=0
		for y in range(0, NUMBER_REP):
			os.system("python simulate_random_IBD.py -c "+str(COV[i])+" "+sigma[i]+" -N "+str(NUMBER_LOCI)+" -M "+str(CM[0])+" > sim/temp.txt")
			File=open("sim/temp.txt")
			E=map(float, File.readline().split()[0:11])
			E2=map(float, File.readline().split()[0:11])
			File.close()

			N+=1
			os.system("cat test2.txt | python filter.py "+str(1.0-0.05*j)+"> sim/test2.txt")
			os.system("cat sim/test2.txt | mapgd relatedness > sim/temp2.txt")
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

			A[i][j][0]=A[i][j][0]*(N-1)/N-(0-O[0])/N
			A[i][j][1]=A[i][j][1]*(N-1)/N-(E[1]-O[1])/N
			A[i][j][2]=A[i][j][2]*(N-1)/N-(E[0]-O[2])/N
			A[i][j][3]=A[i][j][3]*(N-1)/N-(E[2]-O[3])/N
			A[i][j][4]=A[i][j][4]*(N-1)/N-(E[4]-O[4])/N
			A[i][j][5]=A[i][j][5]*(N-1)/N-(E[3]-O[5])/N
			A[i][j][6]=A[i][j][6]*(N-1)/N-(E[5]-O[6])/N
			A[i][j][7]=A[i][j][7]*(N-1)/N-(E[6]-O[7])/N
			A[i][j][8]=A[i][j][8]*(N-1)/N-(E[5]+E[6]-O[6]-O[7])/N

			B[i][j][0]=B[i][j][0]*(N-1)/N+(0-O[0])**2/N
			B[i][j][1]=B[i][j][1]*(N-1)/N+(E[1]-O[1])**2/N
			B[i][j][2]=B[i][j][2]*(N-1)/N+(E[0]-O[2])**2/N
			B[i][j][3]=B[i][j][3]*(N-1)/N+(E[2]-O[3])**2/N
			B[i][j][4]=B[i][j][4]*(N-1)/N+(E[4]-O[4])**2/N
			B[i][j][5]=B[i][j][5]*(N-1)/N+(E[3]-O[5])**2/N
			B[i][j][6]=B[i][j][6]*(N-1)/N+(E[5]-O[6])**2/N
			B[i][j][7]=B[i][j][7]*(N-1)/N+(E[6]-O[7])**2/N
			B[i][j][8]=B[i][j][8]*(N-1)/N+(E[6]+E[5]-O[6]-O[7])**2/N

			A2[i][j][0]=A2[i][j][0]*(N-1)/N-(0-O[0])/N
			A2[i][j][1]=A2[i][j][1]*(N-1)/N-(E2[1]-O[1])/N
			A2[i][j][2]=A2[i][j][2]*(N-1)/N-(E2[0]-O[2])/N
			A2[i][j][3]=A2[i][j][3]*(N-1)/N-(E2[2]-O[3])/N
			A2[i][j][4]=A2[i][j][4]*(N-1)/N-(E2[4]-O[4])/N
			A2[i][j][5]=A2[i][j][5]*(N-1)/N-(E2[3]-O[5])/N
			A2[i][j][6]=A2[i][j][6]*(N-1)/N-(E2[5]-O[6])/N
			A2[i][j][7]=A2[i][j][7]*(N-1)/N-(E2[6]-O[7])/N
			A2[i][j][8]=A2[i][j][8]*(N-1)/N-(E2[5]+E2[6]-O[6]-O[7])/N

			B2[i][j][0]=B2[i][j][0]*(N-1)/N+(0-O[0])**2/N
			B2[i][j][1]=B2[i][j][1]*(N-1)/N+(E2[1]-O[1])**2/N
			B2[i][j][2]=B2[i][j][2]*(N-1)/N+(E2[0]-O[2])**2/N
			B2[i][j][3]=B2[i][j][3]*(N-1)/N+(E2[2]-O[3])**2/N
			B2[i][j][4]=B2[i][j][4]*(N-1)/N+(E2[4]-O[4])**2/N
			B2[i][j][5]=B2[i][j][5]*(N-1)/N+(E2[3]-O[5])**2/N
			B2[i][j][6]=B2[i][j][6]*(N-1)/N+(E2[5]-O[6])**2/N
			B2[i][j][7]=B2[i][j][7]*(N-1)/N+(E2[6]-O[7])**2/N
			B2[i][j][8]=B2[i][j][8]*(N-1)/N+(E2[6]+E2[5]-O[6]-O[7])**2/N
		
			#print y, "running e, fA, fC, r, sAC, sCA, dAC, DAC, ZAC"
			#print y, "MSE  (x10^4) P>0.00", ', '.join(map(str, (round(x*10**4, 2) for x in B2[i][j] ) ) )
			#print y, "Bias (x10^3) P>0.00", ', '.join(map(str, (round(x*10**3, 2) for x in A2[i][j] ) ) )

			#print y, "MSE  (x10^4) P>0.00", ', '.join(map(str, (round(x*10**4, 2) for x in B[i][j] ) ) )
			#print y, "Bias (x10^3) P>0.00", ', '.join(map(str, (round(x*10**3, 2) for x in A[i][j] ) ) )

		N=0
		for y in range(0, NUMBER_REP):
			os.system("python simulate_random_IBD.py -c "+str(COV[i])+" "+sigma[j]+" -I True -N "+str(NUMBER_LOCI)+" -M "+str(CM[0])+" > sim/temp.txt")
			File=open("sim/temp.txt")
			E=map(float, File.readline().split()[0:11])
			E2=map(float, File.readline().split()[0:11])
			File.close()

			N+=1
			os.system("cat test2.txt | python filter.py "+str(0.99)+"> sim/test2.txt")
			os.system("cat sim/test2.txt | mapgd relatedness > sim/temp2.txt")
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

			AI[i][j][0]=AI[i][j][0]*(N-1)/N-(0-O[0])/N
			AI[i][j][1]=AI[i][j][1]*(N-1)/N-(E[1]-O[1])/N
			AI[i][j][2]=AI[i][j][2]*(N-1)/N-(E[0]-O[2])/N
			AI[i][j][3]=AI[i][j][3]*(N-1)/N-(E[2]-O[3])/N
			AI[i][j][4]=AI[i][j][4]*(N-1)/N-(E[4]-O[4])/N
			AI[i][j][5]=AI[i][j][5]*(N-1)/N-(E[3]-O[5])/N
			AI[i][j][6]=AI[i][j][6]*(N-1)/N-(E[5]-O[6])/N
			AI[i][j][7]=AI[i][j][7]*(N-1)/N-(E[6]-O[7])/N
			AI[i][j][8]=AI[i][j][8]*(N-1)/N-(E[5]+E[6]-O[6]-O[7])/N

			BI[i][j][0]=BI[i][j][0]*(N-1)/N+(0-O[0])**2/N
			BI[i][j][1]=BI[i][j][1]*(N-1)/N+(E[1]-O[1])**2/N
			BI[i][j][2]=BI[i][j][2]*(N-1)/N+(E[0]-O[2])**2/N
			BI[i][j][3]=BI[i][j][3]*(N-1)/N+(E[2]-O[3])**2/N
			BI[i][j][4]=BI[i][j][4]*(N-1)/N+(E[4]-O[4])**2/N
			BI[i][j][5]=BI[i][j][5]*(N-1)/N+(E[3]-O[5])**2/N
			BI[i][j][6]=BI[i][j][6]*(N-1)/N+(E[5]-O[6])**2/N
			BI[i][j][7]=BI[i][j][7]*(N-1)/N+(E[6]-O[7])**2/N
			BI[i][j][8]=BI[i][j][8]*(N-1)/N+(E[6]+E[5]-O[6]-O[7])**2/N

			AI2[i][j][0]=AI2[i][j][0]*(N-1)/N-(0-O[0])/N
			AI2[i][j][1]=AI2[i][j][1]*(N-1)/N-(E2[1]-O[1])/N
			AI2[i][j][2]=AI2[i][j][2]*(N-1)/N-(E2[0]-O[2])/N
			AI2[i][j][3]=AI2[i][j][3]*(N-1)/N-(E2[2]-O[3])/N
			AI2[i][j][4]=AI2[i][j][4]*(N-1)/N-(E2[4]-O[4])/N
			AI2[i][j][5]=AI2[i][j][5]*(N-1)/N-(E2[3]-O[5])/N
			AI2[i][j][6]=AI2[i][j][6]*(N-1)/N-(E2[5]-O[6])/N
			AI2[i][j][7]=AI2[i][j][7]*(N-1)/N-(E2[6]-O[7])/N
			AI2[i][j][8]=AI2[i][j][8]*(N-1)/N-(E2[5]+E2[6]-O[6]-O[7])/N

			BI2[i][j][0]=BI2[i][j][0]*(N-1)/N+(0-O[0])**2/N
			BI2[i][j][1]=BI2[i][j][1]*(N-1)/N+(E2[1]-O[1])**2/N
			BI2[i][j][2]=BI2[i][j][2]*(N-1)/N+(E2[0]-O[2])**2/N
			BI2[i][j][3]=BI2[i][j][3]*(N-1)/N+(E2[2]-O[3])**2/N
			BI2[i][j][4]=BI2[i][j][4]*(N-1)/N+(E2[4]-O[4])**2/N
			BI2[i][j][5]=BI2[i][j][5]*(N-1)/N+(E2[3]-O[5])**2/N
			BI2[i][j][6]=BI2[i][j][6]*(N-1)/N+(E2[5]-O[6])**2/N
			BI2[i][j][7]=BI2[i][j][7]*(N-1)/N+(E2[6]-O[7])**2/N
			BI2[i][j][8]=BI2[i][j][8]*(N-1)/N+(E2[6]+E2[5]-O[6]-O[7])**2/N
		
			#print y, "running e, fA, fC, r, sAC, sCA, dAC, DAC, ZAC"
			#print y, "MSE  (x10^4) P>0.00", ', '.join(map(str, (round(x*10**4, 2) for x in BI2[i][j] ) ) )
			#print y, "Bias (x10^3) P>0.00", ', '.join(map(str, (round(x*10**3, 2) for x in AI2[i][j] ) ) )
			#print y, "MSE  (x10^4) P>0.00", ', '.join(map(str, (round(x*10**4, 2) for x in BI[i][j] ) ) )
			#print y, "Bias (x10^3) P>0.00", ', '.join(map(str, (round(x*10**3, 2) for x in AI[i][j] ) ) )
#	print y, "all-2", B

print "\\begin{table}"
print "\\caption{. \\textbf{The effect of trimmming low-frequency alleles on the Bias ($\\times 10^3$) and MSE ($\\times 10^4$) of mapgd} in estimating coancestry ($\\Theta$), cofraternity ($\\Delta_{\\ddot X \\ddot Y}$), and inbreeding ($f_X$) for outbred and inbred siblings. The expected value of all 7 IBD coefficients are listed in table \\ref{tab:ibdexamples}. Observed values of coefficients are calculated based off of actual genome sharing. Results are from "+str(NUMBER_REP)+" simulations of "+str(NUMBER_LOCI)+" SNPs with free recombination.}"
print "\\label{tab:linkage_sim}"
print "\\begin{tabular}{l*{13}{c}r}"
print "\\hline"
print  "\\textbf{Coancestry} &  & \\multicolumn{2}{c}{outbred ($\\Theta=\\frac {1}{4}$)} & & \\multicolumn{2}{c}{inbred ($\\Theta=\\frac {1}{2}$)}&  & \\multicolumn{2}{c}{outbred ($\\Theta=\\frac {1}{4}$)} & & \\multicolumn{2}{c}{inbred ($\\Theta=\\frac {1}{2}$)}\\\\"
print  "\cline {3-4}"
print "\\cline {6-7}"
print "\\cline {9-10}"
print "\\cline {12-13}"
print "Chromosome size (cM) & Cov. &  Bias &  MSE & & Bias  &  MSE & &  Bias &  MSE & & Bias  &  MSE \\\\"
print "\\hline"

for i in [0,2,1]:
	for j in range (0, 3):
		print print1[j], "&", print2[i], "&", round( (A[i][j][3])*10**3,2), "&" , round( (B[i][j][3])*10**4, 2), "&&", round( (AI[i][j][3])*10**3,2), "&" , round( (BI[i][j][3])*10**4, 2), "\\\\"
# "&&", round( (A2[i][j][3])*10**3,2), "&" , round( (B2[i][j][3])*10**4, 2), "&&", round( (AI2[i][j][3])*10**3,2), "&" , round( (BI2[i][j][3])*10**4, 2), "\\\\"
	print "\\hline"

print "\\hline"
print  "\\textbf{Fraternity} &  & \\multicolumn{2}{c}{outbred ($\\Delta=\\frac {1}{4}$)} & & \\multicolumn{2}{c}{inbred ($\\Delta=\\frac {3}{8}$)}&  & \\multicolumn{2}{c}{outbred ($\\Delta=\\frac {1}{4}$)} & & \\multicolumn{2}{c}{inbred ($\\Delta=\\frac {3}{8}$)}\\\\"
print  "\cline {3-4}"
print "\\cline {6-7}"
print "\\cline {9-10}"
print "\\cline {12-13}"
print "Chromosome size (cM) & Cov. &  Bias &  MSE & & Bias  &  MSE & &  Bias &  MSE & & Bias  &  MSE \\\\"
print "\\hline"

for i in [0,2,1]:
	for j in range (0, 3):
		print print1[j], "&", print2[i], "&", round( (A[i][j][7])*10**3,2), "&" , round( (B[i][j][7])*10**4, 2), "&&", round( (AI[i][j][7])*10**3,2), "&" , round( (BI[i][j][7])*10**4, 2), "\\\\" 
#"&&", round( (A2[i][j][7])*10**3,2), "&" , round( (B2[i][j][7])*10**4, 2), "&&", round( (AI2[i][j][7])*10**3,2), "&" , round( (BI2[i][j][7])*10**4, 2), "\\\\"
	print "\\hline"

print "\\hline"
print  "\\textbf{Inbreeding} &  & \\multicolumn{2}{c}{outbred ($f=0$)} & & \\multicolumn{2}{c}{inbred ($f=\\frac {1}{2}$)} & & \\multicolumn{2}{c}{outbred ($f=0$)} & & \\multicolumn{2}{c}{inbred ($f=\\frac {1}{2}$)}\\\\"
print  "\cline {3-4}"
print "\\cline {6-7}"
print "\\cline {9-10}"
print "\\cline {12-13}"
print "Chromosome size (cM) & Cov. &  Bias &  MSE & & Bias  &  MSE & &  Bias &  MSE & & Bias  &  MSE \\\\"
print "\\hline"

for i in [0,2,1]:
	for j in range (0, 3):
		print print1[j], "&", print2[i], "&", round( (A[i][j][1]+A[i][j][2])*10**3/2.,2), "&" , round( (B[i][j][1]+B[i][j][2])/2.*10**4, 2), "&&", round( (AI[i][j][1]+A2[i][j][2])/2.*10**3,2), "&" , round( (BI[i][j][1]+B2[i][j][2])/2.*10**4, 2), "\\\\"
# "&&", round( (A2[i][j][1]+A2[i][j][2])*10**3/2.,2), "&" , round( (B2[i][j][1]+B2[i][j][2])/2.*10**4, 2), "&&", round( (AI2[i][j][1]+AI2[i][j][2])/2.*10**3,2), "&" , round( (BI2[i][j][1]+BI2[i][j][2])/2.*10**4, 2), "\\\\"
	print "\\hline"
print "\\end{tabular}%"
print "\\end{table}"

