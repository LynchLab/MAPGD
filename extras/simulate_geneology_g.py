import random
import numpy
import os
import math

Pfreq=0.85
SIZE=5
bp=1
GEN=5

gdir="graphs/"

File05=open("sim/seq05.gcf", 'w')
File10=open("sim/seq10.gcf", 'w')
File15=open("sim/seq15.gcf", 'w')
File20=open("sim/seq20.gcf", 'w')
File25=open("sim/seq25.gcf", 'w')
File30=open("sim/seq30.gcf", 'w')
File35=open("sim/seq35.gcf", 'w')
File40=open("sim/seq40.gcf", 'w')
File45=open("sim/eq45.gcf", 'w')

def seq(File, P, e, F, N):
	global File2
	global bp
	lnc=math.log(1.0-e)
	lnch=math.log((1.0-e)/2.0+e/6)
	lne=math.log(e/3.)

        File.write( '\t'.join(map(str, [1, bp, "C", "A", P, e] ) ) )
	bp+=1
	for f in F:
	        if f.gen==2:
	                [M, m, E]=numpy.random.multinomial(N,[(1.-e), e/3., 2.*e/3.])
	        elif f.gen==1:
	                [M, m, E]=numpy.random.multinomial(N,[(1.-e)/2.0+e/6., (1.-e)/2.0+e/6, 2.*e/3.])
	        else:
	                [M, m, E]=numpy.random.multinomial(N,[e/3.0, (1.-e), 2.*e/3.])

		MM=lnc*M+lne*(m+E)#-math.log(P**2)
		Mm=lnch*(M+m)+lne*E#-math.log(2*P*Q)
		mm=lnc*m+lne*(M+E)#-math.log(Q**2)

		norm=math.log(math.exp(MM)+math.exp(Mm)+math.exp(mm) )

		File.write('\t'+'\t'.join(map(str, [-MM+norm, -Mm+norm, -mm+norm, N] ) ) )
	File.write('\n')

def f(h, q):
	return ( h-2*q)**2 /(2.*q*(1.-q) )-1. 

def color(num):
	if num>0:
		if (abs(num)>1):
			num=1
		HEX=int(255*abs(1-abs(num) ) )
		return "\"#"+format(HEX, "02x")+format(HEX, "02x")+"ff\""
	else:
		if (abs(num)>1):
			num=1
		HEX=int(255*abs(1-abs(num) ) )
		return "\"#ff"+format(HEX, "02x")+format(HEX, "02x")+"\""

def get_freq (array):
	Sum=0
	for ind in array:
		Sum+=ind.gen
	return Sum

class individual:
	def __init__ (self, name, mother=None, father=None):
		self.mother=mother
		self.father=father
		self.name=name
		if mother is None:
			if random.random()>Pfreq:
				self.maternal_hap=1
			else:
				self.maternal_hap=0
		else:
			if mother.gen==1:
				self.maternal_hap=random.randint(0,1)
			elif mother.gen==2:
				self.maternal_hap=1
			elif mother.gen==0:
				self.maternal_hap=0
		if father is None:
			if random.random()>Pfreq:
				self.paternal_hap=1
			else:
				self.paternal_hap=0
		else:
			if father.gen==1:
				self.paternal_hap=random.randint(0,1)
			elif father.gen==2:
				self.paternal_hap=1
			elif father.gen==0:
				self.paternal_hap=0
				
		self.gen=self.paternal_hap+self.maternal_hap
		self.labeled=False
		self.Num=[ 0 for y in range(SIZE*2) ]
		self.Den=[ 0 for y in range(SIZE*2) ]

	def rand (self, P):
		mother=self.mother
		father=self.father
		if mother is None:
			if random.random()>P:
				self.maternal_hap=1
			else:
				self.maternal_hap=0
		else: 
			if mother.gen==1:
				self.maternal_hap=random.randint(0,1)
			elif mother.gen==2:
				self.maternal_hap=1
			elif mother.gen==0:
				self.maternal_hap=0
		if father is None:
			if random.random()>P:
				self.paternal_hap=1
			else:
				self.paternal_hap=0
		else:
			if father.gen==1:
				self.paternal_hap=random.randint(0,1)
			elif father.gen==2:
				self.paternal_hap=1
			elif father.gen==0:
				self.paternal_hap=0
		self.gen=self.paternal_hap+self.maternal_hap
		return self
				
		

F = [ [] for y in range(GEN)]

digraph_head=[]
digraph_head.append("digraph G {\n")
digraph_head.append("node [shape=circle, width=0.5, ];\n")


for x in range(0, SIZE):
	F[0].append(individual("P") )

for y in range(1, GEN):
	for x in range(0, SIZE):
		w=random.randint(0, SIZE-1)
		z=random.randint(0, SIZE-1)
		F[y].append(individual("F", F[y-1][z], F[y-1][w]))
		if y==1:
			if ( not(F[0][z].labeled) ):
				Q=float(get_freq(F[0]) )/float(2*SIZE)
				F[0][z].labeled=True
			if ( not(F[0][w].labeled) ):
				Q=float(get_freq(F[0]) )/float(2*SIZE)
				F[0][w].labeled=True
		digraph_head.append("F"+str(y-1)+"_"+str(z)+" -> F"+str(y)+"_"+str(x)+";\n")
		digraph_head.append("F"+str(y-1)+"_"+str(w)+" -> F"+str(y)+"_"+str(x)+";\n")
		F[y][x].labeled=True

E=individual("E", F[-1][0], F[-1][1])
D=individual("D", F[-1][1], F[-1][2])
C=individual("C", E, D)
A=individual("C", E, F[-1][3])

count = [[ [0,0,0] for x in range(3)] for y in range(SIZE*2+1)]

a=0

LIM=10

while a<LIM:
	if a%100==0:
		print a
	Pfreq=random.random()
	for y in range(1, GEN):
		for x in range(0, SIZE):
			F[y][x].rand(Pfreq)
	if get_freq(F[-1])==0 or get_freq(F[-1])==SIZE*2:
		continue
	a+=1

	digraph=open(gdir+"temp_gene"+str(a)+".gv", 'w')
	digraph.write(''.join(digraph_head) )
	for y in range(0, GEN):
		for x in range(0, SIZE):
			if F[y][x].labeled:
				digraph.write("\"F"+str(y)+"_"+str(x)+"\" [label=\"\" fillcolor="+color(-F[y][x].gen/2.)+", style=filled];\n" )
	digraph.write("}\n")
	digraph.close()
	os.system("dot -T png "+gdir+"temp_gene"+str(a)+".gv > "+gdir+"file_gene"+str(a)+".png ")
print C.Num
print D.Num
print E.Num
for x in range(1, SIZE*2):
	for y in range (0, 3):
		for z in range(0, 3):
			print count[x][y][z],
	print
