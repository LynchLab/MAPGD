import random
import numpy
import os
import math

Pfreq=0.2
SIZE=1000
bp=1
GEN=8

File=open("seq.gcf", 'w')
#File10=open("seq15.gcf", 'w')
#File15=open("seq15.gcf", 'w')
#File20=open("seq20.gcf", 'w')
#File25=open("seq25.gcf", 'w')
#File30=open("seq30.gcf", 'w')
#File35=open("seq35.gcf", 'w')
#File40=open("seq40.gcf", 'w')
#File45=open("seq45.gcf", 'w')

def seq(File, P, e, F, N):
	global File2
	global bp
	lnc=math.log(1.0-e)
	lnch=math.log((1.0-e)/2.0+e/6)
	lne=math.log(e)

        File.write( '\t'.join(map(str, [1, bp, "C", "A", P, e] ) ) )
	bp+=1
	for f in F:
	        if f.gen==2:
	                [M, m, E]=numpy.random.multinomial(N,[(1-e), e/3, 2*e/3])
	        elif f.gen==1:
	                [M, m, E]=numpy.random.multinomial(N,[(1-e)/2.0+e/6, (1-e)/2.0+e/6, 2*e/3])
	        else:
	                [M, m, E]=numpy.random.multinomial(N,[e/3.0, (1-e), 2*e/3])

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
	def __init__ (self, mother=None, father=None):
		self.mother=mother
		self.father=father
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
	F[0].append(individual() )

for y in range(1, GEN):
	for x in range(0, SIZE):
		z=random.randint(0, SIZE-1)
		w=random.randint(0, SIZE-1)
		F[y].append(individual(F[y-1][z], F[y-1][w]) )
		if y==1:
			if ( not(F[0][z].labeled) ):
				Q=float(get_freq(F[0]) )/float(2*SIZE)
				F[0][z].labeled=True
			if ( not(F[0][w].labeled) ):
				Q=float(get_freq(F[0]) )/float(2*SIZE)
				F[0][w].labeled=True
		digraph_head.append("F"+str(y-1)+"_"+str(z)+" -> F"+str(y)+"_"+str(x)+";\n")
		digraph_head.append("F"+str(y-1)+"_"+str(w)+" -> F"+str(y)+"_"+str(x)+";\n")

C=individual(F[-1][0], F[-1][1])
D=individual(F[-1][0], F[-1][1])
E=individual(F[-1][2], F[-1][2])

count = [[ [0,0,0] for x in range(3)] for y in range(SIZE*2+1)]

a=0

LIM=5000000

while a<LIM:
	Pfreq=random.random()
	for y in range(0, GEN):
		for x in range(0, SIZE):
			F[y][x].rand(Pfreq)
	if get_freq(F[-1])==0 or get_freq(F[-1])==SIZE*2:
		continue
	a+=1

	Q=get_freq(F[-1]) 
	for y in range(0, GEN):
		for x in range(0, SIZE):
			F[y][x].Den[Q]+=1.
			F[y][x].Num[Q]=F[y][x].Num[Q]*(F[y][x].Den[Q]-1)/F[y][x].Den[Q]+f(F[y][x].gen, float(Q)/float(SIZE*2) )*1./F[y][x].Den[Q]

	C.rand(Pfreq)
	D.rand(Pfreq)
	E.rand(Pfreq)	

	C.Den[Q]+=1.
	C.Num[Q]=C.Num[Q]*(C.Den[Q]-1.)/C.Den[Q]+f(C.gen, float(Q)/float(SIZE*2) )*1./C.Den[Q]
	D.Den[Q]+=1.
	D.Num[Q]=D.Num[Q]*(D.Den[Q]-1.)/D.Den[Q]+f(D.gen, float(Q)/float(SIZE*2) )*1./D.Den[Q]
	E.Den[Q]+=1.
	E.Num[Q]=E.Num[Q]*(E.Den[Q]-1.)/E.Den[Q]+f(E.gen, float(Q)/float(SIZE*2) )*1./E.Den[Q]

	count[get_freq(F[-1])][C.gen][D.gen]+=1

	Qf=min(float(Q)/float(SIZE*2), float(SIZE*2-Q)/float(SIZE*2) )

#	if (Qf = 0.05 ):
	seq(File, float(Q)/float(SIZE*2), 0.01, [C, D, E]+F[-1], 50)
#	elif (Qf <=0.10 ):
#		seq(File10, float(Q)/float(SIZE*2), 0.01, [C, D, E]+F[-1], 50)
#	elif (Qf <=0.15 ):
#		seq(File15, float(Q)/float(SIZE*2), 0.01, [C, D, E]+F[-1], 50)
#	elif (Qf <=0.20 ):
#		seq(File20, float(Q)/float(SIZE*2), 0.01, [C, D, E]+F[-1], 50)
#	elif (Qf <=0.25 ):
#		seq(File25, float(Q)/float(SIZE*2), 0.01, [C, D, E]+F[-1], 50)
#	elif (Qf <=0.30 ):
#		seq(File30, float(Q)/float(SIZE*2), 0.01, [C, D, E]+F[-1], 50)
#	elif (Qf <=0.35 ):
#		seq(File35, float(Q)/float(SIZE*2), 0.01, [C, D, E]+F[-1], 50)
#	elif (Qf <=0.40 ):
#		seq(File40, float(Q)/float(SIZE*2), 0.01, [C, D, E]+F[-1], 50)
#	elif (Qf <=0.45 ):
#		seq(File45, float(Q)/float(SIZE*2), 0.01, [C, D, E]+F[-1], 50)

	if a==LIM-1:
		for Q in range(1, SIZE):
			digraph=open("temp"+format(Q, '03')+".gv", 'w')
			digraph.write(''.join(digraph_head) )
			for y in range(0, GEN):
				for x in range(0, SIZE):
					digraph.write("\"F"+str(y)+"_"+str(x)+"\" [label=\"\" fillcolor="+color(F[y][x].Num[Q])+", style=filled];" )
			digraph.write("label=\"q="+"{0:.2f}".format(float(Q)/float(SIZE*2))+"\";\n")
			digraph.write("\"a\" [fillcolor="+color(C.Num[Q])+", style=filled];" )
			digraph.write("\"b\" [fillcolor="+color(D.Num[Q])+", style=filled];" )
			digraph.write("\"c\" [fillcolor="+color(E.Num[Q])+", style=filled];" )
			digraph.write("F"+str(GEN-1)+"_0 -> a;\n")
			digraph.write("F"+str(GEN-1)+"_1 -> a;\n")
			digraph.write("F"+str(GEN-1)+"_0 -> b;\n")
			digraph.write("F"+str(GEN-1)+"_1 -> b;\n")
			digraph.write("F"+str(GEN-1)+"_2 -> c;\n")
			digraph.write("F"+str(GEN-1)+"_2 -> c;\n")
			digraph.write("}\n")
			digraph.close()
			os.system("dot -T png temp"+format(Q, '03')+".gv > file"+format(Q, '03')+".png ")
print C.Num
print D.Num
print E.Num
for x in range(1, SIZE*2):
	for y in range (0, 3):
		for z in range(0, 3):
			print count[x][y][z],
	print
