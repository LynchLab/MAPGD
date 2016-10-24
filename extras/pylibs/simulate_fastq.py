import random
import numpy
import math

#@ Generates a fasta file 
#Duplication rate, 
#Denovo rate,
#Deletion rate,
#Transposition rate,
#Number of classes of transposable elements

def fasta(size):
	ret=[]
	for x in range(0, size):
		ret.append(random.choice("ACGT") )
	return ret 

def spectrum ():
	rand=N**random.random()/N
	while (rand<0.0001):
		rand=N**random.random()/N
	return rand	

class locus:
	def __init__ (self, LG, pos, Freq, ref, alt):
		self.LG=LG
		self.pos=pos
		self.Freq=Freq
		self.states=[]
		self.ref=ref
		self.alt=alt
		self.call=[self.ref, self.alt]
	def random (self):
		return (random.random()<self.Freq)

mutdex=['A', 'C', 'G', 'T']

def mut (base):
	r=random.randint(0, 2)
	if base=='A':
		r=(1+r)%4
		return mutdex[r]
	if base=='C':
		r=(2+r)%4
		return mutdex[r]
	if base=='G':
		r=(3+r)%4
		return mutdex[r]
	if base=='T':
		return mutdex[r]

class Genome:
	def __init__ (self, sizes, N, psizes, pN):
		self.fasta=[]
		self.loci=[]
		normsizes=[]
		refsize=sum(sizes)
		self.LG=len(sizes)
		for size in sizes:
			normsizes.append(float(size)/float(refsize) )
		snpsperLG=numpy.random.multinomial(N, normsizes)
		for x in range(0, len(sizes) ):
			self.fasta.append(fasta_generator(sizes[x]))
		for x in range(0, len(sizes) ):
			pos=random.sample(xrange(sizes[x]), snpsperLG[x])
			pos.sort()
			for p in pos :
				ref=self.fasta[x][p]
				alt=mut(ref)
				self.loci.append(locus(x, p, spectrum(), ref, alt ) )
		normsizes=[]
		psize=sum(psizes)
		for size in psizes:
			x=random.randint(0, self.LG-1)
			p=random.randit(0, len(self.fasta[x])-1)
			self.fasta.append(self.fasta[x][p:p+size])
		for size in psizes:
			normsizes.append(float(size)/float(psize) )
		snpsperPlog=numpy.random.multinomial(pN, normsizes)
		for x in range(0, len(psizes) ):
			pos=random.sample(xrange(psizes[x]), snpsperPlog[x])
			pos.sort()
			for p in pos :
				ref=self.fasta[x+self.LG][p]
				alt=mut(ref)
				self.loci.append(locus(x, p, spectrum(), ref, alt ) )
	def printfasta (self, filename):
		File=open(filename, "w")
		for z in range(0, len(self.fasta) ):	
			File.write(">scaffold_"+str(z)+'\n')	
			for x in range(0, len(self.fasta[z]), 60):
				File.write(''.join(self.fasta[z][x:x+60])+'\n' )
		File.close()	

comp={'A':'T', 'C':'G', 'G':'C', 'T':'A'}

def rcomp(seq):
	ret=[]
	for s in reversed(seq):
		ret.append(comp[s])
	return ret

	
class Diploid:
	def __init__(self, states):
		self.states=states

	@classmethod
	def random(cls, ref):
		states=[]
		for locus in ref.loci:
			states.append([locus, locus.random(), locus.random()])
		return cls(states)

	@classmethod
	def fusion (cls, gamete1, gamete2):
		states=[]
		for x in range(0, len(gamete1) ):
			states.append([gamete1[x][0], gamete1[x][1], gamete2[x][1]] )
		return cls(states)

	def gamete (self):
		states=[]
		thisLG=0
		thispos=0
		phase=0
		for locus in self.states:
			if locus[0].LG!=thisLG:
				phase=(random.random()>0.5)
			if random.random()<0.5*(1-math.exp(-2*(locus[0].pos-thispos)*cMperbp) ):
				phase=not(phase)
			states.append([locus[0], locus[1+phase]])
		return states

	def makefastq (self, ref, NoR, errorrate, filename):
		refM=[]
		refP=[]
		self.M=[]
		self.P=[]
		for fasta in ref.fasta:
			refM.append(fasta[:])
			refP.append(fasta[:])
		for state in self.states:
			refM[state[0].LG][state[0].pos]=state[0].call[state[1]]
			refP[state[0].LG][state[0].pos]=state[0].call[state[2]]
		for x in range (0, len(refM) ):
			self.M.append(''.join(refM[x]) )
		for x in range (0, len(refP) ):
			self.P.append(''.join(refP[x]) )
		File=open(filename, "w")
		Q=int(-10*math.log(errorrate)/math.log(10) )
		for z in xrange(0, NoR):
			x=random.randint(0,len(self.P)-1 )
			y=random.randint(0,len(self.P[x])-150 )
			if random.random()>0.5:
				seq=list(self.P[x][y:y+150])
			else:
				seq=list(self.M[x][y:y+150])
			#ADD ERRORS!
			E=numpy.random.binomial(150, errorrate)
			Epos=random.sample(xrange(0, 150), E)
			for E in Epos:
				seq[E]=mut(seq[E])
			if random.random()>0.5:
				File.write("@INSTRUMENT_NAME:1:1:0:#0\n"+''.join(seq)+"\n+\n"+chr(Q+33)*len(seq)+'\n')
			else:
				File.write("@INSTRUMENT_NAME:1:1:0:#0\n"+''.join(rcomp(seq))+"\n+\n"+chr(Q+33)*len(seq)+'\n')
		File.close()
		del refM
		del refP
		del self.M
		del self.P
		


def mate (p1, p2):
	return Diploid.fusion(p1.gamete(), p2.gamete() )

def simquartet(locus, N, e):
	X=locus[1]+locus[2]
	if X==0:
		[M, m, E1, E2]=numpy.random.multinomial(N,[(1-e), e/3., e/3., e/3.])
	elif X==1:
		[M, m, E1, E2]=numpy.random.multinomial(N,[(1-e)/2.0+e/6., (1.-e)/2.0+e/6, e/3., e/3.])
	else:
		[M, m, E1, E2]=numpy.random.multinomial(N,[e/3., (1.-e), e/3., e/3.])
	return [M, m, E1, E2] 

def simoctet(locus, N, e):
	X=locus[1]+locus[2]
	if X==0:
		[Mf, Mr, mf, mr, E1f, E1r, E2f, E2r]=numpy.random.multinomial(N,[(1.-e)/2., (1.-e)/2., e/6., e/6., e/6., e/6., e/6., e/6.])
	elif X==1:
		[Mf, Mr, mf, mr, E1f, E1r, E2f, E2r]=numpy.random.multinomial(N,[(1-e)/4.0+e/12., (1.-e)/4.0+e/12., (1-e)/4.0+e/12., (1.-e)/4.0+e/12.,  e/6., e/6., e/6., e/6.])
	else:
		[Mf, Mr, mf, mr, E1f, E1r, E2f, E2r]=numpy.random.multinomial(N,[e/6., e/6., (1.-e)/2., (1.-e)/2., e/6., e/6., e/6., e/6.])
	return [Mf, Mr, mf, mr, E1f, E1r, E2f, E2r] 

def printINfile(ref, pop):
	File=open("IN_sim.txt", 'w')
	for x in range(0, len(ref.loci) ):
		out=map(str, [ref.loci[x].LG, ref.loci[x].pos, 'A'])
		for ind in pop:
			out.append('/'.join(map(str, simquartet(ind.states[x], COV, 0.001) ) ) )
		File.write('\t'.join(out)+'\n')
	File.close()

def printPop(ref, pop):
	File=open("sim.txt", 'w')
	for locus in ref.loci:
		e=0.001
		f=0
		if locus.Freq>0.5:
			p=locus.Freq
			q=1-p
			MM=p**2
			Mm=2*p*q
			mm=q**2
			h=2*p*q
			File.write('\t'.join(map(str, [locus.LG, locus.pos, 'A', 'C', 'A', COV*len(pop), p, q, e, 0, f, MM, Mm, mm, h, 100, 100, 0] ) )+'\n' )
		else:
			p=1.-locus.Freq
			q=1-p
			MM=p**2
			Mm=2*p*q
			mm=q**2
			h=2*p*q
			File.write('\t'.join(map(str, [locus.LG, locus.pos, 'A', 'A', 'C', COV*len(pop), p, q, e, 0, f, MM, Mm, mm, h, 100, 100, 0] ) )+'\n' )
	File.close()

def printTruth(ref, pop):
	File=open("sim.txt", 'w')
	T=len(pop)*2.
	for x in range(0, len(ref.loci) ):
		M=0
		MM=0
		Mm=0
		mm=0
		for a in pop:
			M+=a.states[x][1]
			M+=a.states[x][2]
			if (a.states[x][1]) and (a.states[x][2]):
				MM+=2.
			elif (a.states[x][1]) or (a.states[x][2]):
				Mm+=2.
			else:
				mm+=2.
		if (M)/T!=0 and (M)/T!=1.:
			File.write('\t'.join(map(str, [ref.loci[x].LG, ref.loci[x].pos, ref.loci[x].ref, float(M)/T, MM/T, Mm/T, mm/T] ) )+'\n' )
	File.close()
		

def quartet2mpileup(Q):
	N=sum(Q)
	return '\t'.join(map(str, [N, 'a'*Q[0]+'A'*Q[1]+'c'*Q[2]+'C'*Q[3]+'g'*Q[4]+'G'*Q[5]+'t'*Q[6]+'T'*Q[7], 'F'*N] ) )

def printmpileup(ref, pop):
	File=open("sim.mpileup", 'w')
	for x in range(0, len(ref.loci) ):
		out=map(str, [ref.loci[x].LG, ref.loci[x].pos, 'A'])
		for ind in pop:
			out.append(quartet2mpileup(simoctet(ind.states[x], 25, 0.001) ) )
		File.write('\t'.join(out)+'\n')
	File.close()

#ref=Genome([2*10**6, 2*10**6, 2*10**6, 2*10**6, 2*10**6, 2*10**6, 2*10**6, 2*10**6, 2*10**6, 2*10**6, 2*10**6, 2*10**6, 2*10**6, 2*10**6], 200000, [], 0)
ref=Genome([1*10**5, 1*10**5, 1*10**5, 1*10**5], 50000, [], 0)

p=[]

ref.printfasta("ref.fasta")

for x in range (0, 100):

	p1=Diploid.random(ref)
	p2=Diploid.random(ref)

	A1=mate(p1, p2)
	B1=mate(p1, p2)

	C1=mate(p1, p1)
	D1=mate(p1, p1)

	p2=Diploid.random(ref)

	A1.makefastq(ref, NCOV, 0.005, "related/pA"+str(x)+".fastq")
	B1.makefastq(ref, NCOV, 0.005, "related/pB"+str(x)+".fastq")

	C1.makefastq(ref, NCOV, 0.005, "related/pC"+str(x)+".fastq")
	D1.makefastq(ref, NCOV, 0.005, "related/pD"+str(x)+".fastq")


for x in range(0, 96):
	p.append(Diploid.random(ref) )
	p[-1].makefastq(ref, NCOV, 0.005, "unrelated/p"+str(x)+".fastq")

printTruth(ref, p)

#printINfile(ref, [p1,f1,f2,f3,f4,f5,f6,f7, p2, p3, p4, p5, p6, p7]+p)
#printmpileup(ref, [p1,f1,f2,f3,f4,f5,f6,f7, p2, p3, p4, p5, p6, p7]+p)
