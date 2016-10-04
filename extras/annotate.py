import sys

if (len(sys.argv)!=4):
	print "map file, GFF file, FNA file"
	exit()

File=open(sys.argv[1])
GFF=open(sys.argv[2])
FNA=open(sys.argv[3])

isopen=False

CUTOFF=25.0

codon_list={
"TTT":"F", "TCT":"S", "TAT":"Y", "TGT":"C",
"TTC":"F", "TCC":"S", "TAC":"Y", "TGC":"C",
"TTA":"L", "TCA":"S", "TAA":"*", "TGA":"*",
"TTG":"L", "TCG":"S", "TAG":"*", "TGG":"W",

"CTT":"L", "CCT":"P", "CAT":"H", "CGT":"R",
"CTC":"L", "CCC":"P", "CAC":"H", "CGC":"R",
"CTA":"L", "CCA":"P", "CAA":"Q", "CGA":"R",
"CTG":"L", "CCG":"P", "CAG":"Q", "CGG":"R",

"ATT":"I", "ACT":"T", "AAT":"N", "AGT":"S",
"ATC":"I", "ACC":"T", "AAC":"N", "AGC":"S",
"ATA":"I", "ACA":"T", "AAA":"K", "AGA":"R",
"ATG":"M", "ACG":"T", "AAG":"K", "AGG":"R",

"GTT":"V", "GCT":"A", "GAT":"D", "GGT":"G",
"GTC":"V", "GCC":"A", "GAC":"D", "GGC":"G",
"GTA":"V", "GCA":"A", "GAA":"E", "GGA":"G",
"GTG":"V", "GCG":"A", "GAG":"E", "GGG":"G"
}

class gene:
	def __init__ (self, name, start, stop, strand):
		self.name=name
		self.start=start
		self.stop=stop
		self.strand=strand

def is_gene(line):
	try:
		if line[2]=="CDS":
			return True
		else:
			return False
	except:
		return False

start=3
stop=4
name=8
strand=6

gene_sizes={}
genes=[]

for line in GFF:
	line=line.split('\t')
	if is_gene(line):
		for x in line[name].split(';'):
			if x.split('=')[0]=="gene":
				genes.append(gene(x.split('=')[1], int(line[start]), int(line[stop]), line[strand]) )
				gene_sizes[x.split('=')[1]]=(int(line[start])-int(line[stop]))


def rcomp (seq):
	s=[]
	rc={'A':'T', 'T':'A', 'C':'G', 'G':'C'}
	rseq=[]
	for s in seq:
		rseq.insert(0, rc[s])
	return ''.join(rseq)
		
def getmut(start, bp, seq, mut):
	A=int(bp)-(int(bp)-int(start))%3
	B=int(bp)
	return codon_list[seq[A:A+3]]+str((bp-start)/3+1)+codon_list[seq[A:B]+mut+seq[B+1:A+3]]

def getrev(start, bp, seq, mut):
	bp=int(bp)+1
	start=int(start)+1
		
	A=int(bp)+abs(int(start)-int(bp))%3
	B=int(bp)
	return codon_list[rcomp(seq[A-3:A])]+str((start-bp)/3+1)+codon_list[rcomp(seq[A-3:B-1]+mut+seq[B:A])]

def issyn(start, bp, seq, mut):
	A=int(bp)-(int(bp)-int(start))%3
	B=int(bp)
	return codon_list[seq[A:A+3]]==codon_list[seq[A:B]+mut+seq[B+1:A+3]]

def issynrev(start, bp, seq, mut):
	bp=int(bp)+1
	start=int(start)+1
		
	A=int(bp)+abs(int(start)-int(bp))%3
	B=int(bp)
	return codon_list[rcomp(seq[A-3:A])]==codon_list[rcomp(seq[A-3:B-1]+mut+seq[B:A])]

seq=[]
for line in FNA:
        if line[0]!='>':
                line=line.strip('\n')
                seq.append(line)
seq=' '+''.join(seq)

this_gene=genes.pop(0)

for line in File:
	if isopen:
		if line[0]=='@':
			print line,
			isopen=False
		fields=line.split('\t')
		ChiPoly=[]
		ChiFixed=[]
		Freq=[]
		for X in fields[6:]:
			num=X.split('/')
			try:
				ChiPoly.append(float(num[2]))
				ChiFixed.append(float(num[3]))
				P=float(num[0])
				cov=int(num[1])
				if ChiPoly[-1]<CUTOFF:
					P=round(P)
				if cov>0:
					Freq.append(P)
				else:
					Freq.append('N')
			except:
				ChiPoly.append(0)
				ChiFixed.append(0)
				Freq.append('N')
		ref=fields[2]
		alt=fields[3]
		pos=int(fields[1])
		if len(alt)==1:
			#while (pos>min(this_gene.stop, this_gene.start) ):
			while (pos>max(this_gene.start, this_gene.stop ) ):
				this_gene=genes.pop(0)

			if this_gene.strand=='+':
				if (pos>this_gene.start):
					try:
						kind=getmut(this_gene.start, pos, seq, alt)
					except:
						continue
					hit=(this_gene.name+":"+kind)
					if kind[0]==kind[-1]:
						kind='S'
					else:
						kind='N'
				else:
					hit='NC'
					kind='NC'
			elif this_gene.strand=='-':
				if (pos>this_gene.start):
					try:
						kind=getrev(this_gene.stop, pos, seq, alt)
					except:
						continue	
					hit=(this_gene.name+":"+kind)
					if kind[0]==kind[-1]:
						kind='S'
					else:
						kind='N'
				else:
					hit='NC'
					kind='NC'
		if max(ChiPoly)>CUTOFF or max(ChiFixed)>CUTOFF:
			print '\t'.join(map(str, fields[:6]))+'\t'+hit+'\t'+kind+'\t'+'\t'.join(map(str, Freq))
	else:
		print line,
		fields=line.split('\t')
		if fields[0]=="@SCFNAME       ":
			isopen=True
