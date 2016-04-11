import sys

File=open(sys.argv[1])

length=0

poly_counts=[]
polygroup=[]
fixedgroup=[]

minll_poly=30
maxll_poly=0
minll_fixed=30

for line in File:
	if length==0:
		if line[0]=='@':
			if line[0:8]=="@SCFNAME":
			line=line.strip('\n').split()
			names=line[6:]
			length=len(names)
			for x in range(0, length):
				poly_counts.append([])
				for y in range(0, length):
					poly_counts[x].append(0)
	else:
		line=line.strip('\n').split()
		thissnp=line[1]
		for x in range(0, len(line) ):
			ll=line[x].split('/')
			if ll[1]>minll_poly:
				polygroup.append(x)
			elif ll[1]<maxll_poly and ll[2]>minll_fixed:
				fixedgroup.append(x)
		C=polygroup
		if len(C)!=0:
			for x in C:
				for y in C:
					poly_counts[x][y]+=1
				for y in C:
					poly_counts[x][y]+=1
		polygroup=[]
		fixedgroup=[]

for x in range(0, length+1):
	out=[]
	for y in range(0, length):
		out.append(str(poly_counts[x][y]))
	print names[x], ','.join(out)
