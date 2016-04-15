import sys

File=open(sys.argv[1])

length=0

poly_counts=[]
fixed_counts=[]
polygroup=[]
fixedgroup=[]

minll_poly=5
maxll_poly=2
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
					fixed_counts.append([])
					for y in range(0, length):
						poly_counts[x].append(0)
						fixed_counts[x].append(0)
	else:
		if line[0]=='@':
			break
		line=line.strip('\n').split()
		line=line[6:]
		thissnp=line[1]
		for x in range(0, len(line) ):
			ll=line[x].split('/')
			if ll[1]!='...':
				if float(ll[2])>minll_poly:
					polygroup.append(x)
				elif float(ll[2])<maxll_poly and float(ll[3])>minll_fixed:
					fixedgroup.append(x)
		C=polygroup
		if len(C)!=0:
			for x in C:
				for y in C:
					poly_counts[x][y]+=1
				for y in C:
					poly_counts[x][y]+=1
		C=fixedgroup
		if len(C)!=0:
			for x in C:
				for y in C:
					fixed_counts[x][y]+=1
				for y in C:
					fixed_counts[x][y]+=1
		polygroup=[]
		fixedgroup=[]
print "Polymorphisms"
for x in range(0, length):
	out=[]
	for y in range(0, length):
		out.append(str(poly_counts[x][y]/2))
	print names[x], ','.join(out)
print "Fixed"
for x in range(0, length):
	out=[]
	for y in range(0, length):
		out.append(str(fixed_counts[x][y]/2))
	print names[x], ','.join(out)
