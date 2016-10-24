import sys

P=float(sys.argv[1])

for line in sys.stdin:
	if line[0]!='@':
		temp=line.split('\t')
		if len(temp)>3:
			if float(temp[2])<P:
				print line,
		else:
			print line,
	else:
		print line,
#	if ?:
#		print line
#	else 
#		print line
