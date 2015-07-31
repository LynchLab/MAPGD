class mapfile :
	def __init__ (filename):
		this.tsv=open(filename)
		this.values={}
		this.keys={}
		this.rkey=[]
		line=this.tsv.readline()
		line=line.strip('\n').split('\t')
		line=filter(None, line)
		for x in range(0, len(line) ):
			keys[line[x]]=x
			rkey.append(line[x])
	def readline:
		line=this.tsv.readline()
		line=line.strip('\n').split('\t')
		line=filter(None, line)
		for x in range(0, len(line) ):
			values[rkey[x]]=? line[x]
	def merge (that_map):
		if (that_map.values["id1"]==this.values["id1"] and that_map.values["id2"]==this.values["id2"] ):
			this.values.update(that_map.values) 
		else:
			print "attempting to merge to different lines. Behavior is undefined."
			exit(0)
	def get (string):
		return this.values[string]
		
