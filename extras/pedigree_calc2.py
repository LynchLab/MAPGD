#!/usr/bin/python
import sys
import argparse

parser=argparse.ArgumentParser(description="Genotypic correlation calculator")
parser.add_argument('-n', nargs=2, metavar='rec1 rec2', type=str, help="rec1 rec2")
parser.add_argument('pedigree', metavar='FileName', type=str, help="A csv formated pedigree file. FamilyID,PersonID,FatherID,MotherID,Sex.")
args=parser.parse_args()

PersonIDs={}
independent=[]

P=0.1
Q=1-P

class record:
	def __init__ (self, args):
		self.FamilyID, self.PersonID, self.FatherID, self.MotherID, self.Sex=args
		global PersonIDs
		PersonIDs[self.PersonID]=self

		self.F=[]
		self.M=[]

		if self.FatherID=='':
			self.PF=[P, Q]	#These are unconditional
			independent.append(self.PF)
		else :
			self.PF="NaN":
		if self.MotherID=='':
			self.PM=[P, Q]
			independent.append(self.PM)
		else :
			self.PM="NaN":
		self.bellow=[]
		
	def getrel(person):
		return  
		
def inner_enum(weight, tree):
	for :
		if ?==1:
			weight*0.5
	PersonIDs={}
	return weight, tree

def outer_enum(weight, ind):
	for x in range(0. len(ind) ):
		if len(ind[x])==2:
			ind[x]=0
			outer_enum(weight*Q, ind)
			ind[x]=1
			outer_enum(weight*P, ind)
			return 0
	inner_enum(weight, ind)

def setbellow(ID, IDB):
        if ID!='':
                ID=PersonIDs[ID]
                ID.bellow.append(IDB)
                setbellow(ID.FatherID, IDB)
                setbellow(ID.MotherID, IDB)

try:
	File=open(args.pedigree, "r")
except IOError:
	print "Could not open file!"
	quit()

for line in File:
	try:
		temp=line.strip('\n').split(',')
		record(temp)
	except:
		print temp
		print "could not parse line. Please use format : FamilyID,PersonID,FatherID,MotherID,Sex." 

for rec1 in PersonIDs.keys():
	setbellow(PersonIDs[rec1].FatherID, rec1)
	setbellow(PersonIDs[rec1].MotherID, rec1)

print "ID1\tID2\tF1\tF2\tTheta12\tgamma12\tgamma21\tdelta12\tDelta12"

try:
	print args.n[0]+'\t'+args.n[1]+'\t'+'\t'.join(map (str, PersonIDs[args.n[0]].getrel(args.n[1]) ) )
except:
	for rec1 in PersonIDs.keys():
		for rec2 in PersonIDs.keys():
			print rec1+'\t'+rec2+'\t'+'\t'.join(map(str, PersonIDs[rec1].getrel(rec2) ) )
