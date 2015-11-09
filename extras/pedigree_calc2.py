#!/bin/python
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
			self.F=0
			independent.append(self.PF)
		if self.MotherID=='':
			self.PM=[P, Q]
			self.M=0
			independent.append(self.PM)
		self.bellow=[]
	def getP (self, conditions):
		self.getF()
		self.getM()
		self.G=[ (self.F[0]+self.M[0])/2., (self.F[1]+self.M[1])/2. ]

def setbellow(ID, IDB):
        if ID!='':
                ID=PersonIDs[ID]
                ID.bellow.append(IDB)
                setbellow(ID.FatherID, IDB)
                setbellow(ID.MotherID, IDB)

for event in independent:
	trees	


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
