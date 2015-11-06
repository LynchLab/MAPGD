#!/bin/python
import sys
import argparse

parser=argparse.ArgumentParser(description="Command Line argument parser")
parser.add_argument('-n', nargs=2, metavar='rec1 rec2', type=str, help="rec1 rec2")
parser.add_argument('pedigree', metavar='FileName', type=str, help="A csv formated pedigree file. FamilyID,PersonID,FatherID,MotherID,Sex.")
args=parser.parse_args()

PersonIDs={}
Thetas={}

class record:
	def __init__ (self, args):
		self.FamilyID, self.PersonID, self.FatherID, self.MotherID, self.Sex=args
		global PersonIDs
		PersonIDs[self.PersonID]=self
		if self.FatherID=='' and self.MotherID=='':
			self.f=0
		self.bellow=[]
	def getrel (self, rec2):
		global PersonIDs
		try:
			f1=self.f
		except:
			self.f=gettheta(self.FatherID, self.MotherID)
			f1=self.f
		try:
			f2=PersonIDs[rec2].f
		except:
			PersonIDs[rec2].f=gettheta(PersonIDs[rec2].FatherID, PersonIDs[rec2].MotherID)
			f2=PersonIDs[rec2].f	
		theta=gettheta(self.PersonID, rec2)
		d12=getdelta(self.PersonID, rec2)
		z12=getzed(self.PersonID, rec2)
		if rec2==self.PersonID:
			s12=self.f
			s21=self.f
		else:
			s12="NaN"
			s21="NaN"
		return f1, f2, theta, s12, s21, z12, d12

def getzed(ID1, ID2):
	global PersonIDs
	if ID1=='' or ID2=='':
		return 0
	if ID1==ID2:
		return PersonIDs[ID1].f
	ID1=PersonIDs[ID1]
	ID2=PersonIDs[ID2]
	if not(ID1.PersonID in ID2.bellow) and not(ID2.PersonID in ID1.bellow):
		return (gettheta(ID2.MotherID, ID2.FatherID)*gettheta(ID1.MotherID, ID2.FatherID)
				+gettheta(ID1.FatherID, ID1.MotherID)*gettheta(ID1.FatherID, ID2.MotherID) )
	else:
		if (ID2.PersonID in ID1.bellow):
			return (gettheta(ID1.PersonID, ID2.FatherID)*gettheta(ID1.PersonID, ID2.MotherID) )*2
		else:
			return (gettheta(ID2.PersonID, ID1.FatherID)*gettheta(ID2.PersonID, ID1.MotherID) )*2

def getdelta(ID1, ID2):
	global PersonIDs
	if ID1=='' or ID2=='':
		return 0
	if ID1==ID2:
		return 1.-PersonIDs[ID1].f
	ID1=PersonIDs[ID1]
	ID2=PersonIDs[ID2]
	if not(ID1.PersonID in ID2.bellow) and not(ID2.PersonID in ID1.bellow):
		return (gettheta(ID1.FatherID, ID2.MotherID)*gettheta(ID1.MotherID, ID2.FatherID)
				+gettheta(ID1.FatherID, ID2.FatherID)*gettheta(ID1.MotherID, ID2.MotherID) )
	else:
		if (ID2.PersonID in ID1.bellow):
			return (gettheta(ID1.PersonID, ID1.FatherID)*gettheta(ID1.PersonID, ID2.MotherID) )*2
		else:
			return (gettheta(ID2.PersonID, ID1.FatherID)*gettheta(ID2.PersonID, ID1.MotherID) )*2


def gettheta(ID1, ID2):
	global PersonIDs
	if ID1=='' or ID2=='':
		return 0
	if ID1==ID2:
		try:
			return (1.+PersonIDs[ID1].f)/2.
		except:
			print ID1, " isn't set?"
			PersonIDs[ID1].f=gettheta(PersonIDs[ID1].MotherID, PersonIDs[ID1].FatherID)
			return (1.+PersonIDs[ID1].f)/2.
	ID1=PersonIDs[ID1]
	ID2=PersonIDs[ID2]
	if not(ID1.PersonID in ID2.bellow) and not(ID2.PersonID in ID1.bellow):
		return (gettheta(ID1.FatherID, ID2.FatherID)+gettheta(ID1.MotherID, ID2.FatherID)
				+gettheta(ID1.FatherID, ID2.MotherID)+gettheta(ID1.MotherID, ID2.MotherID) )/4.
	else:
		if (ID2.PersonID in ID1.bellow):
			return (gettheta(ID1.PersonID, ID2.FatherID)+gettheta(ID1.PersonID, ID2.MotherID) )/2.
		else:
			return (gettheta(ID2.PersonID, ID1.FatherID)+gettheta(ID2.PersonID, ID1.MotherID) )/2.

def setbellow(ID, IDB):
	if ID!='':
	#	print "putting "+IDB+" in "+ID
		ID=PersonIDs[ID]
		ID.bellow.append(IDB)
		setbellow(ID.FatherID, IDB)	
		setbellow(ID.MotherID, IDB)	
	
	
	#else:
	#	return 0.5

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
