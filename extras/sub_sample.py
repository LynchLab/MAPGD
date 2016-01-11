#!/usr/bin/python
import sys
import argparse
import random

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

parser = argparse.ArgumentParser(description='prints out N random lines from a file.')
parser.add_argument('-N', metavar='--number', type=int, default=50000,
			help='the number of lines printed')
parser.add_argument('filename', metavar='filename', type=str, nargs=1,
                   help='the name of a proFile')

args = parser.parse_args()

File=open(args.filename[0])
length=file_len(args.filename[0])
#File=File.read().split('\n')

nums=[]

for i in range(1, length):
	nums.append(i)
random.shuffle(nums)
out_nums=nums[0:args.N]
x=0
out_nums=sorted(out_nums)
next_num=out_nums.pop(0)
for line in File:
	if (x==next_num):
		print line,
		if len(out_nums)>0:
			next_num=out_nums.pop(0)
		else:
			break
	x+=1
#for x in sorted(out_nums):
#	print File[x]
