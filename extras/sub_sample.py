#!/usr/bin/python
import sys
import argparse
import random

parser = argparse.ArgumentParser(description='prints out N random lines from a file.')
parser.add_argument('-N', metavar='--number', type=int, default=50000,
			help='the number of lines printed')
parser.add_argument('filename', metavar='filename', type=str, nargs=1,
                   help='the name of a proFile')
args = parser.parse_args()

File=open(args.filename[0])
File=File.read().split('\n')
nums=[]
for i in range(0, len(File) ):
	nums.append(i)
random.shuffle(nums)
out_nums=nums[0:args.N]
for x in sorted(out_nums):
	print File[x]
