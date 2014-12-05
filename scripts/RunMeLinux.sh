#!/bin/bash
filename=$1
echo $filename
#if -n $1; then 
#	echo "The mpileup file to be analized:"
#	read filename
#fi
sam2pro-PAF $filename > infile.txt
../bin/MapGD infile.txt outfile.txt
