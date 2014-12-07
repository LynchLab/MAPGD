#!/bin/bash
filename=$1
if [ $# = 0 ]; then 
	echo "please enter the name of the mpileup file to be analyzed."
	exit
fi
if [ ! -f $filename ]; then
	echo "cannot open file" $1 
	exit
fi
if [ ! -f ../bin/sam2pro ]; then
	echo "sam2pro cannot be found. Perhaps you haven't run make?"
	exit
fi

./sam2pro $filename > infile.txt
../bin/mapgd -i ./infile.txt -o ./outfile.txt
