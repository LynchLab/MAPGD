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
if [ ! -f ./sam2pro ]; then
	echo "please compile sam2pro from http://guanine.evolbio.mpg.de/mlRho/ and place the binary in this file"
        echo "sorry for the inconvenience, this step will be uneccisary soon."
	exit
fi

./sam2pro $filename > infile.txt
../bin/mapgd -i ./infile.txt -o ./outfile.txt
