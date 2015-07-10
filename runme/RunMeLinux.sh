#!/bin/bash
files=$@
cm="ei"
if [ -z $files ]; then 
	echo "please enter the name of the file(s) to be analyzed."
	files=`read`
	exit
fi
for file in $files; do
	if [ ! -f $file ]; then
		echo "cannot open file:" $file
		exit
	fi
done;
if [ ! -f ../bin/mapgd ]; then
	echo "mapgd cannot be found. Perhaps you haven't run make?"
	exit
fi

if [ $((`../bin/mapgd $cm -i $files > outfile.txt`)) -eq 0 ]; then
	echo "mapgd complete execution, results are now in the outfile.txt."
	notify-send "mapgd complete execution, results are now in the outfile.txt."
fi
