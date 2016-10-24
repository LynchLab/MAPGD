#!/bin/bash

bro="spitze-pro.out"

#long=`wc -l $lmp | cut -d ' ' -f 1`     #number of lines in the long test file
mapgd="../../bin/mapgd"
a="allele"
b="pro"
export OMP_NUM_THREADS=1
time1=`(/usr/bin/time -f "%e" $mapgd allele -i $bro) 2>&1 > /dev/null`
export OMP_NUM_THREADS=8
time2=`(/usr/bin/time -f "%e" $mapgd allele -i $bro) 2>&1 > /dev/null`
rm -f $a.out
echo "(1) $time1 (8) $time2 "
d=`echo "scale=2; $long/$time2" |bc -l`
echo "Processing ~$d lines/sec"
