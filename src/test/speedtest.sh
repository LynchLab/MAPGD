#!/bin/bash

bro="profile.brt"
tmf="takahiro.txt"
lmp="mpileup-l.txt"
smp="mpileup-s.txt"
pro="profile.txt"

long=`wc -l $lmp | cut -d ' ' -f 1`     #number of lines in the long test file
mapgd="../../bin/mapgd"
a="ei"
b="pro"
export OMP_NUM_THREADS=1
time1=`(/usr/bin/time -f "%e" $mapgd ei -i $bro -o $a.out) 2>&1 > /dev/null`
export OMP_NUM_THREADS=8
time2=`(/usr/bin/time -f "%e" $mapgd ei -i $bro -o $a.out) 2>&1 > /dev/null`
rm -f $a.out
echo "(2) $time1 (8) $time2 "
d=`echo "scale=2; $long/$time2" |bc -l`
echo "Processing ~$d lines/sec"

time2=`(/usr/bin/time -f "%e" cat $pro | $mapgd proview -b | $mapgd ei ) 2>&1 > /dev/null`
echo -n "pro file $time2 "
d=`echo "scale=2; $long/$time2" |bc -l`
echo "Processing ~$d lines/sec"

time2=`(/usr/bin/time -f "%e" cat $bro | $mapgd proview -b | $mapgd ei ) 2>&1 > /dev/null`
echo -n "bro file $time2 "
d=`echo "scale=2; $long/$time2" |bc -l`
echo "Processing ~$d lines/sec"

time2=`(/usr/bin/time -f "%e" cat $tmf | $mapgd proview -b | $mapgd ei ) 2>&1 > /dev/null`
echo -n "takahiro file $time2 "
d=`echo "scale=2; $long/$time2" |bc -l`
echo "Processing ~$d lines/sec"

time2=`(/usr/bin/time -f "%e" cat $lmp | $mapgd proview -b | $mapgd ei ) 2>&1 > /dev/null`
echo -n "mpileup file $time2 "
d=`echo "scale=2; $long/$time2" |bc -l`
echo "Processing ~$d lines/sec"

