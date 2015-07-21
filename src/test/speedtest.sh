#!/bin/bash
bro="profile.brt"
lmp="mpileup-l.txt"
smp="mpileup-s.txt"

long=`wc -l $lmp | cut -d ' ' -f 1`     #number of lines in the long test file
mapgd="../../bin/mapgd"
a="ei"
b="pro"
export OMP_NUM_THREADS=2
time1=`(/usr/bin/time -f "%e" $mapgd ei -i $bro -o $a.out) 2>&1 > /dev/null`
export OMP_NUM_THREADS=8
time2=`(/usr/bin/time -f "%e" $mapgd ei -i $bro -o $a.out) 2>&1 > /dev/null`
rm -f $a.out
echo "(2) $time1 (8) $time2 ."
d=`echo "scale=2; $long/$time2" |bc -l`
echo "Processing ~$d lines/sec"

pro="profile.txt"
time2=`(/usr/bin/time -f "%e" $mapgd ei -i $pro -o $a.out) 2>&1 > /dev/null`
rm -f $a.out
echo " .pro file $time2 ."
d=`echo "scale=2; $long/$time2" |bc -l`

time2=`(/usr/bin/time -f "%e" $mapgd ei -i $lmp -o $a.out) 2>&1 > /dev/null`
rm -f $a.out
echo " mpileup file $time2 ."

time2=`(/usr/bin/time -f "%e" $mapgd -i $lmp $smp | $mapgd ei -o $a.out) 2>&1 > /dev/null`
rm -f $a.out
echo " merge mpileup file $time2 ."

time2=`(/usr/bin/time -f "%e" $mapgd -i $lmp $smp | $mapgd ei -o $a.out -p $b.out) 2>&1 > /dev/null`
rm -f $a.out
rm -f $b.out
echo " merge mpileup and print cleansed .pro file $time2 ."
