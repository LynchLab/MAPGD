#!/bin/bash

bro="spitze-pro.out.gz"
mapgd="../../bin/mapgd"
long=`$mapgd proview -p $bro | wc -l`     #number of lines in the long test file
a="allele"
b="pro"
export OMP_NUM_THREADS=1
time1=`(/usr/bin/time -f "%e" $mapgd allele -i $bro) 2>&1 > /dev/null`
export OMP_NUM_THREADS=8
time8=`(/usr/bin/time -f "%e" $mapgd allele -i $bro) 2>&1 > /dev/null`
rm -f $a.out
echo "Multi-threading with OMP_NUM_THREADS: (1) $time1 (8) $time8 "
d=`echo "scale=2; $long/$time8" |bc -l`
echo "Processing ~$d lines/sec"

export OMP_NUM_THREADS=1
time2=`(/usr/bin/time -f "%e" $mapgd allele -t 2 -i $bro) 2>&1 > /dev/null`
echo "2 $time2"
time3=`(/usr/bin/time -f "%e" $mapgd allele -t 3 -i $bro) 2>&1 > /dev/null`
echo "3 $time3"
time4=`(/usr/bin/time -f "%e" $mapgd allele -t 4 -i $bro) 2>&1 > /dev/null`
echo "4 $time4"
time5=`(/usr/bin/time -f "%e" $mapgd allele -t 5 -i $bro) 2>&1 > /dev/null`
echo "5 $time5"
time6=`(/usr/bin/time -f "%e" $mapgd allele -t 6 -i $bro) 2>&1 > /dev/null`
echo "6 $time6"
time7=`(/usr/bin/time -f "%e" $mapgd allele -t 7 -i $bro) 2>&1 > /dev/null`
echo "7 $time7"
time9=`(/usr/bin/time -f "%e" $mapgd allele -t 9 -i $bro) 2>&1 > /dev/null`
echo "9 $time9"
time10=`(/usr/bin/time -f "%e" $mapgd allele -t 10 -i $bro) 2>&1 > /dev/null`
echo "10 $time10"
time11=`(/usr/bin/time -f "%e" $mapgd allele -t 11 -i $bro) 2>&1 > /dev/null`
echo "11 $time11"
time12=`(/usr/bin/time -f "%e" $mapgd allele -t 12 -i $bro) 2>&1 > /dev/null`
echo "12 $time12"

rm -f $a.out
echo "Multi-threading with -t: (1) $time1 (2) $time2 (3) $time4 (4) $time8 "
echo "Multi-threading with -t: (5) $time16 (6) $time32 (7) $time64 (8) $time128 "
echo "Multi-threading with -t: (9) $time16 (10) $time32 (11) $time64 (12) $time128 "
d=`echo "scale=2; $long/$time8" |bc -l`

#export OMP_NUM_THREADS=1
#time1=`(/usr/bin/time -f "%e" mpirun -np 1 $mapgd allele -i $bro) 2>&1 > /dev/null`
#time2=`(/usr/bin/time -f "%e" mpirun -np 2 $mapgd allele -i $bro) 2>&1 > /dev/null`
#time4=`(/usr/bin/time -f "%e" mpirun -np 4 $mapgd allele -i $bro) 2>&1 > /dev/null`
#time8=`(/usr/bin/time -f "%e" mpirun -np 8 $mapgd allele -i $bro) 2>&1 > /dev/null`

#echo "mpirun: (1) $time1 (2) $time2 (4) $time4 (8) $time8 "
#d=`echo "scale=2; $long/$time8" |bc -l`
#echo "Processing ~$d lines/sec"
