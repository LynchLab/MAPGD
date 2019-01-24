#!/bin/bash

bro="spitze-pro.out.gz"
mapgd="../../bin/mapgd"
long=`$mapgd proview -p $bro | wc -l`     #number of lines in the long test file
a="allele"
b="pro"
export OMP_NUM_THREADS=1
time1=`(/usr/bin/time -f "%e" $mapgd allele -i $bro) 2>&1 > /dev/null`
export OMP_NUM_THREADS=2
time2=`(/usr/bin/time -f "%e" $mapgd allele -i $bro) 2>&1 > /dev/null`
export OMP_NUM_THREADS=4
time4=`(/usr/bin/time -f "%e" $mapgd allele -i $bro) 2>&1 > /dev/null`
export OMP_NUM_THREADS=8
time8=`(/usr/bin/time -f "%e" $mapgd allele -i $bro) 2>&1 > /dev/null`
rm -f $a.out
echo "Multi-threading: (1) $time1 (2) $time2 (4) $time4 (8) $time8 "
d=`echo "scale=2; $long/$time8" |bc -l`
echo "Processing ~$d lines/sec"

export OMP_NUM_THREADS=1
time1=`(/usr/bin/time -f "%e" mpirun -n 1 $mapgd allele -i $bro) 2>&1 > /dev/null`
time2=`(/usr/bin/time -f "%e" mpirun -n 2 $mapgd allele -i $bro) 2>&1 > /dev/null`
time4=`(/usr/bin/time -f "%e" mpirun -n 4 $mapgd allele -i $bro) 2>&1 > /dev/null`
time8=`(/usr/bin/time -f "%e" mpirun -n 8 $mapgd allele -i $bro) 2>&1 > /dev/null`

echo "mpirun: (1) $time1 (2) $time2 (4) $time4 (8) $time8 "
d=`echo "scale=2; $long/$time8" |bc -l`
echo "Processing ~$d lines/sec"
