#format test

#echo "WARNING: STATISTICAL TEST HAVE NOT YET BEEN IMPLEMENTED. TAKE RESULTS WITH A GRAIN OF SALT UNTILL I GET SIMULATIONS UP AND RUNNING."


testa() 
{
if [ $? -ne 0 ]; then
	echo "[$a] produced unexpected behavior."
	exit
fi
if [ `wc $a.out -l | cut -d  ' ' -f 1` -eq $size ]; then
	echo " SUCCESS "
	rm -f $a.out
else
	echo  `wc $a.out -l `
	echo  "expected $size"
	echo "[$a] $msg FAIL"
	exit
fi
}

header="pA1.header"

mp2="mpileup-2.txt"
mp1="mpileup-1.txt"
tmf="takahiro.txt"
lmp="mpileup-l.txt"
smp="mpileup-s.txt"
pro="profile.txt"
bro="profile.brt"

formats=($lmp $pro $bro $tmf)

#??
pop=96					#number of individuals in population
proheader=3				#Size of pro header
vcfheader=1				#Size of vcf header
gcfheader=$((1+$pop))			#Size of gcf header
short=`wc -l $smp | cut -d ' ' -f 1`	#number of lines in the short test file
long=`wc -l $lmp | cut -d ' ' -f 1`	#number of lines in the long test file

mapgd="../../bin/mapgd"

if [ -z $mapgd ]; then 
	echo "cannot find mapgd. Did you type 'make' in the src directory?"
	exit
fi

msg="    "

#commands 
commands=("proview pool allele")  

a="proview"

for format in ${formats[@]}; do
msg="$format"
rm -f $a.out
size=$(($long+$proheader))
echo -n "cat $format | $mapgd $a > $a.out    			"
timeout 20s bash -c "cat $format | $mapgd $a > $a.out"
testa

done

a="allele"

for format in ${formats[@]}; do
msg="$format"
rm -f $a.out
size=$(($long+$gcfheader))
echo -n "cat $format | $mapgd $a > $a.out				"
timeout 30s bash -c "cat $format | $mapgd $a > $a.out"
testa

done

a="pool"

for format in ${formats[@]}; do
msg="$format"
rm -f $a.out
size=$(($long+$vcfheader))
echo -n "cat $format | $mapgd $a > $a.out				"
timeout 10s bash -c "cat $format | $mapgd $a > $a.out"
testa

done

#a="cp"
#for format in ${formats[@]}; do
#msg="$format"
#rm -f $a.out
#size=$(($long+$vcfheader))
#echo -n "cat $format | $mapgd $a > $a.out					"
#timeout 10s bash -c "cat $format | $mapgd $a > $a.out"
#testa
#
#done

msg="-i  "
rm -f $a.out
size=$(($short+$proheader))
echo -n "$mapgd proview -i $smp > $a.out				"
timeout 5s $mapgd proview -i $smp > $a.out
testa

msg="-o  "
rm -f $a.out
size=$(($short+$proheader))
echo -n "$mapgd proview -i $smp -o $a.out				"
timeout 5s $mapgd proview -i $smp -o $a.out
testa

msg="-f  "
size=$(($long+$proheader))
rm -f $a.out
echo -n "$mapgd proview -i $lmp -o $a.out				"
timeout 5s $mapgd proview -i $lmp -o $a.out
testa

msg="-b  "
size=$(($long+$proheader))
rm -f $a.out
echo -n "$mapgd proview -i $lmp -b | $mapgd proview -o $a.out	"
timeout 5s bash -c "$mapgd proview -i $lmp -b | $mapgd proview -o $a.out" 
testa

msg="-c 5"
size=$(($short+$proheader+1))
rm -f $a.out
echo -n "$mapgd proview -i $smp -c 5 -o $a.out			"
timeout 5s $mapgd proview -i $smp -c 5 -o $a.out
testa

msg="-c 6"
size=$(($short+$proheader))
rm -f $a.out
echo -n "$mapgd proview -i $smp -c 6 -o $a.out			"
timeout 5s $mapgd proview -i $smp -c 6 -o $a.out
testa

msg="-c 7"
size=$(($short+$proheader))
rm -f $a.out
echo -n "$mapgd proview -i $smp -c 7 -o $a.out			"
timeout 5s $mapgd proview -i $smp -c 7 > $a.out
testa

msg="-B"
size=$(($short+1))
rm -f $a.out
echo -n "$mapgd proview -i $smp -B > $a.out				"
timeout 5s $mapgd proview -i $smp -B > $a.out
testa


msg="	"

a="pool"
size=$(($short+$vcfheader))
rm -f $a.out
echo -n "$mapgd $a -i $smp -o $a.out				"
timeout 5s $mapgd $a -i $smp -o $a.out
testa

#a="cp"
#size=$(($short+$vcfheader))
#rm -f $a.out
#echo -n "$mapgd cp -i $smp -o $a.out					"
#timeout 5s $mapgd cp -i $smp -o $a.out
#testa

a="allele"
msg="-p "
size=$(($short+$proheader))
rm -f $a.out
echo -n "$mapgd $a -i $smp -p $a.out -o temp.out		"
timeout 10s $mapgd $a -i $smp -p $a.out -o temp.out
testa
rm -f temp.out

msg="	"
a="allele"
size=$(($short+$gcfheader))
rm -f $a.out
echo -n "$mapgd $a -i $smp -o $a.out				"
timeout 10s $mapgd $a -i $smp -o $a.out
testa

msg="merge "
a="proview"
size=$(($long+$proheader))
rm -f $a.out
echo -n "$mapgd $a -i $smp $lmp > $a.out		"
timeout 5s $mapgd $a -i $smp $lmp > $a.out 
testa

msg="error catching"
a="proview"
size=9
rm -f $a.out
echo -n "$mapgd $a -i $mp2 $mp1 > $a.out		"
timeout 5s $mapgd $a -i $mp2 $mp1 > $a.out 
testa

msg="merge -H"
a="proview"
size=33
rm -f $a.out
echo -n "$mapgd $a -i $mp2 $mp1 -H $header > $a.out		"
timeout 5s $mapgd $a -i $mp2 $mp1 -H $header > $a.out 
testa

msg="genotype"
a="genotype"
size=34
rm -f $a.out
echo -n "$mapgd allele -i $mp1 | $mapgd $a > $a.out		"
timeout 5s $mapgd allele $mp1 | $mapgd $a > $a.out
testa

msg="convert"
a="convert"
size=34
rm -f $a.out
echo -n "$mapgd $a -i  > $a.out		"
timeout 5s $mapgd $a -i $mp2 $mp1 > $a.out 
testa

msg="linkage"
a="linkage"
size=33
rm -f $a.out
echo -n "$mapgd $a -i $ldf  > $a.out		"
timeout 5s $mapgd $a -i $ldf > $a.out 
testa




echo -n "SPEED TEST ... "

export OMP_NUM_THREADS=2
time1=`(/usr/bin/time -f "%e" $mapgd proview -i $bro -b | $mapgd ei ) 2>&1 > /dev/null`
export OMP_NUM_THREADS=8
time2=`(/usr/bin/time -f "%e" $mapgd proview -i $bro -b | $mapgd ei ) 2>&1 > /dev/null`
rm -f $a.out
echo "(2) $time1 (8) $time2 ."
d=`echo "scale=2; $long/$time2" |bc -l`
echo "Processing ~$d lines/sec"
echo -n "LEAK TEST ... "

for comm in ${commands[@]}; do
	valgrind --leak-check=full --show-reachable=yes --track-origins=yes --log-file="mapgd-$comm.val" ../../bin/mapgd $comm -i $bro -o > /dev/null
done

echo "done"
