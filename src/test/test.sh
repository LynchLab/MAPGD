#format test

testa() 
{
if [ $? -ne 0 ]; then
	echo "[$a] produced unexpected behavior."
	exit
fi
if [ `wc $a.out -l | cut -d  ' ' -f 1` -eq $size ]; then
	echo "[$a] $msg PASS"
	rm -f $a.out
else
	echo  `wc $a.out -l `
	echo  "expected $size"
	echo "[$a] $msg FAIL"
	exit
fi
}


tmf="test.tmf"
lmp="test.mpileup"
smp="test.mpileup-f"
pro="test.pro"
bro="test.bro"

formats=($lmp $pro $bro $tmf)

#??
pop=96					#number of individuals in population
proheader=3				#Size of pro header
vcfheader=1				#Size of vcf header
gcfheader=$((1+$pop))			#Size of gcf header
short=`wc -l $smp | cut -d ' ' -f 1`	#number of lines in the short test file
long=`wc -l $lmp | cut -d ' ' -f 1`	#number of lines in the long test file

mapgd="../../bin/mapgd"

msg="    "

#commands 
commands=("proview ep cp ei")  

a="proview"

for format in ${formats[@]}; do
msg="$format"
rm -f $a.out
size=$(($long+$proheader))
timeout 5s bash -c "cat $format | $mapgd $a > $a.out"
testa

done

a="ei"

for format in ${formats[@]}; do
msg="$format"
rm -f $a.out
size=$(($long+$gcfheader))
timeout 10s bash -c "cat $format | $mapgd $a > $a.out"
testa

done

a="ep"

for format in ${formats[@]}; do
msg="$format"
rm -f $a.out
size=$(($long+$vcfheader))
timeout 10s bash -c "cat $format | $mapgd $a > $a.out"
testa

done

a="cp"

for format in ${formats[@]}; do
msg="$format"
rm -f $a.out
size=$(($long+$vcfheader))
timeout 10s bash -c "cat $format | $mapgd $a > $a.out"
testa

done

msg="-i  "
rm -f $a.out
size=$(($short+$proheader))
timeout 5s $mapgd proview -i $smp > $a.out
testa

msg="-o  "
rm -f $a.out
size=$(($short+$proheader))
timeout 5s $mapgd proview -i $smp -o $a.out
testa

msg="-f  "
size=$(($long+$proheader))
rm -f $a.out
timeout 5s $mapgd proview -i $lmp -o $a.out
testa

msg="-b  "
size=$(($long+$proheader))
rm -f $a.out
timeout 5s bash -c "$mapgd proview -i $lmp -b | $mapgd proview -o $a.out" 
testa

msg="-c 5"
size=$(($short+$proheader))
rm -f $a.out
timeout 5s $mapgd proview -i $smp -c 5 -o $a.out
testa

msg="-c 6"
size=$(($short+$proheader))
rm -f $a.out
timeout 5s $mapgd proview -i $smp -c 6 -o $a.out
testa

msg="-c 7"
size=$(($short+$proheader))
rm -f $a.out
timeout 5s $mapgd proview -i $smp -c 7 > $a.out
testa

msg="	"

a="ep"
size=$(($short+$vcfheader))
rm -f $a.out
timeout 5s $mapgd ep -i $smp -o $a.out
testa

a="cp"
size=$(($short+$vcfheader))
rm -f $a.out
timeout 5s $mapgd cp -i $smp -o $a.out
testa

a="ei"
msg="-p "
size=$(($short+$proheader))
rm -f $a.out
timeout 10s $mapgd ei -i $smp -p $a.out -o temp.out
testa
rm -f temp.out

msg="	"
a="ei"
size=$(($short+$gcfheader))
rm -f $a.out
timeout 10s $mapgd ei -i $smp -o $a.out
testa

msg="merge "
a="proview"
size=$(($long+$proheader))
rm -f $a.out
timeout 5s $mapgd proview -i $smp $lmp > $a.out 
testa

echo -n "SPEED TEST ... "

export OMP_NUM_THREADS=2
time1=`(/usr/bin/time -f "%e" $mapgd ei -i $bro -o $a.out) 2>&1 > /dev/null`
export OMP_NUM_THREADS=8
time2=`(/usr/bin/time -f "%e" $mapgd ei -i $bro -o $a.out) 2>&1 > /dev/null`
rm -f $a.out
echo "(2) $time1 (8) $time2 ."
d=`echo "scale=2; $long/$time2" |bc -l`
echo "Processing ~$d lines/sec"
