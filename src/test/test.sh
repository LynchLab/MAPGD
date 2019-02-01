#format test

timeout=timeout
wc=wc
type timeout >/dev/null 2>&1 || { echo >&2 "timeout required but not installed.  Checking for gtimeout."; 
timeout=gtimeout
wc=gwc
type gtimeout >/dev/null 2>&1 || { echo >&2 "not found, aborting."; exit 1; } }


testa() 
{
if [ $? -ne 0 ]; then
	echo "[$a] produced unexpected behavior."
	exit 1
fi
if [ `$wc $a.out -l | cut -d  ' ' -f 1` -eq $size ]; then
	echo " PASS"
    if [[ ${commands[$a]} == "UNTESTED" ]]
    then
        commands[$a]="PASS"
    fi
    if [[ "$CANNON" == true ]]
    then 
	cp $a.out ./output/${msg}.out
	rm -f $a.out
    else
	rm -f $a.out
    fi
else
    echo `$wc $a.out -l `
    echo "expected $size"
    echo "[$a] $msg FAIL"
    commands[$a]="FAIL"
    if [[ "$HALT" == true ]]
    then
        echo "If you would like to continue testing after this error call 'export HALT=false'"
        cat $a.out
        exit 1
    else
        echo "If you would like to halt at this error call 'export HALT=true"
        rm -f $a.out
    fi
fi
}

testb() 
{
if [ $? -ne 0 ]; then
	echo "[$a] produced unexpected behavior."
	exit 1
fi

if ! [ "$(diff -q $a.out ./output/$msg.out)" ]; then
	echo " PASS"
    if [[ ${commands[$a]} == "UNTESTED" ]]
    then
        commands[$a]="PASS"
    fi
    rm -f $a.out
else
    echo "$(diff $a.out ./output/$msg.out)" 
    echo "[$a] $msg FAIL"
    commands[$a]="FAIL"
    if [[ "$HALT" == true ]]
    then
        echo "If you would like to continue testing after this error call 'export HALT=false'"
        cat $a.out
        exit 1
    else
        echo "If you would like to halt at this error call 'export HALT=true'"
        rm -f $a.out
    fi
fi
}

testc()
{
testa
}

HALT="${HALT:-true}"
CANNON="${CANNON:-false}"

header="spitze-header.txt"
mpileup="spitze-mpileup.txt"
pro="spitze.pro"
idx="spitze.idx"
idx_pro="spitze.txt"
pro_binary="spitze.bin"
name="name-file.txt"
unicode="name-file-unicode.txt"

formats=($mpileup)

pop=8		#number of individuals in population
idx_header=3	#Size of pro header
pro_header=3	#Size of pro header
vcf_header=3	#Size of vcf header
gcf_header=3	#Size of gcf header
gof_header=3	#Size of gcf header
pro_size=30		#number of lines in the pro file
good_size=10	#number of lines post filtering
good_pool=4		#number of lines post filtering
idx_size=3		#number of lines in the idx file
gof_size=8		#number of lines in the idx file
smp_size=0		#number of lines in the smp file

mapgd="../../bin/mapgd"

if [ -z $mapgd ]; then 
	echo "cannot find mapgd. Did you type 'make' in the src directory?"
	exit 1
fi

msg="    "

declare -A commands

for com in `$mapgd -a`
do
    commands[$com]="UNTESTED"
done

a="proview"
msg="proview_1"
rm -f $a.out
size=$(($pro_size+$pro_header+$idx_size+$idx_header+$smp_size))
echo -n "cat $mpileup | $mapgd $a -H $header > $a.out							"
$timeout 5s bash -c "cat $mpileup | $mapgd $a -H $header > $a.out"
testc

a="proview"
msg="proview_1s"
rm -f $a.out
size=$(($pro_size+$pro_header+$idx_size+$idx_header+$smp_size-2))
echo -n "cat $mpileup | $mapgd $a -s -H $header > $a.out							"
$timeout 5s bash -c "cat $mpileup | $mapgd $a -s -H $header > $a.out"
testc

a="proview"
msg="proview_2"
rm -f $a.out
size=$(($pro_size+$pro_header+$idx_size+$idx_header+$smp_size))
echo -n "$mapgd $a -n $name -H $header > $a.out								"
$timeout 5s bash -c "$mapgd $a -n $name -H $header > $a.out"
testc

a="proview"
msg="proview_3"
rm -f $a.out
size=$(($pro_size+$pro_header))
echo -n "$mapgd $a -n $name -H $header -o $a; ($a.pro)							"
$timeout 5s bash -c "$mapgd $a -n $name -H $header -o $a"
mv $a.pro $a.out
testc

msg="proview_4"
size=$(($idx_size+$idx_header))
echo -n "$mapgd $a -n $name -H $header -o $a; ($a.idx)							"
mv $a.idx $a.out
testc

a="allele"
msg="allele_1"
rm -f $a.out
size=$(($pro_size+$pro_header+$idx_size+$idx_header+$gof_header+$gof_size))
echo -n "$mapgd proview -n $name -H $header | $mapgd $a -c 1 > $a.out				"
$timeout 5s bash -c "$mapgd proview -n $name -H $header | $mapgd $a -c 1 > $a.out"
testc

a="allele"
msg="allele_2"
rm -f $a.out
mkfifo out.pro
size=$(($pro_size+$pro_header))
echo -n "$mapgd proview -n $name -H $header | $mapgd $a -c 1 -p out > /dev/null				"
$timeout 5s bash -c "$mapgd proview -n $name -H $header | $mapgd $a -c 1 -p out > /dev/null &"
cat out.pro > $a.out
rm -rf out.pro
testc

a="pool"
msg="pool"
rm -f $a.out
size=$(($pro_size+$pro_header+$idx_size+$idx_header))
echo -n "$mapgd proview -n $name -H $header | $mapgd $a > $a.out						"
$timeout 5s bash -c "$mapgd proview -n $name -H $header | $mapgd $a > $a.out"
testc

a="filterpool"
msg="filterpool"
rm -f $a.out
size=$(($good_pool+$pro_header+$idx_size+$idx_header))
echo -n "$mapgd proview -n $name -H $header | $mapgd pool | $mapgd $a > $a.out	"
$timeout 5s bash -c "$mapgd proview -n $name -H $header | $mapgd pool | $mapgd $a > $a.out"
testc

a="filter"
msg="filter_1"
rm -f $a.out
size=$(($good_size+$pro_header+$idx_size+$idx_header))
echo -n "$mapgd proview -n $name -H $header | $mapgd allele -c 1 | $mapgd $a > $a.out	"
$timeout 5s bash -c "$mapgd proview -n $name -H $header | $mapgd allele -c 1 | $mapgd $a > $a.out			"
testc

a="filter"
msg="filter_2"
rm -f $a.out
size=$(($good_size+$pro_header+$idx_size+$idx_header-1))
echo    "$mapgd proview -n $name -H $header | \ "
echo -n "$mapgd allele -c 1 -b | $mapgd $a -X 0.05 > $a.out								"
$timeout 5s bash -c "$mapgd proview -n $name -H $header | $mapgd allele -c 1 -b | $mapgd $a -X 0.05 > $a.out			"
testc

a="sam2idx"
msg="sam2idx_1"
rm -f $a.out
size=$(($idx_size+$idx_header))
echo -n "cat $header | $mapgd $a > $a.out										"
$timeout 5s bash -c "cat $header | $mapgd $a > $a.out"
testc

a="sam2idx"
msg="sam2idx_2"
rm -f $a.out
size=$(($idx_size+$idx_header))
echo -n "$mapgd $a -H $header > $a.out										"
$timeout 5s bash -c "$mapgd $a -H $header > $a.out"
testc

a="genotype"
msg="genotype_1"
rm -f $a.out
size=$(($idx_size+$idx_header+$good_size+$gcf_header))
echo "$mapgd proview -H $header -i $mpileup -o temp	 								"
$mapgd proview -H $header -i $mpileup -o temp 
echo "$mapgd allele -i temp.pro -o temp -c 1		 								"
$mapgd allele -i temp.pro -o temp -c 1
echo "$mapgd filter -i temp.map -o temp-filtered	 								"
$mapgd filter -i temp.map -o temp-filtered 
echo -n "$mapgd $a -p temp.pro -m temp-filtered.map > $a.out								"
$timeout 5s bash -c "$mapgd $a -p temp.pro -m temp-filtered.map > $a.out"
testc
rm temp*

a="genotype"
msg="genotype_2"
rm -f $a.out
size=$(($idx_size+$idx_header+$good_size+$gcf_header))
echo "$mapgd proview -H $header -i $mpileup > temp_pro.out 								"
$mapgd proview -H $header -i $mpileup > temp_pro.out
echo "$mapgd allele -i temp_pro.out -c 1 > temp_map.out 								"
$mapgd allele -i temp_pro.out -c 1 > temp_map.out
echo "$mapgd filter -i temp_map.out > temp_map_filtered.out 								"
$mapgd filter -i temp_map.out > temp_map_filtered.out
echo -n "$mapgd $a -p temp_pro.out -m temp_map_filtered.out > $a.out							"
$timeout 5s bash -c "$mapgd $a -p temp_pro.out -m temp_map_filtered.out > $a.out"
testc
rm temp*

echo -n "$mapgd $a -p pro -m map > $a.out (fifopipes)									"
msg="genotype_3"
rm -f map
rm -f pro
mkfifo map
mkfifo pro
size=$(($idx_size+$idx_header+$good_size+$gcf_header))
$mapgd proview -H $header -i $mpileup | tee pro | $mapgd allele -c 1 | $mapgd filter > map &
$mapgd genotype -p pro -m map > $a.out
rm -f map
rm -f pro
testc

a="linkage"
msg="linkage"
size=17
echo -n "cat linkage | $mapgd $a > $a.out 											"
$mapgd proview -H $header -n $unicode > temp_pro
$mapgd allele -i temp_pro -c 1 > temp_map
$mapgd filter -i temp_map -p 5 > temp_filtered_map
$mapgd genotype -p temp_pro -m temp_filtered_map > temp_genotype
$mapgd linkage -i temp_genotype -M 4 > linkage.out
testc
rm -f temp*

rm -f map
rm -f pro
a="relatedness"
msg="relatedness_1"
echo -n "cat genotype.out | $mapgd $a > $a.out 									"
mkfifo map
mkfifo pro
size=$(($pop*($pop-1)/2+3))
$mapgd proview -H $header -i $mpileup | tee pro | $mapgd allele -c 1 | $mapgd filter > map &
$mapgd genotype -p pro -m map > genotype
cat genotype | $mapgd relatedness > $a.out 2> /dev/null
testc
rm -f genotype
rm -f map
rm -f pro

a="relatedness"
msg="relatedness_2"
echo -n "cat genotype.out | $mapgd $a -o $a.out 									"
rm -f map
rm -f pro
mkfifo map
mkfifo pro
size=$(($pop*($pop-1)/2+3))
$mapgd proview -H $header -i $mpileup | tee pro | $mapgd allele -c 1 | $mapgd filter > map &
$mapgd genotype -p pro -m map -o genotype 
$mapgd relatedness -i genotype.gcf -o $a.out 2> /dev/null
mv $a.out.rel $a.out
testc
rm -f $a.idx
rm -f genotype.gcf
rm -f genotype.idx
rm -f map
rm -f pro

a="unicode"
msg="unicode"
size=$(($pop*($pop-1)/2+3))
echo -n "cat genotype.out | $mapgd $a > $a.out 										"
$mapgd proview -H $header -n $unicode > temp_pro.out
$mapgd allele -i temp_pro.out -c 1 > temp_allele.out
$mapgd genotype -p temp_pro.out -m temp_allele.out > temp_genotype.out
$mapgd relatedness -i temp_genotype.out > $a.out 2> /dev/null
testc
rm -f temp*

a="readphen"
msg="readphen"
size=$(($pop+3))
echo -n "$mapgd $a -i plink.pheno > $a.out 											"
$mapgd $a -i plink.pheno > $a.out
testc
rm -f temp*

#a="keyinfo"
#msg="keyinfo"
#size=$((2))
#echo -n "$mapgd $a KEY > $a.out 												"
#$mapgd $a KEY > $a.out
#testa
#rm -f temp*

#a="help"
#msg="help"
#size=$((76))
#echo -n "$mapgd $a allele > $a.out 													"
#$mapgd $a allele > $a.out
#testa
#rm -f temp*

a="writevcf"
msg="writevcf"
size=$((42))
echo -n "$mapgd $a -g -m > $a.out 												"
$mapgd proview -H $header -n $name > temp.pro
$mapgd allele -i temp.pro -c 1 > temp.map
$mapgd genotype -p temp.pro -m temp.map > temp.gcf
$mapgd $a -g temp.gcf -m temp.map > $a.out
testa
rm -f temp*

#a="readvcf"
#msg="readvcf"
#size=$(($idx_size+$idx_header+$good_size+$gcf_header))
#echo -n "$mapgd $a -g -m > $a.out 												"
#$mapgd proview -H $header -n $name > temp.pro
#$mapgd allele -i temp.pro -c 1 > temp.map
#$mapgd genotype -p temp.pro -m temp.map > temp.gcf
#$mapgd writevcf -g temp.gcf -m temp.map > temp.vcf
#$mapgd $a -i temp.vcf  > $a.out
#testa
#rm -f temp*

#a="reltest"
#msg="reltest"
#size=$(($pop*($pop-1)/2+3))
#echo -n "$mapgd temp_genotype.out temp_relatedness.out relatedness_test.rel > $a.out 								"
#$mapgd proview -H $header -n $unicode > temp_pro.out
#$mapgd allele -i temp_pro.out -c 1 > temp_allele.out
#$mapgd genotype -p temp_pro.out -m temp_allele.out > temp_genotype.out
#$mapgd relatedness -i temp_genotype.out > temp_relatedness.out
#$mapgd $a temp_genotype.out temp_relatedness.out relatedness_test.rel > $a.out
#testa
#rm -f temp*

#  fastview	[ ]	Quickly displays contents of a file
#  filter	[x]	Filter sites in '.map' files
#  filterpool	[x]	Filter sites in '.pol' files
#  filtergcf	[ ]	Filter sites in '.gcf' files
#  genotype	[x]	Calculate genotype probabilities for individuals
#  linkage	[x]	Estimates linkage disequilibrium between loci
#  pool		[x]	Estimates allele frequencies using pooled data
#  proview	[x]	Prints data in the '.pro' file quartet format
#  relatedness	[x]	Estimates the pairwise relatedness of individuals
#  reltest	[ ]	Test for sig dif between relatedness estiamtes
#  sam2idx	[x]	Reformats a sam header file to a idx file
#  keyinfo	[x]	Displays information regarding keys (i.e. column names)
#  writevcf	[x]	Prints as a vcf file
#  writevcf2[ ]	Prints as a vcf file
#  readvcf	[ ]	Reads a vcf file
#  readbed	[ ]	Reads a bed file
#  readphen	[x]	Reads plink's pheno file
#  help		[x]	Prints helpful information
#  read		[ ]	Reads data from the SQL database
#  write	[ ]	Writes data to the SQL database


a="write"
msg="write"
size=6
rm -f test.db
echo "$mapgd sam2idx -H spitze-header.txt | $mapgd write -d test.db 								"
$mapgd sam2idx -H spitze-header.txt | $mapgd write -d test.db
echo -n "$mapgd read -d test.db -t REGIONS										"
$mapgd read -d test.db -t REGIONS > $a.out
rm -f test.db
testc

<<<<<<< HEAD
#for com in `$mapgd -a`
#do
#    printf '%-16s\t\t%s\n' $com ${commands[$com]}
#done
=======
for com in `$mapgd -a`
do
    printf '%-16s\t\t%s\n' $com ${commands[$com]}
done

echo -e "\nIf something is untested, it is probably because it isn't working...\n"
>>>>>>> cad0416f540e2b60dbeaad0647ad2e456ec7dbc3
