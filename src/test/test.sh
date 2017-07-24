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
	rm -f $a.out
else
	echo `$wc $a.out -l `
	echo "expected $size"
	echo "[$a] $msg FAIL"
	exit 1
fi
}

header="spitze-header.txt"
mpileup="spitze-mpileup.txt"
pro="spitze.pro"
idx="spitze.idx"
idx_pro="spitze.txt"
pro_binary="spitze.bin"
name="name-file.txt"
unicode="name-file-unicode.txt"

formats=($mpileup)

pop=8			#number of individuals in population
idx_header=3		#Size of pro header
pro_header=3		#Size of pro header
vcf_header=3		#Size of vcf header
gcf_header=3		#Size of gcf header
gof_header=3		#Size of gcf header
pro_size=30		#number of lines in the pro file
good_size=9		#number of lines post filtering
idx_size=3		#number of lines in the idx file
gof_size=8		#number of lines in the idx file

mapgd="../../bin/mapgd"

if [ -z $mapgd ]; then 
	echo "cannot find mapgd. Did you type 'make' in the src directory?"
	exit 1
fi

msg="    "

commands=("proview pool allele filter genotype sam2idx relatedness linkage")  

a="proview"
msg="proview"
rm -f $a.out
size=$(($pro_size+$pro_header+$idx_size+$idx_header))
echo -n "cat $mpileup | $mapgd $a -H $header > $a.out							"
$timeout 5s bash -c "cat $mpileup | $mapgd $a -H $header > $a.out"
testa

a="proview"
msg="proview"
rm -f $a.out
size=$(($pro_size+$pro_header+$idx_size+$idx_header))
echo -n "$mapgd $a -n $name -H $header > $a.out								"
$timeout 5s bash -c "$mapgd $a -n $name -H $header > $a.out"
testa

a="proview"
msg="proview"
rm -f $a.out
size=$(($pro_size+$pro_header))
echo -n "$mapgd $a -n $name -H $header -o $a; ($a.pro)							"
$timeout 5s bash -c "$mapgd $a -n $name -H $header -o $a"
mv $a.pro $a.out
testa

size=$(($idx_size+$idx_header))
echo -n "$mapgd $a -n $name -H $header -o $a; ($a.idx)							"
mv $a.idx $a.out
testa

a="allele"
msg="allele"
rm -f $a.out
size=$(($pro_size+$pro_header+$idx_size+$idx_header+$gof_header+$gof_size))
echo -n "$mapgd proview -n $name -H $header | $mapgd $a -c 1 > $a.out				"
$timeout 5s bash -c "$mapgd proview -n $name -H $header | $mapgd $a -c 1 > $a.out"
testa

a="pool"
msg="pool"
rm -f $a.out
size=$(($pro_size+$pro_header+$idx_size+$idx_header))
echo -n "$mapgd proview -n $name -H $header | $mapgd $a > $a.out						"
$timeout 5s bash -c "$mapgd proview -n $name -H $header | $mapgd $a > $a.out"
testa

a="filter"
msg="filter"
rm -f $a.out
size=$(($good_size+$pro_header+$idx_size+$idx_header))
echo -n "$mapgd proview -n $name -H $header | $mapgd allele -c 1 | $mapgd $a > $a.out	"
$timeout 5s bash -c "$mapgd proview -n $name -H $header | $mapgd allele -c 1 | $mapgd $a > $a.out			"
testa

a="sam2idx"
msg="sam2idx"
rm -f $a.out
size=$(($idx_size+$idx_header))
echo -n "cat $header | $mapgd $a > $a.out										"
$timeout 5s bash -c "cat $header | $mapgd $a > $a.out"
testa

a="sam2idx"
msg="sam2idx"
rm -f $a.out
size=$(($idx_size+$idx_header))
echo -n "$mapgd $a -H $header > $a.out										"
$timeout 5s bash -c "$mapgd $a -H $header > $a.out"
testa

a="genotype"
msg="genotype"
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
testa
rm temp*

a="genotype"
msg="genotype"
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
testa
rm temp*

echo -n "$mapgd $a -p pro -m map > $a.out (fifopipes)									"
rm -f map
rm -f pro
mkfifo map
mkfifo pro
size=$(($idx_size+$idx_header+$good_size+$gcf_header))
$mapgd proview -H $header -i $mpileup | tee pro | $mapgd allele -c 1 | $mapgd filter > map &
$mapgd genotype -p pro -m map > $a.out
rm -f map
rm -f pro
testa

#a="linkage"
#msg="linkage"
#size=15
#echo -n "cat linkage | $mapgd $a > $a.out 											"
#$mapgd proview -H $header -n $unicode > temp_pro
#$mapgd allele -i temp_pro -c 1 > temp_map
#$mapgd filter -i temp_map -p 5 > temp_filtered_map
#$mapgd genotype -p temp_pro -m temp_filtered_map > temp_genotype
#$mapgd linkage -i temp_genotype > linkage.out
#testa
#rm -f temp*

rm -f map
rm -f pro
a="relatedness"
msg="relatedness"
echo -n "cat genotype.out | $mapgd $a > $a.out 									"
mkfifo map
mkfifo pro
size=$(($pop*($pop-1)/2+3))
$mapgd proview -H $header -i $mpileup | tee pro | $mapgd allele -c 1 | $mapgd filter > map &
$mapgd genotype -p pro -m map > genotype
cat genotype | $mapgd relatedness > $a.out
testa
rm -f genotype
rm -f map
rm -f pro

a="relatedness"
msg="relatedness"
echo -n "cat genotype.out | $mapgd $a -o $a.out 									"
rm -f map
rm -f pro
mkfifo map
mkfifo pro
size=$(($pop*($pop-1)/2+3))
$mapgd proview -H $header -i $mpileup | tee pro | $mapgd allele -c 1 | $mapgd filter > map &
$mapgd genotype -p pro -m map -o genotype 
$mapgd relatedness -i genotype.gcf -o $a.out
mv $a.out.rel $a.out
testa
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
$mapgd relatedness -i temp_genotype.out > $a.out
testa
rm -f temp*

exit 0

a="write"
msg="write/read"
rm -f test.db
echo "cat spitze.idx | $mapgd write -d test.db 									"
cat spitze.idx | $mapgd write -d test.db
echo -n "$mapgd read -d test.db -t INDEX										"
$mapgd read -d test.db -t INDEX
#rm -f test.db
testa

