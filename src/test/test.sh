#format test

testa() 
{
if [ $? -ne 0 ]; then
	echo "[$a] produced unexpected behavior."
	exit
fi
if [ `wc $a.out -l | cut -d  ' ' -f 1` -eq $size ]; then
	echo " PASS"
	rm -f $a.out
else
	echo `wc $a.out -l `
	echo "expected $size"
	echo "[$a] $msg FAIL"
	exit
fi
}

header="spitze-header.txt"
mpileup="spitze-mpileup.txt"
pro="spitze.pro"
idx="spitze.idx"
idx_pro="spitze.txt"
pro_binary="spitze.bin"
name="name-file.txt"

formats=($mpileup)

pop=96			#number of individuals in population
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
	exit
fi

msg="    "

commands=("proview pool allele filter genotype sam2idx relatedness linkage")  

a="proview"
msg="proview"
rm -f $a.out
size=$(($pro_size+$pro_header+$idx_size+$idx_header))
echo -n "cat $mpileup | $mapgd $a -H $header > $a.out							"
timeout 5s bash -c "cat $mpileup | $mapgd $a -H $header > $a.out"
testa

a="proview"
msg="proview"
rm -f $a.out
size=$(($pro_size+$pro_header+$idx_size+$idx_header))
echo -n "$mapgd $a -n $name -H $header > $a.out								"
timeout 5s bash -c "$mapgd $a -n $name -H $header > $a.out"
testa

a="proview"
msg="proview"
rm -f $a.out
size=$(($pro_size+$pro_header))
echo -n "$mapgd $a -n $name -H $header -o $a; ($a.pro)							"
timeout 5s bash -c "$mapgd $a -n $name -H $header -o $a"
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
echo -n "$mapgd proview -n $name -H $header | $mapgd $a -M 1 > $a.out				"
timeout 5s bash -c "$mapgd proview -n $name -H $header | $mapgd $a -M 1 > $a.out"
testa

a="pool"
msg="pool"
rm -f $a.out
size=$(($pro_size+$pro_header+$idx_size+$idx_header))
echo -n "$mapgd proview -n $name -H $header | $mapgd $a > $a.out						"
timeout 5s bash -c "$mapgd proview -n $name -H $header | $mapgd $a > $a.out"
testa

a="filter"
msg="filter"
rm -f $a.out
size=$(($good_size+$pro_header+$idx_size+$idx_header))
echo -n "$mapgd proview -n $name -H $header | $mapgd allele -M 1 | $mapgd $a > $a.out	"
timeout 5s bash -c "$mapgd proview -n $name -H $header | $mapgd allele -M 1 | $mapgd $a > $a.out			"
testa

a="sam2idx"
msg="sam2idx"
rm -f $a.out
size=$(($idx_size+$idx_header))
echo -n "cat $header | $mapgd $a > $a.out										"
timeout 5s bash -c "cat $header | $mapgd $a > $a.out"
testa

a="sam2idx"
msg="sam2idx"
rm -f $a.out
size=$(($idx_size+$idx_header))
echo -n "$mapgd $a -H $header > $a.out										"
timeout 5s bash -c "$mapgd $a -H $header > $a.out"
testa

a="genotype"
msg="genotype"
rm -f $a.out
size=$(($idx_size+$idx_header+$good_size+$gcf_header))
$mapgd proview -H $header -i $mpileup -o temp
$mapgd allele -i temp.pro -o temp -M 1
$mapgd filter -i temp.map -o temp-filtered
echo -n "$mapgd $a -p temp.pro -m temp-filtered.map > $a.out								"
timeout 5s bash -c "$mapgd $a -p temp.pro -m temp-filtered.map > $a.out"
testa
rm temp*

a="genotype"
msg="genotype"
rm -f $a.out
size=$(($idx_size+$idx_header+$good_size+$gcf_header))
$mapgd proview -H $header -i $mpileup > temp_pro.out
$mapgd allele -i temp_pro.out -M 1 > temp_map.out
$mapgd filter -i temp_map.out > temp_map_filtered.out
echo -n "$mapgd $a -p temp_pro.out -m temp_map_filtered.out > $a.out							"
timeout 5s bash -c "$mapgd $a -p temp_pro.out -m temp_map_filtered.out > $a.out"
testa
rm temp*

echo -n "$mapgd $a -p pro -m map > $a.out (fifopipes)									"
mkfifo map
mkfifo pro
size=$(($idx_size+$idx_header+$good_size+$gcf_header))
$mapgd proview -H $header -i $mpileup | tee pro | $mapgd allele -M 1 | $mapgd filter > map &
$mapgd genotype -p pro -m map > $a.out
rm -f map
rm -f pro
testa
