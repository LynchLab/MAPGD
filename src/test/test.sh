#format test
clean()
if [ -a $a.out ]; then
	rm $a.out
fi

testa() 
{
if [ $? -ne 0 ]; then
	echo "[$a] produced unexpected behavior."
	exit
fi
if [ `wc $a.out -l | cut -d  ' ' -f 1` -eq $size ]; then
	echo "[$a] $msg PASS"
	clean
else
	echo  `wc $a.out -l `
	echo "[$a] $msg FAIL"
	exit
fi
}



msg="    "

a="proview"

size=10002

clean
cat test.mpileup | ../../bin/mapgd proview > proview.out 
testa

msg="-i  "
clean
timeout 5s ../../bin/mapgd proview -i test.mpileup > proview.out
testa

msg="-o  "
clean
timeout 5s ../../bin/mapgd proview -i test.mpileup -o proview.out
testa

msg="-f  "
size=102
clean
cat test.mpileup-f | ../../bin/mapgd proview > proview.out 
testa

msg="-b  "
size=102
clean
cat test.mpileup-f | ../../bin/mapgd proview -b | ../../bin/mapgd proview > proview.out 
testa

msg="-c 5"
clean
cat test.mpileup-f | ../../bin/mapgd proview -c 5 > proview.out
testa

msg="-c 6"
clean
cat test.mpileup-f | ../../bin/mapgd proview -c 6 > proview.out 
testa

msg="-c 7"
clean
cat test.mpileup-f | ../../bin/mapgd proview -c 7 > proview.out
testa

msg="	"

a="ep"
size=100

clean
timeout 5s ../../bin/mapgd ep -i test.mpileup-f -o ep.out
testa

a="cp"
size=98	#ERROR!

clean
timeout 5s ../../bin/mapgd cp -i test.mpileup-f -o cp.out
testa

a="ei"
size=196

msg="-p "

clean
timeout 10s ../../bin/mapgd ei -i test.mpileup-f -p pro.out -o ei.out
testa

msg="	"
a="ei"
clean
timeout 10s ../../bin/mapgd ei -i test.mpileup-f -o ei.out
testa

size=10002
msg="merge "
a="proview"
clean
timeout 5s ../../bin/mapgd proview -i test.mpileup test.mpileup-f > proview.out 
testa

echo "SPEED TEST"

export OMP_NUM_THREADS=1
time1=`(time ../../bin/mapgd ei -i pro.out -o ei.out) 2>&1 >/dev/null | grep real | cut -d '	' -f 2 | cut -d 'm' -f 1 ` 
export OMP_NUM_THREADS=4
time2=`(time ../../bin/mapgd ei -i pro.out -o ei.out) 2>&1 >/dev/null | grep real | cut -d '	' -f 2 | cut -d 'm' -f 1 ` 
if [ $(($time1)) -gt $(($time2*2)) ]; then
echo "PASSED"
else
echo "FAILED"
fi

