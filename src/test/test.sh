#format test
clean()
if [ -a $a.out ]; then
	rm $a.out
fi

testa() 
{
sleep 1
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

size=10001

clean
cat test.mpileup | ../../bin/mapgd proview > proview.out 
testa

msg="-i  "
clean
../../bin/mapgd proview -i test.mpileup > proview.out
testa

msg="-o  "
clean
../../bin/mapgd proview -i test.mpileup -o proview.out
testa

msg="-f  "
size=101
clean
cat test.mpileup-f | ../../bin/mapgd proview > proview.out &
testa

msg="-b  "
size=101
clean
cat test.mpileup-f | ../../bin/mapgd proview -b | ../../bin/mapgd proview > proview.out 
testa

msg="-c 5"
clean
cat test.mpileup-f | ../../bin/mapgd proview -c 5 > proview.out &
testa

msg="-c 6"
clean
cat test.mpileup-f | ../../bin/mapgd proview -c 6 > proview.out &
testa

msg="-c 7"
clean
cat test.mpileup-f | ../../bin/mapgd proview -c 7 > proview.out &
testa

msg="	"

a="ep"
size=100

clean
../../bin/mapgd ep -i test.mpileup-f -o ep.out &
testa

a="cp"
size=98	#ERROR!

clean
../../bin/mapgd cp -i test.mpileup-f -o cp.out &
testa

a="ei"
size=100

msg="-p "

clean
../../bin/mapgd ei -i test.mpileup-f -p pro.out -o ei.out &
testa

msg="	"
a="ei"
clean
../../bin/mapgd ei -i test.mpileup-f -o ei.out &
testa

msg="merge files"
a="proview"
clean
../../bin/mapgd proview -i test.mpileup test.mpileup-f > proview.out & 
testa
