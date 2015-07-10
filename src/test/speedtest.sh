o "1"
export OMP_NUM_THREADS=1
echo `time ../../bin/mapgd ei -i test.pro -o ei.out`
echo "2"
export OMP_NUM_THREADS=2
echo `time ../../bin/mapgd ei -i test.pro -o ei.out`
echo "4"
export OMP_NUM_THREADS=4
echo `time ../../bin/mapgd ei -i test.pro -o ei.out`
echo "8"
export OMP_NUM_THREADS=8
echo `time ../../bin/mapgd ei -i test.pro -o ei.out`
echo "16"
export OMP_NUM_THREADS=16
echo `time ../../bin/mapgd ei -i test.pro -o ei.out`
