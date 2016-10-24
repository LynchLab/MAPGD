IMGS=`ls *.pdf | cut -d '.' -f 1`

for img in $IMGS
do
	echo $img
	convert -density 500 $img.pdf $img.jpg
	rm $img.pdf
done
