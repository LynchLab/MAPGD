#!/bin/bash

#running these simulations requires ~? of disk space.

#create the simulated reads

#two loops to align the reads to the reference

#create one vcf file for each pair of related individuals

for j in $(seq 1 100); do
	mapgd filter -p 0.05 -P 0.45 -l 20 > data/related-o-$j.gcf
	~/src/MAPGD/bin/convert data/related-o-$j.gcf data/related-o-$j.bcf 
	python ~/src/MAPGD/extras/relatedness.py --bcf data/related-o-$j.bcf -a 0 -b 1 -A 0 -B 1 | tail -1 | cut -d "," -f 11-24 >> related-o-mapgd.txt

	python ~/src/MAPGD/extras/mapgdutils.py -p 0.05 -P 0.45 -l 20 --pro data/related-i-$j.pro --map data/related-i-$j.map --mode G> data/related-i-$j.gcf
	~/src/MAPGD/bin/convert data/related-i-$j.gcf data/related-i-$j.bcf 
	python ~/src/MAPGD/extras/relatedness.py --bcf data/related-i-$j.bcf -a 0 -b 1 -A 0 -B 1 | tail -1 | cut -d "," -f 11-24 >> related-i-mapgd.txt
done
