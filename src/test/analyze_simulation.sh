#!/bin/bash

#running these simulations requires ~? of disk space.

#samtools view -H fastq/pA0.fastq.sort.bam > simulated-header.txt
#samtools mpileup fastq/*.sort.bam | mapgd proview -I simulated-header.txt > output/simulated.pro
#cat output/simulated.pro | mapgd allele > output/simulated.map
#cat output/simulated.map | mapgd filter -q 0.05 -Q 0.45 -p 20 > output/simulated-filtered.map
mapgd genotype -m output/simulated-filtered.map -p output/simulated.pro > output/simulated.gcf
cat output/simulatd.gcf | mapgd relatedness > output/simulated.rel
