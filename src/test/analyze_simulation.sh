#!/bin/bash

#running these simulations requires ~? of disk space.

samtools view -H fastq/pA0.fastq.sort.bam > output/simulated-header.txt
samtools mpileup fastq/*.sort.bam -f output/ref.fasta | mapgd proview -H output/simulated-header.txt > output/simulated-pro.out
cat output/simulated-pro.out | mapgd allele > output/simulated-map.out
cat output/simulated-map.out | mapgd filter -q 0.05 -Q 0.45 -p 20 > output/simulated-filtered-map.out
mapgd genotype -m output/simulated-filtered-map.out -p output/simulated-pro.out > output/simulated-gcf.out
cat output/simulatd-gcf.out | mapgd relatedness > output/simulated.rel
