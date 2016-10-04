#!/bin/bash

#running these simulations requires ~? of disk space.
mapgd=../../bin/mapgd

samtools view -H fastq/pA0.fastq.sort.bam > output/simulated-header.txt

samtools mpileup fastq/*.sort.bam -f output/ref.fasta | $mapgd proview -H output/simulated-header.txt > output/simulated-pro.out
cat output/simulated-pro.out | $mapgd allele | $mapgd filter -q 0.05 -Q 0.45 -p 20 > output/simulated-filtered-map.out
$mapgd genotype -m output/simulated-filtered-map.out -p output/simulated-pro.out > output/simulated-gcf.out
cat output/simulated-gcf.out | $mapgd relatedness > mapgd.relate


#IBD calculations
plink2
SAMtools
king

#Pooled snp calling
BRESEQ

#LD calculations


#SNP calling
angdsd
GATK
SOAPsnp
SAMtools
SNVer
