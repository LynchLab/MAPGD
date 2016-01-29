#!/bin/bash

#running these simulations requires ~? of disk space.

#create the simulated reads
rm related/*
rm unrelated/*

python simulate_fastq.py

#two loops to align the reads to the reference
bwa index ref.fasta

for i in related/*.fastq; do
	echo $i
	bwa aln ref.fasta $i > $i.sai
	bwa samse ref.fasta $i.sai $i > $i.sam 
	samtools view -bS $i.sam > $i.bam
	rm $i.sam
	samtools sort $i.bam $i.sort 
	rm $i.bam
	samtools index $i.sort.bam
	gzip $i
done

for i in unrelated/*.fastq; do
	echo $i
	bwa aln ref.fasta $i > $i.sai
	bwa samse ref.fasta $i.sai $i > $i.sam 
	samtools view -bS $i.sam > $i.bam
	rm $i.sam
	samtools sort $i.bam $i.sort 
	rm $i.bam
	samtools index $i.sort.bam
	gzip $i
done

./analyze_simulation.sh
