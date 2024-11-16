#!/bin/bash

mkdir merged
mkdir merged/fastas

for i in *R1.fastq.gz
do
	base=$(basename $i _R1.fastq.gz)
	bbmerge.sh in=${base}_R1.fastq.gz in2=${base}_R2.fastq.gz out=./merged/${base}_merged_loose.fastq.gz adapters=adapters.fasta loose=t 
done

for j in ./merged/*.fastq.gz
do
	base2=$(basename $j .fastq.gz)
	any2fasta $j > ./merged/fastas/${base2}.fasta
done
