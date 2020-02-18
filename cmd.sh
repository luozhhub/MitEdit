#!/usr/bin/bash
for i in `ls /home/zyang/Project/mitochondria/blood_tissue/sra/*.sra`
do 
j=` basename $i`
echo "fastq-dump --gzip --outdir /home/zyang/Project/mitochondria/blood_tissue/fastq/ --split-3 $i" > /home/zyang/Project/mitochondria/blood_tissue/pbs/${j%.*}.pbs
done
