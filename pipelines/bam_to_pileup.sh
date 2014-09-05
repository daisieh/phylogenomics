#!/bin/bash

# sample file has host sample path
SAMPLE=$1

# pairwise fasta for "genome"
FASTA=$2

# output dir
OUTDIR=$3

mkdir $OUTDIR

bwa index $FASTA

cd $OUTDIR
python $REPOS/phylogenomics/python/bwa_to_bam.py -i ../$SAMPLE -r ../$FASTA -p 8 -n 100000000


# sort the bam files
while read line
do
arr=($line);
samplename=${arr[1]};
samtools sort $samplename.bam $samplename.sorted
done < $SAMPLE
cd ..
