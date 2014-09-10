#!/bin/bash

REF=$2;
if [ -z $2 ];
then
REF="$HOME/Populus/reference_seqs/populus.trichocarpa.cp.fasta";
fi

echo "reference is $REF";

for f in $1 #or whatever files contain the formatted input files for bwa_to_bam.py
do
	if [ -f $f ]
	then
	echo "looking in $f"
	# this cuts out the first ~2.2GB of the unsorted bam file: 10,000,000 (10M) reads is equivalent to >1000 reads/bp
	python ~/phylogenomics/python/bwa_to_bam.py -i $f -r $REF -p 8 -n 10000000
	python ~/phylogenomics/python/bam_to_vcf.py -i $f -r $REF -p 8
	fi
done
