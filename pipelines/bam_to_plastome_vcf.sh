#!/bin/bash

for f in balsam.txt tricho.txt #or whatever files contain the formatted input files for bwa_to_bam.py
do
	f="samples/$f"
	echo $f

	# this cuts out the first ~2.2GB of the unsorted bam file: equivalent to >1000 reads/bp
	python ~/phylogenomics/python/bwa_to_bam.py -i $f -r ~/Populus/reference_seqs/populus.trichocarpa.cp.fasta -p 8 -n 10000000
	python ~/phylogenomics/python/bam_to_vcf.py -i $f -r ~/Populus/reference_seqs/populus.trichocarpa.cp.fasta -p 8
done
