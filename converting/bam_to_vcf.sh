#!/bin/bash

sample=$1
refname=$2

echo "bam_to_vcf $1 $2"
samtools faidx $refname
samtools mpileup -B -E -C50 -f $refname -u $sample.sorted.bam | bcftools view -c > $sample.vcf
