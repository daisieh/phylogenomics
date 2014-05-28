#!/bin/bash

#PBS -l nodes=1:ppn=1
#PBS -l walltime=24:00:00
#PBS -l mem=8000mb
#PBS -j oe
#PBS -l epilogue=/home/daisieh/epilogue.script

refdir=/home/daisieh/refs/pop_refs

for chr in Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11 Chr12 Chr13 Chr14 Chr15 Chr16 Chr17 Chr18 Chr19 ChrT

do
#echo $chr
# filter single genes.
#perl /home/daisieh/phylogenomics/filtering/filter_atram_single_copy_only.pl /home/daisieh/complete_chrs/$chr/validate/results.txt /home/daisieh/complete_chrs/$chr/best /home/daisieh/complete_chrs/$chr/single
# output is in $chr/single/genelist.txt

# blast list
#perl /home/daisieh/phylogenomics/parsing/blast_list.pl -ref $refdir -gene /home/daisieh/complete_chrs/$chr/single/genelist.txt -fasta /home/daisieh/complete_chrs/$chr/single/ -out /home/daisieh/complete_chrs/$chr/blast

# filter list for best hits
# cd /home/daisieh/complete_chrs/$chr/
# perl /home/daisieh/phylogenomics/filtering/cdscutoff.pl /home/daisieh/complete_chrs/$chr/single/genelist.txt /home/daisieh/complete_chrs/$chr/blast/ 90

# slice files for reference
perl /home/daisieh/phylogenomics/parsing/slice_gff_from_fasta.pl -gff /home/daisieh/refs/Ptrichocarpa_210_gene_exons.gff3 -fasta /home/daisieh/refs/Chrs/$chr.fasta -out $refdir -gene /home/daisieh/complete_chrs/$chr/cutoff.90

# merge gff
# perl /home/daisieh/phylogenomics/parsing/merge_to_gff.pl -gff /home/daisieh/refs/Ptrichocarpa_210_gene_exons.gff3 -gene /home/daisieh/complete_chrs/$chr/cutoff.90 -fasta /home/daisieh/complete_chrs/$chr/single/ -blast /home/daisieh/complete_chrs/$chr/blast -out /home/daisieh/complete_chrs/$chr/gff

done

# rename all to Ser.ph:
# FILES=/home/daisieh/complete_chrs/gff/*
# for file in $FILES
# do
# sed 's/Potri/Ser.ph.Potri/g' $file > Ser.ph.$file
# rm -f $file
# done
#
# mv /home/daisieh/complete_chrs/gff/ /home/daisieh/complete_chrs/Ser_Potri_ph

# make list of these gff files:
ls -1 /home/daisieh/complete_chrs/Ser_Potri_ph > /home/daisieh/complete_chrs/gfflist.txt

# set up the pairwise comparison:
while read p
do
echo $p
perl /home/daisieh/phylogenomics/pipelines/pairwise_ser2pop_cds.pl -gff /home/daisieh/complete_chrs/Ser_Potri_ph/$p -ref $refdir -out /home/daisieh/complete_chrs/pairwise
done < /home/daisieh/complete_chrs/gfflist.txt

