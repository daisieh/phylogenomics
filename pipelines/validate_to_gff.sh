#!/bin/bash

#PBS -l nodes=1:ppn=8
#PBS -l walltime=48:00:00
#PBS -l mem=8000mb
#PBS -j oe

module load perl/5.14.2
module load blast/ncbi-2.2.28+

# $1 is the input directory, $2 is the output directory.
for f in $1
do
echo $f
perl /home/daisieh/aTRAM/Postprocessing/ValidateGenes.pl -ref /home/daisieh/potri_cds_refseqs.fasta -processes 8 -output $2/validate -in $1

# filter single genes.
perl /home/daisieh/phylogenomics/filtering/filter_atram_single_copy_only.pl $2/validate/results.txt $1 $2/single
# output is in $1/single/genelist.txt

# blast list
perl /home/daisieh/phylogenomics/parsing/blast_list.pl -ref /home/daisieh/refs/pop_refs -gene $2/single/genelist.txt -fasta $2/single/ -out $2/blast

# merge gff
perl /home/daisieh/phylogenomics/parsing/merge_to_gff.pl -gff /home/daisieh/refs/Ptrichocarpa_210_gene_exons.gff3 -gene $2/single/genelist.txt -fasta $2/single/ -blast $2/blast -out $2/gff

done
