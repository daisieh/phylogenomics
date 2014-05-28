files needed:
# Workflow:
# 1) aTRAM basic pipeline.
# 2) validate genes
# 3) filter just single copy genes
# 4) slice corresponding reference gff to fasta for blast
# 5) blast single copies against reference fastas
# 6) merge_to_gff creates gff version using aTRAM best copies and reference gene features


2. validate.sh:

#!/bin/bash

#PBS -l nodes=1:ppn=8
#PBS -l walltime=2:00:00
#PBS -l mem=8000mb
#PBS -j oe
#PBS -l epilogue=/home/daisieh/epilogue.script

for chr in Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11 Chr12 Chr13 Chr14 Chr15 Chr16 Chr17 Chr18 Chr19 ChrT
do
echo $chr
perl /home/daisieh/aTRAM/Postprocessing/ValidateGenes.pl -db /home/daisieh/complete_chrs/validate/potri_cds_refseqs.blastdb -ref /home/daisieh/potri_cds_refseqs.fasta -processes 8 -output /home/daisieh/complete_chrs/$chr/validate -in /home/daisieh/complete_chrs/$chr/best
done

3. filter_single.sh

#!/bin/bash

#PBS -l nodes=1:ppn=1
#PBS -l walltime=24:00:00
#PBS -l mem=8000mb
#PBS -j oe
#PBS -l epilogue=/home/daisieh/epilogue.script

for chr in Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11 Chr12 Chr13 Chr14 Chr15 Chr16 Chr17 Chr18 Chr19 ChrT
do
echo $chr
# filter single genes.
perl /home/daisieh/phylogenomics/filtering/filter_atram_single_copy_only.pl $chr/validate/results.txt $chr/best $chr/single
# output is in $chr/single/genelist.txt

# slice files for reference
perl /home/daisieh/phylogenomics/parsing/slice_gff_from_fasta.pl -gff /home/daisieh/refs/Ptrichocarpa_210_gene_exons.gff3 -fasta /home/daisieh/refs/Chrs/$chr.fasta -out /home/daisieh/refs/pop_refs -gene /home/daisieh/complete_chrs/$chr/single/genelist.txt

# blast list
perl /home/daisieh/phylogenomics/parsing/blast_list.pl -ref /home/daisieh/refs/pop_refs -gene /home/daisieh/complete_chrs/$chr/single/genelist.txt -fasta /home/daisieh/complete_chrs/$chr/single/ -out /home/daisieh/complete_chrs/$chr/blast

# merge gff
perl /home/daisieh/phylogenomics/parsing/merge_to_gff.pl -gff /home/daisieh/refs/Ptrichocarpa_210_gene_exons.gff3 -gene /home/daisieh/complete_chrs/$chr/single/genelist.txt -fasta /home/daisieh/complete_chrs/$chr/single/ -blast /home/daisieh/complete_chrs/$chr/blast -out /home/daisieh/complete_chrs/$chr/gff

done

GetOptions ('gfffile=s' => \$gff_file,
			'blastdir=s' => \$blastfile,
			'fastadir=s' => \$fastafile,
			'outfile=s' => \$outfile,
			'genefile=s' => \$genefile,
