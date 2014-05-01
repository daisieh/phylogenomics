# Workflow:
# 1) aTRAM basic pipeline.
# 2) validate genes
# 3) filter just single copy genes
# 4) slice corresponding reference gff to fasta for blast
# 5) blast single copies against reference fastas
# 6) merge_to_gff creates gff version using aTRAM best copies and reference gene features

# filter out just the single copies.
# `perl $FindBin::Bin/../filtering/filter_atram_single_copy_only.pl $validatefile $contigdir $outdir`;
#
# my $genefile = "genelist.txt";

# make the fasta reference files to blast:
# GetOptions ('gfffile=s' => \$gff_file,
# 			'fastafile=s' => \$fastafile,
# 			'outfile=s' => \$outfile,
# 			'genefile=s' => \$genefile,


