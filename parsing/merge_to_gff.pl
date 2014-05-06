#!/usr/bin/env perl
use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use File::Spec;
use File::Path qw (make_path);
use FindBin;
use lib "$FindBin::Bin/..";
use Subfunctions qw (split_seq parse_fasta);
use lib "$FindBin::Bin";
use GFF qw (feature_to_seq parse_gff_block parse_attributes export_gff_block read_gff_block write_gff_file set_gff_sequence);

my $gff_file = "";
my $gene = "";
my $fastafile = "";
my $outfile = "";
my $genefile = "";
my $blastfile = "";
my $help = 0;

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

GetOptions ('gfffile=s' => \$gff_file,
			'blastdir=s' => \$blastfile,
			'fastadir=s' => \$fastafile,
			'outfile=s' => \$outfile,
			'genefile=s' => \$genefile,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

my @genes = ();
my $outdir = File::Spec->rel2abs($outfile);
unless (-d $outdir) {
	make_path($outdir);
}
open GENE_FH, "<", $genefile or die "couldn't open gene file.";
while (my $line = readline GENE_FH) {
	if ($line =~ /^(.+?)\s*/) {
		push @genes, $line;
	}
}

my @sorted_genes = sort @genes;

open my $fh, "<", $gff_file;

foreach my $gene (@sorted_genes) {
	my $gff_block = read_gff_block($fh, $gene);

	if ($gff_block eq "") {
		print "No gene named $gene found.\n";
		exit;
	}

	my $gff_hash = parse_gff_block ($gff_block);
	(my $seqhash, undef) = parse_fasta (File::Spec->catfile ($fastafile, "$gene.fasta"));

	set_gff_sequence ($gff_hash, $seqhash->{$gene});
	$gff_hash->{"seqid"} = $gene;
	$gff_hash->{"source"} = "Ser_aTRAM";

	# first clear out all the start and end values within each mRNA:
	for (my $i=1; exists $gff_hash->{"mRNA"}->{$i}; $i++) {
		foreach my $type (keys %{$gff_hash->{"mRNA"}->{$i}}) {
			my $mRNA_hash = $gff_hash->{"mRNA"}->{$i};
			$gff_hash->{"mRNA"}->{$i}->{"start"} = 0;
			$gff_hash->{"mRNA"}->{$i}->{"end"} = 0;
			if ((ref $mRNA_hash->{$type}) =~ /HASH/) {
				for (my $j=1; exists $mRNA_hash->{$type}->{$j}; $j++) {
					$gff_hash->{"mRNA"}->{$i}->{$type}->{$j}->{"start"} = 0;
					$gff_hash->{"mRNA"}->{$i}->{$type}->{$j}->{"end"} = 0;
				}
			}
		}
	}


	open BLASTFH, "<", File::Spec->catfile ($blastfile, $gene);
	foreach my $line (<BLASTFH>) {
		if ($line =~ /$gene\.gene/) {
			next;
		}
# 		Potri.001G000200.1.exon.1	525	608
		my ($longname,$start,$end,undef) = split (/\t/, $line, 4);
		chomp $end;
		if ($longname =~ /$gene\.(\d+)\.(.+?)\.(\d+)/) {
			my $mRNA = $1;
			my $type = $2;
			my $num = $3;
			$gff_hash->{"mRNA"}->{$mRNA}->{$type}->{$num}->{"start"} = $start;
			$gff_hash->{"mRNA"}->{$mRNA}->{$type}->{$num}->{"end"} = $end;
		} elsif ($longname =~ /$gene\.(\d+)$/) {
			print "mRNA $longname\n";
			$gff_hash->{"mRNA"}->{$1}->{"start"} = $start;
			$gff_hash->{"mRNA"}->{$1}->{"end"} = $end;
		}
	}

	write_gff_file ($gff_hash, File::Spec->catfile($outdir,"$gene.gff"));
}

close $fh;

__END__

=head1 NAME

parse_blast

=head1 SYNOPSIS

GetOptions ('gfffile=s' => \$gff_file,
			'blastdir=s' => \$blastfile,
			'fastadir=s' => \$fastafile,
			'outfile=s' => \$outfile,
			'genefile=s' => \$genefile,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

=head1 DESCRIPTION

Parses an "outfmt 3" formatted blastn file to generate a list of regions to be used in
Genbank annotations.

=cut

