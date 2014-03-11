#!/usr/bin/env perl
use strict;
use File::Basename;
use File::Temp qw(tempfile);
use Getopt::Long;
use Pod::Usage;

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my $reffile = "";
my $outputfile = "";
my $bamfile = "";
my $help = 0;
my $max_reads = 0;

GetOptions ('reference=s' => \$reffile,
			'bamfile=s' => \$bamfile,
            'output=s' => \$outputfile,
            'maxreads=i' => \$max_reads,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

if ($reffile !~ /\.fa/) {
    pod2usage(-msg => "Must specify a fasta reference file.");
}
if ($bamfile !~ /\.bam/) {
    pod2usage(-msg => "Must specify a bam archive.");
}
if ($outputfile eq "") {
    pod2usage(-msg => "Must specify a destination path.");
}

`which bwa`;
if (($? >> 8) != 0) {
	die "Couldn't find bwa.";
}

`which samtools`;
if (($? >> 8) != 0) {
	die "Couldn't find samtools.";
}

my $cmd = "";
# make the output file an absolute path, just to be safe.
$outputfile = File::Spec->rel2abs($outputfile);
my $output_path = dirname ($outputfile);

# check to see if the path for the outputfile exists; if not, create it.
unless (-d $output_path) {
	make_path ($output_path);
}
my $i = 0;
print $i++ . "\n";
my $bam_input = $bamfile;
# make a bam input file with the corresponding number of reads, if requested.
if ($max_reads > 0) {
	open FH, ">", "$outputfile.bam";
	close FH;
	$cmd = "samtools view -h $bamfile | head -n $max_reads | samtools view -bhS - > $bam_input";
	print $i++ . ": $cmd\n";
	`$cmd`;
}

# index the reference
unless (-e "$reffile.bwt") {
	$cmd = "bwa index $reffile";
	print $i++ . ": $cmd\n";
	`$cmd`;
}

# bwa align the two paired ends
(undef, my $sai1) = tempfile ();
$cmd = "bwa aln -b1 $reffile $bam_input > $sai1";
print $i++ . ": $cmd\n";
`$cmd`;

(undef, my $sai2) = tempfile ();
$cmd = "bwa aln -b2 $reffile $bam_input > $sai2";
print $i++ . ": $cmd\n";
`$cmd`;

# generate paired end alignments
(undef, my $sam) = tempfile ();
$cmd = "bwa sampe $reffile $sai1 $sai2 $bam_input $bam_input > $sam";
print $i++ . ": $cmd\n";
`$cmd`;


# make bam output
$cmd = "samtools view -Sub $sam > $outputfile";
print $i++ . ": $cmd\n";
`$cmd`;


__END__

=head1 NAME

pileup_bam_to_ref.pl

=head1 SYNOPSIS

pileup_bam_to_ref.pl -ref reference.fasta -bam file.bam -out outfile [-max maxreads]

Aligns paired end reads in BAM format to a reference sequence in fasta format.

=head1 OPTIONS

 -ref:      reference to align bam reads to.
 -bam:      bam file with paired reads to align.
 -output:   name of output file.
 -maxreads: optional: if the read file is very large, can specify max number of reads to pileup.

=cut

