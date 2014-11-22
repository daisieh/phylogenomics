#!/usr/bin/env perl
use strict;
use File::Temp qw/ tempfile tempdir /;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use FindBin;
use lib "$FindBin::Bin/../lib";
use Subfunctions qw(split_seq reverse_complement meld_matrices parse_fasta trim_to_ref align_to_ref);
if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my $runline = "running " . basename($0) . " " . join (" ", @ARGV) . "\n";

my $ref_file = 0;
my $align_file = 0;
my $out_file = 0;
my $help = 0;
my $ref_out = 0;
my $align = 0;

GetOptions ('fasta|input=s' => \$align_file,
            'outputfile=s' => \$out_file,
            'reference=s' => \$ref_file,
            'include_ref' => \$ref_out,
            'align!' => \$align,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

unless($ref_file and $align_file and $out_file) {
    pod2usage(-msg => "Must specify reference file, alignment file, and output file.");
}

print $runline;

my (undef, $alnfile) = tempfile(UNLINK => 1);
my $aligned;
if ($align) {
	$aligned = align_to_ref ($ref_file, $align_file, $alnfile, "mafft");
} else {
	($aligned, undef) = parse_fasta ($align_file);
	$aligned->{"reference"} = delete $aligned->{$ref_file};
	delete $aligned->{"length"};
}

trim_to_ref ($aligned, "reference");

open OUT_FH, ">", $out_file;

my $ref = delete $aligned->{"reference"};

if ($ref_out == 1) {
	print OUT_FH ">reference\n$ref\n";
}

foreach my $seq (keys %$aligned) {
	print OUT_FH ">$seq\n$aligned->{$seq}\n";
}

close OUT_FH;

__END__

=head1 NAME

trim_to_ref

=head1 SYNOPSIS

trim_to_ref -fasta fastafile -reference reffile -output outputfile

=head1 OPTIONS

  -fasta|input:     fasta file of aligned sequences.
  -reference:       fasta file with sequences of interest.
  -outputfile:      output file name.
  -include_ref:     optional: if specified, includes the reference sequence in the fasta output.

=head1 DESCRIPTION

Takes a fasta file and finds aligned regions in each sequence in the fasta file that
match the reference sequence(es). Returns a fasta file of aligned regions of similarity.

=cut

