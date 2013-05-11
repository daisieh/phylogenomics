#!/usr/bin/perl
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Bio::SeqIO;
use Bio::Seq;
use Bio::Align::Utilities qw(cat);
require "subfuncs.pl";

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my $runline = "running " . basename($0) . " " . join (" ", @ARGV) . "\n";

my ($fastafile, $resultfile, $seq_name, $help) = 0;
GetOptions ('fasta=s' => \$fastafile,
            'outputfile=s' => \$resultfile,
            'sequence=s' => \$seq_name,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

print $runline;

unless (($fastafile && $resultfile)) {
    my $msg = qq{Error: an option was mis-specified:
    fasta=$fastafile
    outputfile=$resultfile
    };

    pod2usage(-msg => $msg, -exitval => 2);
}

open FH, "<", "$fastafile" or die "couldn't open $fastafile";
my $line = readline FH;

while ($line !~ /$seq_name/) {
	$line = readline FH;
}

open RESULT_FH, ">", "$resultfile"."\.fasta";
print RESULT_FH $line;
$line = readline FH;

while ($line !~ />/) {
	print RESULT_FH $line;
	$line = readline FH;
}

close RESULT_FH;
close FH;

__END__

=head1 NAME

slice_fasta_file

=head1 SYNOPSIS

slice_fasta_file [-fasta fa_file] [-genbank gb_file | -slices slice_file | -start -end] [-outputfile output_file]

=head1 OPTIONS

  -fasta:           fasta file of aligned sequences
  -genbank|gb_file: genbank file with CDS coordinates
  -start:           start position of single slice
  -end:             end position of single slice
  -outputfile:      output file name

=head1 DESCRIPTION

Given a fasta file of aligned sequences and a corresponding genbank file
with the CDS coordinates, will create a fasta file with each
CDS corresponding to a separate sequence block.

=cut
