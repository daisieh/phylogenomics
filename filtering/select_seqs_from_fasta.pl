use strict;
use File::Temp qw/ tempfile tempdir /;
use FindBin;
use lib "$FindBin::Bin/..";
use Subfunctions qw(parse_fasta find_sequences);

if (@ARGV < 2) {
	die "Usage: sequenceretrieval.pl fastafile sequencelist\n";
}

my $fastafile = shift;
my $sequencelist = shift;

my $outfile = "$sequencelist.fasta";

unless (-e $sequencelist) {
	die "File $sequencelist does not exist.\n";
}

unless (-e $fastafile) {
	die "File $fastafile does not exist.\n";
}

my @seqs = ();
open LIST_FH, "<", "$sequencelist";
foreach my $seq (<LIST_FH>) {
	chomp $seq;
	push @seqs, $seq;
}
close LIST_FH;

my $found = find_sequences ($fastafile, \@seqs);

open OUT_FH, ">", "$outfile";
foreach my $seq (@seqs) {
	print OUT_FH ">$seq\n$found->{$seq}\n";
}
close OUT_FH;

