#!/usr/bin/env perl
use strict;
use FindBin;
use lib "$FindBin::Bin/..";
use Subfunctions qw(parse_fasta disambiguate_str get_iupac_code consensus_str);

if (@ARGV < 1) {
	die "Usage: consensus2.pl fastafile\n";
}

my $fastafile = shift;

unless (-e $fastafile) {
	die "File $fastafile does not exist.\n";
}

my ($taxa, $taxanames) = parse_fasta ($fastafile);
my @taxarray = ();
foreach my $taxon (@$taxanames) {
	push @taxarray, $taxa->{$taxon};
}

my $result = consensus_str(\@taxarray);
print ">$fastafile\n";
print "$result\n";
