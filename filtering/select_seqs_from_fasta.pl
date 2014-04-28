use strict;
use File::Temp qw/ tempfile tempdir /;
use FindBin;
use lib "$FindBin::Bin/..";
use Subfunctions qw(parse_fasta);

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

my $found = findsequences ($fastafile, \@seqs);

open OUT_FH, ">", "$outfile";
foreach my $seq (@seqs) {
	print OUT_FH ">$seq\n$found->{$seq}\n";
}
close OUT_FH;

sub findsequences {
	my $fastafile = shift;
	my $names = shift;

	unless (-e $fastafile) {
		die "File $fastafile does not exist.\n";
	}

	my $hashed_seqs = {};
	my ($taxa, $taxanames) = parse_fasta ($fastafile);

	foreach my $name (@$names) {
		if (exists $taxa->{$name}) {
			$hashed_seqs->{$name} = $taxa->{$name};
		}
	}
	return $hashed_seqs;
}
