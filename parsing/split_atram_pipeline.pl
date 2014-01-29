#!/usr/bin/perl
use File::Spec;
require "subfuncs.pl";

my $fastafile = shift;
my $outdir = shift;

my ($seqs, $seqnames) = parse_fasta($fastafile);
$outdir = File::Spec->rel2abs($outdir);

open LISTFH, ">", File::Spec->catfile( $outdir, "fastalist.txt" );

foreach my $t (@$seqnames) {
	my $clean_name = $t;
	$clean_name =~ s/\|.*$//;
	my $seqfile = File::Spec->catfile( $outdir, "$clean_name.fasta" );
	print LISTFH "$clean_name\t$seqfile\n";
	open FH, ">", $seqfile;
	print FH ">$t\n$seqs->{$t}\n";
	close FH;
}

close LISTFH;
