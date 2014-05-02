#!/usr/bin/env perl

use strict;
use Getopt::Long;
use Pod::Usage;
use File::Spec qw (rel2abs);
use File::Temp qw (tempfile);
use File::Path qw (make_path);
use FindBin;
use lib "$FindBin::Bin/..";
use Subfunctions qw(parse_fasta);

my $validatefile = shift;
my $contigdir = shift;
my $outdir = shift;

$outdir = File::Spec->rel2abs($outdir);
make_path($outdir);


my $singlecopy_contigs = {};
my $singlecopy_seqs = {};
open VALFH, "<", $validatefile;
foreach my $line (<VALFH>) {
#Potri.001G166200.1	1e-90	457	810	81.18	self	4.22_len_1868_cov_11.1
	if ($line =~ /(.+?)\t.+?\tself\t(.*)$/) {
		$singlecopy_contigs->{$1} = $2;
	}
}
close VALFH;
print "found " . (keys %$singlecopy_contigs) . " single copy genes\n";

foreach my $key (keys %$singlecopy_contigs) {
	my $bestfile = File::Spec->catfile($contigdir, "$key.best.fasta");
	my ($contigs, $contigarray) = parse_fasta($bestfile);
	$singlecopy_seqs->{$key} = $contigs->{$singlecopy_contigs->{$key}};
	print "finding $singlecopy_contigs->{$key} in $key.best.fasta\n";
}

open OUTFH, ">", File::Spec->catfile($outdir, "genelist.txt");
foreach my $key (keys %$singlecopy_seqs) {
	print OUTFH "$key\n";
	my $fh_name = File::Spec->catfile($outdir, "$key.fasta");
	my $aln_name = File::Spec->catfile($outdir, "$key");
	open FH, ">", $fh_name;
	print FH ">$key\n$singlecopy_seqs->{$key}\n";
	close FH;
}
close OUTFH;
