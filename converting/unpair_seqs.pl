#!/usr/bin/perl
use strict;
use File::Basename;


my $usage = "perl " . basename($0);
$usage .=	" <interleaved.fasta> <result_prefix>\n\n";

my $fastafile = shift or die "$usage";
my $resultprefix = shift or die "$usage";

open FH, "<", $fastafile or die "couldn't open fasta file";

open OUT1_FH, ">", "$resultprefix.1.fasta" or die "couldn't create result file";
open OUT2_FH, ">", "$resultprefix.2.fasta" or die "couldn't create result file";

my $fs = readline FH;
while ($fs) {
	print OUT1_FH $fs;
	$fs = readline FH;
	print OUT1_FH $fs;
	$fs = readline FH;
	print OUT2_FH $fs;
	$fs = readline FH;
	print OUT2_FH $fs;
	$fs = readline FH;
}

close FH;
close OUT1_FH;
close OUT2_FH;
