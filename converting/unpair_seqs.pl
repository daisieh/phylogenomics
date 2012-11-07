#!/usr/bin/perl
use strict;
use File::Basename;


my $usage = "perl " . basename($0);
$usage .=	" <file.txt> <result>\n\n";

my $fastafile = shift or die "$usage";
my $resultfile = shift or die "$usage";

open my $F, "<$fastafile" or die "couldn't open fasta file";

my $fs = readline $F;
chomp $fs;
my $result1 .= ">$fs\n";
$fs = readline $F;
chomp $fs;
my $result2 .= ">$fs\n";
$fs = readline $F;
$fs = readline $F;

while ($fs) {
	$result1 .= "$fs";
	$fs = readline $F;
	$result2 .= "$fs";
	$fs = readline $F;
	$fs = readline $F;
}
open my $outfile1, ">$resultfile.fasta" or die "couldn't create result file";
truncate $outfile1, 0;

print $outfile1 "$result1\n$result2";

close $F;
close $outfile1;
