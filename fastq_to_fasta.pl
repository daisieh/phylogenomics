#!/usr/bin/perl 
use strict;
use File::Basename;

my $usage = "perl " . basename($0);
$usage .=	" <fastqfile>\n\n";

my $fastqfile = shift or die "$usage";

open my $F, "<$fastqfile" or die "couldn't open fastq file";
my $fs = readline $F;

$fastqfile =~ /(.*)\.[fastq|fq]/;
if ($1 ne "") {
	$fastqfile = $1;
}

open my $fastafile, ">$fastqfile.fasta" or die "couldn't open file for writing";
truncate $fastafile, 0;

while ($fs ne "") {	
	#first line:
	my $first_header = "$fs";
	$first_header =~ s/@/>/;
	$fs = readline $F;

	#second line:
	my $sequence = "$fs";
	$fs = readline $F;

	#third line:
	$fs = readline $F;
	
	#fourth line:
	$fs = readline $F;
	
	my $result = "$first_header$sequence\n";
	print $fastafile $result;
}

close $fastafile;
close $F;
