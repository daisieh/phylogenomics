#!/usr/bin/perl
use strict;
use File::Basename;


my $usage = "perl " . basename($0);
$usage .=	" <file.fastq> <result>\n\n";

my $fastafile = shift or die "$usage";
my $resultfile = shift or die "$usage";

my $result1 = "";
my $result2 = "";
my $ref_seq = "";
my $second_seq = "";
my $ref_taxon = "";
my $second_taxon = "";

open my $F, "<$fastafile" or die "couldn't open fasta file";
my $fs = readline $F;
my $ref_taxon = "$fs";
$fs = readline $F;

while ($fs !~ m/>/) {
	chomp $fs;
	$ref_seq .= "$fs";
	$fs = readline $F;
}

$second_taxon = "$fs";
$fs = readline $F;

while ($fs) {
	chomp $fs;
	$second_seq .= "$fs";
	$fs = readline $F;
}

print "$ref_taxon\n";
$ref_seq =~ m/(-+)/g;
my $x = $1;
my $start_pos = 0;
my $stop_pos = (pos $ref_seq)-(length $x);
open my $outfile1, ">$resultfile.fasta" or die "couldn't create result file";
truncate $outfile1, 0;


while (pos $ref_seq) {
	print "$start_pos $stop_pos\n";
	$result1 .= substr ($ref_seq, $start_pos, ($stop_pos - $start_pos));
	$result2 .= substr ($second_seq, $start_pos, ($stop_pos - $start_pos));
# 	$result1 .= "\n";
# 	$result2 .= "\n";

	$start_pos = pos $ref_seq;
	$x = $1;
	$ref_seq =~ m/(-+)/g;
	$stop_pos = (pos $ref_seq)-(length $1);
}

	$start_pos = $start_pos + (length $x);
	$stop_pos = length $ref_seq;
	print "#$start_pos $stop_pos\n";
# 	$result1 .= "\n";
# 	$result2 .= "\n";
	$result1 .= substr ($ref_seq, $start_pos, ($stop_pos - $start_pos));
	$result2 .= substr ($second_seq, $start_pos, ($stop_pos - $start_pos));

print "done\n";

close $F;
print $outfile1 "$ref_taxon$result1\n$second_taxon$result2\n";
close $outfile1;
