#!/usr/bin/perl
use strict;
use Bio::SeqIO;
use Bio::Align::Utilities qw(cat);

require "subfuncs.pl";

my $inputfile = @ARGV[0];
my $fa_file = @ARGV[1];

my $whole_aln = make_aln_from_fasta_file ($fa_file);

open my $fh, "<", $inputfile or die "couldn't open $inputfile";

my $line = readline $fh;
my $header;
my $curr_slice;
my $start = 0;
my $stop = 0;
my $val;
my $type;
my $flag = 1;
my @gene_alns;
while ($line = readline $fh) {
	if ($line =~ m/^(\d+?)\t(\d+?)\t(.+)$/) {	# 74	1	gene
		$start = $1;
		$stop = $2;
		$type = $3;
		my $region;
		if ($start > $stop) { # complement
			$region = "complement($stop..$start)";
			$curr_slice = $whole_aln->slice($stop, $start);
		} else {
			$region = "$start..$stop";
			$curr_slice = $whole_aln->slice($start, $stop);
		}
		$flag = 1;
	} elsif ($line =~ m/^\t\t\t(.+?)\t(.+)$/) { # 			gene	trnH-GUG
		my $type = $1;
		$val = $2;
		$val =~ s/\r|\n//g;
		$curr_slice->description("$val $start-$stop");
		$flag = 0;

	}
	if ($flag == 0) {
		if ($type =~ m/gene/) {
			push @gene_alns, $curr_slice;
		}
	}
}

foreach my $aln (@gene_alns) {
	foreach my $seq ( $aln->each_seq()) {
		print ">" . $aln->description . "\n" . $seq->seq() . "\n";
	}
}
