#!/usr/bin/env perl

use strict;
use File::Temp qw/ tempfile /;

my $fastafile = shift;
my $regionsfile = shift;

my (undef, $tempfile) = tempfile(UNLINK => 1);
flattenfasta ($fastafile, $tempfile, "\t");

open FH, "<", $tempfile;
my @seqs = ();
my @names = ();
foreach my $line (<FH>) {
	$line =~ />(.*?)\t(.*)$/;
	push @names, $1;
	push @seqs, $2;
}
close FH;

open FH, "<", $regionsfile;
my $location = 0;
my @newseqs = ();
foreach my $line (<FH>) {
	$line =~ /(.*?)\t(.*?)\t(.*)$/;
	my $start = $2;
	my $end = $3;
	my $region = $1;
	my $regionsize = $end - $start + 1;

	for (my $i=0; $i<@seqs; $i++) {
		my ($a, $b, $c) = split_seq ($seqs[$i],$start,$end);
		$newseqs[$i] .= "$b";
	}
	$location = $end+1;

}

for (my $i=0; $i<@seqs; $i++) {
	print ">$names[$i]\n$newseqs[$i]\n";
}

sub split_seq {
	my $seq = shift;
	my $start = shift;
	my $end = shift;

	my $max = 30000;
	my $seqlen = length ($seq);
	my $startseq = "";
	my $regionseq = "";
	my $endseq = "";

	my $currstart = $start-1;
	my $currend = $end;
	while ($currstart > $max) {
		$seq =~ /^(.{$max})(.*)$/;
		$startseq .= $1;
		$seq = $2;
		$currstart -= $max;
		$currend -= $max;
	}
	if ($currstart > 0) {
		$seq =~ /^(.{$currstart})(.*)$/;
		$startseq .= $1;
		$seq = $2;
		$currstart = 1;
		$currend = $end - (length ($startseq));
	}

	my $regionsize = $end - $start + 1;
	while ($regionsize > $max) {
		$seq =~ /^(.{$max})(.*)$/;
		$regionseq .= $1;
		$seq = $2;
		$currstart -= $max;
		$currend -= $max;
		$regionsize -= $max;
	}
	if ($regionsize > 0) {
		$seq =~ /^(.{$regionsize})(.*)$/;
		$regionseq .= $1;
		$endseq = $2;
	}
	return ($startseq, $regionseq, $endseq);
}

sub flattenfasta {
	my $fastafile = shift;
	my $outfile = shift;
	my $separator = shift;

	unless ($separator) {
		$separator = '\n';
	}

	my (undef, $tempfile) = tempfile(UNLINK => 1);
	system ("gawk '{if (NF==0) next; s = \"\"; for (i=2;i<=NF;i++) s = s\$i; print \$1\",\"s}' RS=\">\" $fastafile | gawk '{print \">\" \$1 \"$separator\" \$2}' FS=\",\" > $outfile");
}
