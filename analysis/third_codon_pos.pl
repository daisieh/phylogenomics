#!/usr/bin/env perl

use strict;
use File::Temp qw/ tempfile /;

my $fastafile = shift;

my (undef, $tempfile) = tempfile(UNLINK => 1);
# my $tempfile = "$fastafile.temp";
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
my $location = 0;
my @ns_seqs = ();
my @s_seqs = ();
for (my $i=0; $i<@seqs; $i++) {
	push @ns_seqs, "";
	push @s_seqs, "";
}
print "hi\n";
my $count = 0;
while ($seqs[0]) {
	for (my $i=0; $i<@seqs; $i++) {
		if ($seqs[$i] =~ /^(.{2})(.)(.*)$/) {
			$ns_seqs[$i] .= "$1";
			$s_seqs[$i] .= "$2";
			$seqs[$i] = "$3";
		} else {
			die "not a multiple of 3, $count\n";
# 			$seqs[$i] = ;
		}
	}
	$count++;
}

if ($fastafile =~ /(.*)\.fa.*/) {
	$fastafile = $1;
}
open NS_FH, ">", "$fastafile.ns.fasta";
open S_FH, ">", "$fastafile.s.fasta";
for (my $i=0; $i<@seqs; $i++) {
	print NS_FH ">$names[$i]\n$ns_seqs[$i]\n";
	print S_FH ">$names[$i]\n$s_seqs[$i]\n";
}
close NS_FH;
close S_FH;


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
