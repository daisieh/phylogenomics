#!/usr/bin/env perl

use strict;
use File::Temp qw/ tempfile /;

my $fastafile = shift;

my (undef, $tempfile) = tempfile(UNLINK => 1);
flattenfasta ($fastafile, $tempfile, "\t");

open FH, "<:crlf", $tempfile;
my @seqs = ();
my @names = ();
foreach my $line (<FH>) {
	$line =~ />(.*?)\t(.*)$/;
	push @names, $1;
	push @seqs, $2;
}
close FH;

my $total = 0;
my $gaps = 0;
my $diffs = 0;
my $ambiguities = 0;
my $length = length ($seqs[0]);
print "seqs are $length long\n";
while ($seqs[0]) {
	my $column = "";
	for (my $i=0; $i<@seqs; $i++) {
		if ($seqs[$i] =~ /^(.)(.*)$/) {
			$column .= uc("$1");
			$seqs[$i] = "$2";
		}
	}

	my $col_gaps = $column =~ tr/\-N//d;
	if ($col_gaps > 0) {
		$gaps++;
		next;
	} else {
		$total++;
	}
	my $col_ambig = $column =~ tr/ACGT//c;
	if ($col_ambig > 0) {
		$ambiguities++;
	}

	$column = join ('', sort (split(//, $column)));
	$column =~ /(.)/;
	my $char = $1;
	$column =~ s/$char+//;
	if ($column ne "") {
		$diffs++;
	}
}

print "file $fastafile:\n";
print "total length: $length\n";
print "total gaps: $gaps\n";
print "total compared: $total\n";
print "total diffs: $diffs\n";
print "total ambiguities: $ambiguities\n";
print "percent divergence: " . (($diffs / $total) * 100) . "\n";

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
