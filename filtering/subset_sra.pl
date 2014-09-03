#!/usr/bin/env perl

my $bamfile = shift;
my $frac = shift;

my $totallines = `samtools view -c $bamfile`;
my $fraclines = $totallines * $frac;
$fraclines =~ s/\.*//;

print `samtools view $bamfile | head -n $fraclines`;
