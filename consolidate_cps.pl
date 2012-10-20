#!/usr/bin/perl
use strict;
use File::Basename;

my $filelistname = shift;
open FH, "<", $filelistname or die "couldn't open $filelistname";
my @filelist = <FH>;
close FH;

foreach my $fastafile (@filelist) {
	open my $filehandle, "<$fastafile" or die "couldn't open $fastafile";
	my $cp_entry = ">$fastafile";
	my $line = readline $filehandle;
	while ($line ne "") {
		$line = readline $filehandle;
		$cp_entry .= "$line";
	}
	print "$cp_entry\n";
}

