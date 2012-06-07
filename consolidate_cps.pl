#!/usr/bin/perl 
use strict;
use File::Basename;

my @filelist = @ARGV;

foreach my $fastafile (@filelist) {
	open my $filehandle, "<$fastafile" or die "couldn't open $fastafile";
	my $cp_entry = ">$fastafile\n";
	my $line = readline $filehandle;
	while ($line ne "") {
		$line = readline $filehandle;
		$cp_entry .= "$line";
	}
	print "$cp_entry\n";
}

