#!/usr/bin/perl
use strict;
use File::Basename;
use Getopt::Long;

my ($filelistname, $namefile, $out_file) = 0;

GetOptions ('files|input=s' => \$filelistname,
            'names:s' => \$namefile,
            'outputfile:s' => \$out_file) or die "options misspecified";

open FH, "<:crlf", $filelistname or die "couldn't open $filelistname";
my @filelist = <FH>;
close FH;

my @namelist = @filelist;

if ($namefile) {
    open FH, "<:crlf", $namefile or die "couldn't open $namefile";
    @namelist = <FH>;
    close FH;
}

if (scalar @namelist != scalar @filelist) {
    die "number of names doesn't match number of files.";
}

my $result = "";
for (my $i; $i < @filelist; $i++) {
    my $fastafile = @filelist[$i];
	open my $filehandle, "<$fastafile" or die "couldn't open $fastafile";
	my $cp_entry = ">@namelist[$i]";
	my $line = readline $filehandle;
	while ($line ne "") {
		$line = readline $filehandle;
		$cp_entry .= "$line";
	}
	$result .= "$cp_entry\n";
}

if ($out_file) {
    open OUTFH, ">", $out_file;
    print OUTFH $result;
    close OUTFH;
} else {
    print $result;
}
