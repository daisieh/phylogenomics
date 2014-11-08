#!/usr/bin/env perl

use strict;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my $help = 0;
my $outfile = "";
my $lookupfile = "";
my $subjectfile = "";
my $match = 0;

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

GetOptions ('lookup=s' => \$lookupfile,
			'subject=s' => \$subjectfile,
			'match=i' => \$match,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help){
    pod2usage(-verbose => 1);
}

my @items_to_find = ();

open FIND_FH, "<", $lookupfile;
foreach my $line (<FIND_FH>) {
	chomp $line;
	push @items_to_find, $line;
}
close FIND_FH;

my $dictionary = {};
open FH, "<", $subjectfile;
foreach my $line (<FH>) {
	my @items = split (/\t/, $line);
	my $key = $items[$match-1];
	$dictionary->{$key} = $line;
}
close FH;
foreach my $item (@items_to_find) {
	print $dictionary->{$item};
}
