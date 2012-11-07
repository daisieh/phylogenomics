#!/usr/bin/perl
use strict;

my $inputfile = @ARGV[0];

open my $fh, "<", $inputfile or die "couldn't open $inputfile";

open my $outfh, ">", "$inputfile.gb" or die "couldn't open $inputfile.gb";

my $line = readline $fh;
my $header;
$line =~ />Features .+\[organism=(.*)\].*\[location=(.*)\].*\] (.*)$/;
$header = qq{DEFINITION\t$3
SOURCE\t$2
ORGANISM\t$1
FEATURES\tLocation/Qualifiers
};
print $outfh $header;

while ($line = readline $fh) {
	if ($line =~ m/^(\d+?)\t(\d+?)\t(.+)$/) {	# 74	1	gene
		my $start = $1;
		my $stop = $2;
		my $type = $3;
		$type =~ s/\r|\n//g;
		my $region;
		if ($start > $stop) { # complement
			$region = "complement($stop..$start)";
		} else {
			$region = "$start..$stop";
		}
		print $outfh "\t$type\t\t$region\n";
	} elsif ($line =~ m/^\t\t\t(.+?)\t(.+)$/) { # 			gene	trnH-GUG
		my $type = $1;
		my $val = $2;
		$val =~ s/\r|\n//g;
		print $outfh qq{\t\t\t/$type="$val"\n};
	}
}

close $fh;
close $outfh;
