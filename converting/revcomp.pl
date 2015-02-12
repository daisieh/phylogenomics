#!/usr/bin/env perl

use strict;
use FindBin;
use lib "$FindBin::Bin/../lib";
use Subfunctions qw (reverse_complement);
use File::Temp qw (tempfile tempdir);

my $input_file = shift;
if (-e $input_file) {
	open FH, "<:crlf", $input_file;

	# my (undef, $tempfile) = tempfile(UNLINK => 1);

	# system ("gawk '{if (NF == 0) next; s = \"\"; for (i=2;i<=NF;i++) s = s\$i; print ">"\$1\"\\n\"s}' RS=\">\" $input_file > $tempfile");
	my $seq = "";
	my $name = readline FH;
	$name =~ s/>//;
	chomp $name;
	while (my $line = readline FH) {
		chomp $line;
		if ($line =~ />(.*)$/) {
			if ($seq ne "") {
				print ">$name\n".reverse_complement($seq)."\n";
				$name = $1;
			}
		} else {
			$seq .= $line;
		}
	}
	print ">$name\n".reverse_complement($seq)."\n";
} else {
	print reverse_complement($input_file) . "\n";
}
