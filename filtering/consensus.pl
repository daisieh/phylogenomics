use strict;
require "bioperl_subfuncs.pl";

my $input_file = shift;

my $aln = make_aln_from_fasta_file($input_file);
my $consensus = $aln->consensus_string();

print ">$input_file\n$consensus\n";
