#!/usr/bin/perl

use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/../lib";
use Genbank qw(parse_genbank write_features_as_fasta write_features_as_table);

my $gbfile = shift;

my $gene_array = parse_genbank($gbfile);

open FASTA_FH, ">", "$gbfile.fasta";
print FASTA_FH write_features_as_fasta ($gene_array);
close FASTA_FH;

open FEATURES_FH, ">", "$gbfile.features";
print FEATURES_FH write_features_as_table ($gene_array);
close FEATURES_FH;

__END__

=head1 NAME

genbank.pl

=head1 SYNOPSIS

genbank gbfile

=head1 DESCRIPTION

Parses a genbank file into a fasta-formatted list of sequences and a features table.

=cut

