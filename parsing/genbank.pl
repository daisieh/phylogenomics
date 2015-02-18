#!/usr/bin/perl

use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/../lib";
use Genbank;

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my $gbfile = shift;
my $outfile = shift;

if (!defined $outfile) {
	$outfile = $gbfile;
}

my $gene_array = Genbank::parse_genbank($gbfile);

open FASTA_FH, ">", "$outfile.fasta";
print FASTA_FH Genbank::write_features_as_fasta ($gene_array);
close FASTA_FH;

open TBL_FH, ">", "$outfile.tbl";
print TBL_FH Genbank::write_sequin_tbl($gene_array, "test");
close TBL_FH;

__END__

=head1 NAME

genbank.pl

=head1 SYNOPSIS

genbank.pl gbfile

=head1 DESCRIPTION

Parses a genbank file into a fasta-formatted list of sequences and a features table.

=cut

