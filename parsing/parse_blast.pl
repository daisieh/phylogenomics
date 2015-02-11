#!/usr/bin/env perl

use Getopt::Long;
use Pod::Usage;
use File::Temp qw (tempfile tempdir);
use FindBin;
use lib "$FindBin::Bin/../lib";
use Subfunctions qw (parse_fasta write_fasta blast_to_genbank);
use Genbank qw (parse_genbank write_sequin_tbl);
use Data::Dumper;

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my $help = 0;
my $outfile = "";
my $gbfile = "";
my $fastafile = "";


GetOptions ('reference=s' => \$gbfile,
			'fastafile=s' => \$fastafile,
			'outfile=s' => \$outfile,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

if ($gbfile !~ /\.gb$/) {
	print "reference file needs to be a fully annotated Genbank file.\n";
	exit;
}

my ($ref_hash, $ref_array) = blast_to_genbank ($gbfile, $fastafile);

my $gene_array = Subfunctions::align_regions_to_reference ($ref_hash, $ref_array, $gbfile);

my (undef, $fastaarray) = parse_fasta($fastafile);
# there should be only one key, so just one name.
my $name = @$fastaarray[0];

open FH, ">", "$outfile.tbl";
print FH write_sequin_tbl ($gene_array, $name);
close FH;

__END__

=head1 NAME

parse_blast

=head1 SYNOPSIS

parse_blast -reference genbank.gb -fasta fastafile [-outputfile output_file]

=head1 OPTIONS

  -fastafile:       fasta sequence to blast
  -reference:       genbank file to blast against
  -outputfile:      name of output file

=head1 DESCRIPTION

Given a Genbank-formatted reference file and a fasta-formatted sequence of a putative homolog
of that reference, outputs a Sequin-formatted feature table with the homolog aligned to the
gene features in the reference file.

=cut
