#!/usr/bin/env perl

use Getopt::Long;
use Pod::Usage;
use File::Temp qw (tempfile tempdir);
use FindBin;
use lib "$FindBin::Bin/../lib";
use Subfunctions qw (parse_fasta write_fasta blast_to_genbank);
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

my ($result_hash, $result_array) = blast_to_genbank ($gbfile, $fastafile, $outfile);

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


=cut
