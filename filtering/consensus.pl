#!/usr/bin/env perl
use strict;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use FindBin;
use lib "$FindBin::Bin/../lib";
use Subfunctions qw(parse_fasta pad_seq_ends debug set_debug consensus_str);

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my $inputfile = "";
my $outname = "";
my $help = 0;
my $debug = 0;

if (@ARGV == 1) {
	$inputfile = shift @ARGV;
} else {
	GetOptions ('fasta|input=s' => \$inputfile,
				'outputfile=s' => \$outname,
				'debug' => \$debug,
				'help' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);
}
if ($help) {
    pod2usage(-verbose => 1);
}

set_debug ($debug);

my ($seqmatrix, $seqids)  = parse_fasta ( $inputfile );
delete $seqmatrix->{"length"};

my @rows = values(%$seqmatrix);


my $final_str = consensus_str(\@rows);

if ($outname eq "") {
	print "$final_str\n";
} else {
	open FH, ">", "$outname.fasta";
		print FH ">consensus\n$final_str\n";
	close FH;
}

# recursive function:
# if the number of Ns is less than $max_col_ambigs, this group is fine, return the array.
# if the matrix is a single column with more than $max_col_ambigs, return 0.
# else divide the matrix into halves, recurse on both halves, merge results.


__END__

=head1 NAME

consensus

=head1 SYNOPSIS

consensus -input fastafile -output outputfile

=head1 OPTIONS

  -fasta|input:     fasta file of aligned sequences.
  -outputfile:      output file name.

=head1 DESCRIPTION

returns the consensus sequence of an array of sequences.

=cut

