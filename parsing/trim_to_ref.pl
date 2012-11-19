#!/usr/bin/perl
use strict;
use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use Getopt::Long;
use Pod::Usage;
require "subfuncs.pl";

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my $runline = "running " . basename($0) . " " . join (" ", @ARGV) . "\n";

my ($fastafile, $outfile, $reference, $ref_seq) = 0;
my $window_size = 1000;
my $help = 0;
my $keepfiles = 0;

GetOptions ('fasta:s' => \$fastafile,
            'outputfile=s' => \$outfile,
            'reference:s' => \$reference,
            'keepfiles!' => \$keepfiles,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

print $runline;

my $master_alignment = make_aln_from_fasta_file ($fastafile);

if ($reference ne "") {
    $master_alignment = $master_alignment->set_new_reference($reference);
} else {
    $reference = $master_alignment->get_seq_by_pos(1)->id();
}

# look over the reference sequence: find regions that are not gaps
my @regions = ();
$ref_seq = $master_alignment->get_seq_by_pos(1)->seq();
my $curr_pos = 0;
my $position;
my @sequences = ();
foreach my $seqio ($master_alignment->each_seq) {
	my $seq = "";
	push @sequences, \$seq;
}

while ($ref_seq =~ m/-/gc) {
	$position = pos($ref_seq);
	my @reg_pair = ($curr_pos, $position);
	if (!$position) {
		@reg_pair[1] = length ($ref_seq);
		push @regions, \@reg_pair;
		last;
	}
	push @regions, \@reg_pair;
	$ref_seq =~ m/-*/gc;
	$curr_pos = pos($ref_seq);
}
# then, for all of those regions, push the slices into a new SimpleAlign.
foreach my $reg_ptr (@regions) {
	my @reg_pair = @$reg_ptr;
	$curr_pos = @reg_pair[0];
	$position = @reg_pair[1];
	for (my $i=0; $i<@sequences; $i++) {
		${@sequences[$i]} .= substr($master_alignment->get_seq_by_pos($i+1)->seq(),$curr_pos,($position-$curr_pos-1));
	}
}

open FH, ">", $outfile or die "couldn't make output file $outfile";
foreach (my $i=0; $i<@sequences; $i++) {
	my $name = $master_alignment->get_seq_by_pos($i+1)->id();
	print FH ">$name\n".${@sequences[$i]}."\n";
}
close FH;

__END__

=head1 NAME

pairwise_circle_graphs

=head1 SYNOPSIS

trim_to_ref -fasta -output [-reference] [--keepfiles]

=head1 OPTIONS
  -fasta:           fasta file of aligned sequences
  -outputfile:      prefix of output files
  -reference:       optional: name of sequence to be used as reference (default is first seq)
  --keepfiles:		optional: keep temp files (default is no)

=head1 DESCRIPTION

Given a fasta file of aligned sequences, removes gaps from the reference sequence and
trims the rest of the alignment to match.

=cut

