#!/usr/bin/perl
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Bio::SeqIO;
use Bio::Seq;
use Bio::Align::Utilities qw(cat);
require "subfuncs.pl";

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my $runline = "running " . basename($0) . " " . join (" ", @ARGV) . "\n";

my ($fastafile, $resultfile, $gb_file, $start, $end, $oneslice, $help, $slice_list) = 0;
GetOptions ('fasta=s' => \$fastafile,
            'outputfile=s' => \$resultfile,
            'genbank|gb_file:s' => \$gb_file,
            'slices|slicelist:s' => \$slice_list,
            'start:i' => \$start,
            'end:i' => \$end,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

print $runline;

if ($start && $end) {
	$oneslice = 1;
}

unless (($fastafile && $resultfile)) {
    my $msg = qq{Error: an option was mis-specified:
    fasta=$fastafile
    outputfile=$resultfile
    };

    pod2usage(-msg => $msg, -exitval => 2);
}

if ($oneslice) {
    my $whole_aln = make_aln_from_fasta_file ($fastafile);
	my $curr_slice = $whole_aln->slice($start, $end);

	open my $fh, ">", $resultfile;
	my $aln_out = Bio::AlignIO->new(-fh => $fh, -format => "fasta");
	$aln_out->write_aln($curr_slice);
	close $fh;
} elsif ($gb_file) {
    slice_fasta_to_exons ($fastafile, $gb_file, $resultfile);
} elsif ($slice_list) {
    my $seqio = Bio::SeqIO->new(-file => $fastafile, -format => "fasta");
    my %seqhash;
    while (my $seq = $seqio->next_seq()) {
        $seqhash{$seq->id} = $seq;
    }
    open FH, "<", $slice_list or die "Couldn't open slice list $slice_list";
    print "hashed fasta file\n";
    my @slices = <FH>;
    close FH;
    open FH, ">", $resultfile or die "Couldn't open output file $resultfile";
    foreach my $slice (@slices) {
        if ($slice =~ /(.+?)\t(.+?)\t(.+)$/) {
            my $name = $1;
            my $start = $2;
            my $end = $3;
            my $seq = $seqhash{$name};
            bless $seq, "Bio::Seq";
            my $result = substr($seq->seq,$start,$end-$start);
            print FH ">$name:$start-$end\n$result\n";
        }
    }
} else {
    pod2usage(-msg => "Must specify either a genbank file or a start and end point", -exitval => 2);
}

__END__

=head1 NAME

slice_fasta_file

=head1 SYNOPSIS

slice_fasta_file [-fasta fa_file] [-genbank gb_file | -start -end] [-outputfile output_file]

=head1 OPTIONS

  -fasta:           fasta file of aligned sequences
  -genbank|gb_file: genbank file with CDS coordinates
  -start:           start position of single slice
  -end:             end position of single slice
  -outputfile:      output file name

=head1 DESCRIPTION

Given a fasta file of aligned sequences and a corresponding genbank file
with the CDS coordinates, will create a fasta file with each
CDS corresponding to a separate sequence block.

=cut
