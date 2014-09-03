#!/usr/bin/env perl

use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
require "bioperl_subfuncs.pl";

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my $runline = "running " . basename($0) . " " . join (" ", @ARGV) . "\n";

my ($fastafile, $resultfile, $gb_file, $start, $end, $oneslice, $help, $slice_list, $concat, $multiple) = 0;
my $type = "CDS";
GetOptions ('fasta=s' => \$fastafile,
            'outputfile=s' => \$resultfile,
            'genbank|gb_file:s' => \$gb_file,
            'slices|slicelist:s' => \$slice_list,
            'start:i' => \$start,
            'end:i' => \$end,
            'concatenate' => \$concat,
            'locus|type:s' => \$type,
            'multiple' => \$multiple,
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
my $whole_aln = make_aln_from_fasta_file ($fastafile);

my @alns = ();
my $gene_alns = \@alns;

if ($oneslice) {
	my $curr_slice = $whole_aln->slice($start, $end, 1);
    push @alns, $curr_slice;
} elsif ($gb_file) {
    $gene_alns = slice_fasta_from_genbank_file ($fastafile, $gb_file, $type);
} elsif ($slice_list) {
    open FH, "<", $slice_list or die "Couldn't open slice list $slice_list";
    my @slices = <FH>;
    close FH;

    foreach my $slice (@slices) {
        if ($slice =~ /(.+?)\t(.+?)\t(.+)$/) {
            my $name = $1;
            my $start = $2;
            my $end = $3;
            my $curr_slice = $whole_aln->slice($start, $end, 1);
            $curr_slice->description($name);
            push @alns, $curr_slice;
        } elsif ($slice =~ /^(.+?)-(.+)$/) {
            my $curr_slice = $whole_aln->slice($1, $2, 1);
            push @alns, $curr_slice;
        } else {
            pod2usage(-msg => "Slice file improperly specified.", -exitval => 2);
        }
    }
} else {
    pod2usage(-msg => "Must specify either a genbank file or a start and end point", -exitval => 2);
}

if ($concat) {
    # construct a new SimpleAlign by going seq by seq through gene_alns and adding it to the new concatenated SimpleAlign
    my $concat_aln = new Bio::SimpleAlign();
    my @seqs = @$gene_alns[0]->each_seq();
    for (my $i=0;$i<@seqs;$i++) {
        my $concat_seq = "";
        foreach my $aln (@$gene_alns) {
            $concat_seq .= $aln->get_seq_by_pos($i+1)->seq();
        }
        $concat_aln->add_seq(new Bio::LocatableSeq(-seq => $concat_seq, -id  => $seqs[$i]->id()));
    }
    open FH, ">", $resultfile;
    foreach my $seq ($concat_aln->each_seq()) {
        my $name = $seq->id();
        print FH ">$name\n";
        print FH $seq->seq() . "\n";
    }
	close FH;
} elsif ($multiple) {
    my $i = 1;
    foreach my $aln (@$gene_alns) {
        my $gene_name = $aln->description();
        unless ($gene_name) {
            $gene_name = $i;
            $i++;
        }
        my $filename;
        if (-d $resultfile) {
            $filename = $resultfile . $gene_name . ".fasta";
        } else {
            $filename = $resultfile . "_" . $gene_name . ".fasta";
        }
        open FH, ">", $filename or die "Couldn't open output file $filename";
        foreach my $seq ($aln->each_seq()) {
            my $name = $seq->id();
            print FH ">$name\n";
            print FH $seq->seq() . "\n";
        }
        close FH;
    }
} else {
    open FH, ">", $resultfile or die "Couldn't open output file $resultfile";
    foreach my $aln (@$gene_alns) {
        my $gene_name = $aln->description();
        if ($gene_name) {
            $gene_name = "_".$gene_name;
        }
        foreach my $seq ($aln->each_seq()) {
            my $name = $seq->id() . "$gene_name";
            print FH ">$name\n";
            print FH $seq->seq() . "\n";
        }
    }
    close FH;
}



__END__

=head1 NAME

slice_fasta_file

=head1 SYNOPSIS

slice_fasta_file [-fasta fa_file] [-genbank gb_file [-locus locus_type] | -slices slice_file | -start -end] [-outputfile output_file] [-concatenate]

=head1 OPTIONS

  -fasta:            fasta file of aligned sequences
  -genbank|gb_file:  genbank file with CDS coordinates
  -slices|slicelist: tab-delimited file with CDS coordinates
  -start:            start position of single slice
  -end:              end position of single slice
  -outputfile:       output file name
  -locus|type:       [optional] if present, specifies the type of genbank tag to slice
  -concatenate:      if flag is present, concatenate the slices into one long alignment
  -multiple:         if flag is present, save each slice as a separate file

=head1 DESCRIPTION

Given a fasta file of aligned sequences and a corresponding genbank file
with the CDS coordinates, will create a fasta file with each
CDS corresponding to a separate sequence block.

=cut
