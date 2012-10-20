#!/usr/bin/perl
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Bio::SeqIO;
use Bio::Align::Utilities qw(cat);

require "subfuncs.pl";

print "running " . basename($0) . " " . join (" ", @ARGV) . "\n";

my ($fastafile, $resultfile, $gb_file, $start, $end, $oneslice) = 0;
GetOptions ('fasta=s' => \$fastafile,
            'outputfile=s' => \$resultfile,
            'genbank|gb_file=s' => \$gb_file,
            'start:i' => \$start,
            'end:i' => \$end) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

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
my $cat_aln = 0;
my $result_str = "";

if ($oneslice) {
	my $curr_slice = $whole_aln->slice($start, $end);

	open my $fh, ">$resultfile";
	my $aln_out = Bio::AlignIO->new(-fh => $fh, -format => "fasta");
	$aln_out->write_aln($curr_slice);
	close $fh;
} elsif ($gb_file) {
    print "$gb_file\n";
    my @gene_alns;

    my $seqio_object = Bio::SeqIO->new(-file => $gb_file);
    my $seq_object = $seqio_object->next_seq;

    my $result_str = "";
    while ($seq_object) {
        for my $feat_object ($seq_object->get_SeqFeatures) {
            if ($feat_object->primary_tag eq "CDS") {
                my $name = main_name_for_gb_feature($feat_object);
                my @locations = $feat_object->location->each_Location;
                my $cat_aln = 0;
                my $strand = 0;
                my $last_end = 0;
                foreach my $loc (@locations) {
                    $strand = $loc->strand;
                    my $start = $loc->start;
                    my $end = $loc->end;
                    $last_end = $end;
                    my $curr_slice = $whole_aln->slice($start, $end);
                    if ($cat_aln == 0) {
                        $cat_aln = $curr_slice;
                    } else {
                        $cat_aln = cat($cat_aln, $curr_slice);
                    }
                    if ($result_str eq "") { $result_str = "$name\t$start"; }
                }
                if ($strand < 0) {
                    # must flip each seq in the curr_slice
                    my $flipped_aln = Bio::SimpleAlign->new();
                    foreach my $seq ( $cat_aln->each_seq() ) {
                        $seq = $seq->revcom();
                        $flipped_aln->add_seq($seq);
                    }
                    $cat_aln = $flipped_aln;
                }

                $cat_aln = $cat_aln->slice(1, $cat_aln->length()-3);
                $cat_aln->description($name);
                push @gene_alns, $cat_aln;
                $result_str = "";
            }
        }
        $seq_object = $seqio_object->next_seq;
    }

    foreach my $aln (@gene_alns) {
        my $gene = $aln->description;
        open geneFH, ">", "$resultfile\/$gene.fa";
        foreach my $seq ( $aln->each_seq()) {
            print geneFH ">" . $seq->id() . "\n" . $seq->seq() . "\n";
        }
        close geneFH;
    }
}

__END__

=head1 NAME

tree_omega

=head1 SYNOPSIS

slice_fasta_file [options]

=head1 OPTIONS

    -fasta:     fasta file of aligned sequences
	-outputfile:    name of output file or directory for output files (if using -genbank)
	-genbank|gb_file:	genbank file specifying genes
	-start:	start position of single slice
	-end:   end position of single slice

=head1 DESCRIPTION

=cut
