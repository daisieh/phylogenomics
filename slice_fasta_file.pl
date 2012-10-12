#!/usr/bin/perl
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Bio::SeqIO;
use Bio::Align::Utilities qw(cat);

require "subfuncs.pl";

my $usage = "perl " . basename($0);
$usage .=	" <fastafile> <cp.fasta> <resultfile>\n\n";
print "running " . basename($0) . " " . join (" ", @ARGV) . "\n";

my ($fastafile, $resultfile, $gb_file, $start, $end, $oneslice) = 0;
GetOptions ('fasta=s' => \$fastafile,
            'outputfile=s' => \$resultfile,
            'genbank|gb_file' => \$gb_file,
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
    my @gene_alns;

    my $seqio_object = Bio::SeqIO->new(-file => $gb_file);
    my $seq_object = $seqio_object->next_seq;
    print $seq_object;

    my $result_str = "";
    while ($seq_object) {
        for my $feat_object ($seq_object->get_SeqFeatures) {
            print $feat_object->primary_tag;
            if ($feat_object->primary_tag eq "CDS") {
                my $name = main_name_for_gb_feature($feat_object);
                my @locations = $feat_object->location->each_Location;
                my $cat_aln = 0;
                my $strand = 0;
                my $last_end = 0;
                foreach $loc (@locations) {
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
                    foreach $seq ( $cat_aln->each_seq() ) {
                        $seq = $seq->revcom();
                        $flipped_aln->add_seq($seq);
                    }
                    $cat_aln = $flipped_aln;
                }

                $cat_aln = $cat_aln->slice(1, $cat_aln->length()-3);
                $cat_aln->description($name);
                push @gene_alns, $cat_aln;
                print "$result_str\t$last_end\n";
                $result_str = "";
            }
        }
        $seq_object = $seqio_object->next_seq;
    }

    #open my $gene_file, ">", "$result_dir\/$gb_file.fa";
    foreach my $aln (@gene_alns) {
        foreach my $seq ( $aln->each_seq()) {
            print ">" . $aln->description . "\n" . $seq->seq() . "\n";
        }
    }
}

__END__

=head1 NAME

tree_omega

=head1 SYNOPSIS

filter_cp [options]

=head1 OPTIONS

    -fasta:     fasta file to filter
    -reference: fasta file of the reference to be filtered against
    -task:      "discard-hits" means to discard hits to the reference
                "keep-hits" means to keep hits
    -evalue:    sets the evalue for blastn (default is 10)
    --blast/noblast:    runs blastn or uses previously generated filtered list (default is to run blastn)

=head1 DESCRIPTION

=cut
