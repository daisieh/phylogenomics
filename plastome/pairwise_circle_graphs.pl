#!/usr/bin/env perl
use CircleGraph;
use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use Getopt::Long;
use Pod::Usage;
require "bioperl_subfuncs.pl";
require "circlegraphs.pl";
use FindBin;
use lib "$FindBin::Bin/../lib";
use Subfunctions;

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my $runline = "running " . basename($0) . " " . join (" ", @ARGV) . "\n";

my ($fastafile, $outfile, $reference, $ref_seq, $gb_file, $labelfile) = 0;
my $window_size = 1000;
my $help = 0;
my $keepfiles = 0;

GetOptions ('fasta:s' => \$fastafile,
            'outputfile=s' => \$outfile,
            'window:i' => \$window_size,
            'reference:s' => \$reference,
            'genbank|gb:s' => \$gb_file,
            'labels:s' => \$labelfile,
            'keepfiles!' => \$keepfiles,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

print $runline;

unless ($outfile) {
    $fastafile =~ /(.*?)\.fa(sta)*/;
    $outfile = $1;
    print "No output name specified, using $outfile.\n";
}

my $master_alignment = make_aln_from_fasta_file ($fastafile);

if ($reference ne "") {
    $master_alignment = $master_alignment->set_new_reference($reference);
} else {
    $reference = $master_alignment->get_seq_by_pos(1)->id();
}

my @tempfiles;

for (my $i=2; $i<=$master_alignment->num_sequences(); $i++) {
    my $comp_seq = $master_alignment->get_seq_by_pos($i);
    my $comp_aln = $master_alignment->select_noncont(1,$i);
    print $comp_seq->id() . "\n";
    my $start_pos = 1;
    my $stop_pos = $window_size;
    my $filename = $outfile.".".$comp_seq->id().".diffs";
    push @tempfiles, $filename;
    open FH, ">", $filename;
    print FH "pos\t".$comp_seq->id()."\n";
    my $val = perc_diff_partition ($comp_aln, $start_pos, $stop_pos, 1);
    my $firstval = $val;
    while ($val >= -100) {
        if ($val >= -100) {
            print FH "$start_pos\t$val\n";
        }
        $start_pos = $stop_pos;
        $stop_pos = $start_pos + $window_size;
        $val = perc_diff_partition ($comp_aln, $start_pos, $stop_pos, 1);
    }
    my $circle_size = $comp_aln->length();
    print FH "$circle_size\t$firstval\n";
    close FH;
}

my $diff_matrix = combine_files(\@tempfiles, 1, 1);

my $filename = $outfile.".total.diffs";
push @tempfiles, $filename;
open FH, ">", $filename;
print FH $diff_matrix;
close FH;

print "drawing graphs...\n";
my $circlegraph_obj = draw_circle_graph($filename, 0);

# draw the gene map
if ($gb_file) {
	if ($gb_file =~ /\.gb$/) {
        my ($fh, $filename) = tempfile(UNLINK => 1);
		print $fh get_locations_from_genbank_file($gb_file);
		close $fh;
	    $gb_file = "$filename";
	}
	draw_gene_map ($gb_file, $circlegraph_obj);
}

$circlegraph_obj->append_to_legend("Pairwise comparisons to $reference");
my $newlegend = $circlegraph_obj->legend;

$circlegraph_obj->draw_legend_text({size => 10, height => 15});

open OUT, ">", $outfile.".ps" or die "couldn't make output file $outfile.ps";
print OUT $circlegraph_obj->output_ps . "\n";
close OUT;

if (!$keepfiles) {
    foreach my $file (@tempfiles) {
        unlink $file or warn "could not delete $file";
    }
}

__END__

=head1 NAME

pairwise_circle_graphs

=head1 SYNOPSIS

pairwise_circle_graphs -fasta -window -output [-reference] [-genbank] [-labels]

=head1 OPTIONS
  -fasta:           fasta file of aligned sequences
  -window:          size of sliding window (if not specified, 1000)
  -outputfile:      prefix of output files
  -reference:       optional: name of sequence to be used as reference (default is first seq)
  -genbank|gb:      optional: genbank file to generate a map along the graph
  -labels:          optional: tab-delimited list of samples and their corresponding labels

=head1 DESCRIPTION

Given a fasta file of aligned sequences, draws graphs comparing them pairwise to the reference seq.

=cut

