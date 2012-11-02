use CircleGraph;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
require "subfuncs.pl";
require "circlegraphs.pl";

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my $runline = "running " . basename($0) . " " . join (" ", @ARGV) . "\n";

my ($fastafile, $datafile, $out_file, $gb_file, $window_size, $help) = 0;

GetOptions ('fasta:s' => \$fastafile,
            'datafile|inputfile:s' => \$datafile,
            'outputfile=s' => \$out_file,
            'genbank|gb_file:s' => \$gb_file,
            'window:i' => \$window_size,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

print $runline;

my $circle_size;


if ($fastafile) {   # if we were given a fasta file, we should create the diffs file.
    unless ($window_size > 0) { pod2usage(-msg => "no window size specified.", -exitval => 2); }

    my $start_pos = 1;
    my $stop_pos = $window_size;
    my $len;
    my $val = 1;

    my $curr_aln = make_aln_from_fasta_file ($fastafile);
    $datafile = "$out_file.diffs";

    open FH, ">", "$datafile" ;
    print FH "pos\t$fastafile\n";
    $val = perc_diff_partition ($curr_aln, $start_pos, $stop_pos);
    my $firstval = $val;
    while ($val >= 0) {
        my $gene_name = "$start_pos";
        if ($val >= 0) {
            print FH "$start_pos\t$val\n";
        }
        $start_pos = $stop_pos;
        $stop_pos = $start_pos + $window_size;
        $val = perc_diff_partition ($curr_aln, $start_pos, $stop_pos);
    }
    $circle_size = $curr_aln->length();
    print FH "$circle_size\t$firstval\n";
    close FH;
}

# using the data file to make a graph
my $circlegraph_obj = draw_circle_graph($datafile);

if ($window_size > 0) {
    $circlegraph_obj->append_to_legend("Sliding window of $window_size bp");
}

$circlegraph_obj->set_font("Helvetica", 12, "black");
$circlegraph_obj->draw_legend_text;


if ($gb_file) {
	if ($gb_file =~ /\.gb$/) {
		open FH, ">", "$out_file.genes";
		print FH get_locations_from_genbank_file($gb_file);
		close FH;
	    $gb_file = "$out_file.genes";
	}
	draw_gene_map ($gb_file, $circlegraph_obj);
}

$circlegraph_obj->output_ps();
open OUT, ">", "$out_file.ps" or die "couldn't make output file $out_file";
print OUT $circlegraph_obj->output_ps . "\n";
close OUT;

__END__

=head1 NAME

sliding_window

=head1 SYNOPSIS

sliding_window [-fasta -window | -data] [-genbank] [-output]

=head1 OPTIONS

  -fasta:           fasta file of aligned sequences
  -window:          if -fasta is specified, size of sliding window
  -data:            previously generated differences file
  -genbank|gb_file: (optional) genbank file to make gene map
  -outputfile:      prefix of output files

=head1 DESCRIPTION

Given a fasta file and a sliding window size, calculates
the percent difference for each sliding window and maps it to a circular graph.

Can graph a previously-generated differences file as well.
Optionally, if a corresponding genbank file is specified,
will draw a gene map along the inside of the circle.

=cut
