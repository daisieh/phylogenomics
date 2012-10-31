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

my ($datafile, $outfile, $gb_file, $help) = 0;

GetOptions ('datafile|inputfile=s' => \$datafile,
            'outputfile=s' => \$outfile,
            'genbank|gb_file:s' => \$gb_file,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}
print $runline;


unless (($datafile) && ($outfile)) {
print "$datafile, $outfile\n";
    pod2usage(-msg => "GetOptions failed.", -exitval => 2);
}

my $circlegraph_obj = draw_circle_graph($datafile);
$circlegraph_obj->draw_legend_text ({"size"=>10, "height"=>15});

if ($gb_file) {
	if ($gb_file =~ /\.gb$/) {
		open FH, ">", "$outfile.genes";
		print FH get_locations_from_genbank_file($gb_file);
		close FH;
	    $gb_file = "$outfile.genes";
	}
	draw_gene_map ($gb_file, $circlegraph_obj);
}

open OUT, ">", "$outfile.ps" or die "couldn't make output file $outfile.ps";
print OUT $circlegraph_obj->output_ps . "\n";
close OUT;

__END__

=head1 NAME

draw_circle_graph

=head1 SYNOPSIS

draw_circle_graph [-data] [-genbank] [-output]

=head1 OPTIONS

  -data:            table to graph, with positions in the first column and values to
                    graph in subsequent columns.
  -genbank|gb_file: (optional) genbank file to make gene map
  -outputfile:      prefix of output files

=head1 DESCRIPTION

Graphs values, given a file with a table of values plotted
along positions on a circle.

Optionally, if a corresponding genbank file is specified,
will draw a gene map along the inside of the circle.

=cut
