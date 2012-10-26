use CircleGraph;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
require "subfuncs.pl";
require "circlegraphs.pl";

my ($datafile, $outfile) = 0;
GetOptions ('data=s' => \$datafile,
            'outputfile=s' => \$outfile) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if (($datafile==0) && ($outfile==0)) {
    pod2usage(-msg => "GetOptions failed.", -exitval => 2);
}

my $x = draw_circle_graph($datafile);
$x->set_font("Helvetica", 12, "black");
$x->draw_legend_text;

open OUT, ">", "$outfile.ps" or die "couldn't make output file $outfile.ps";
print OUT $x->output_ps . "\n";
close OUT;

__END__

=head1 NAME

tree_omega

=head1 SYNOPSIS

draw_circle_graph [options]

=head1 OPTIONS
    -datafile:      name of tab-delimited file with positions as first column, graphs as subsequent columns.
	-outputfile:    name of output file.

=head1 DESCRIPTION

=cut
