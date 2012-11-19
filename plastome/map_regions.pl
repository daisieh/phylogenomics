use CircleGraph;
use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use Getopt::Long;
use Pod::Usage;
require "subfuncs.pl";
require "circlegraphs.pl";


if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my $runline = "running " . basename($0) . " " . join (" ", @ARGV) . "\n";

my ($infile, $outfile, $reference, $ref_seq, $gb_file, $labelfile) = 0;
my $help = 0;
my $keepfiles = 0;

GetOptions ('input=s' => \$infile,
            'outputfile=s' => \$outfile,
            'genbank|gb:s' => \$gb_file,
            'labels:s' => \$labelfile,
            'keepfiles!' => \$keepfiles,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

print $runline;

unless ($outfile) {
    $infile =~ /(.*?)\./;
    $outfile = $1;
    print "No output name specified, using $outfile.\n";
}

#parse infile into separate chromosome files
open FH, "<", $infile;
my @lines = <FH>;
close FH;
shift @lines;
(undef, undef, my $circle_size, undef) = split /\t/, pop @lines;
$circle_size =~ s/\n//;
my @sortedlines = sort(@lines);

my $currfile = @sortedlines[0];
my @files = ();

print "drawing graphs...\n";
my $circlegraph_obj = new CircleGraph();

# draw the gene map
if ($gb_file) {
	if ($gb_file =~ /\.gb$/) {
        my ($fh, $filename) = tempfile(UNLINK => 1);
		print $fh get_locations_from_genbank_file($gb_file);
		close $fh;
	    $gb_file = "$filename";
	}
	$circlegraph_obj->inner_radius($circlegraph_obj->inner_radius + 50);
	draw_gene_map ($gb_file, $circlegraph_obj);
}

my $radius = 60;
my $i=0;
@sortedlines[0] =~ /(.*?)\t/;
my $currname = $1;
my @radii = ($currname);
my $j=0;
$circlegraph_obj->draw_circle($radius + ($i*15)-3);
foreach my $line (@sortedlines) {
    $line =~ /(.*?)\t/;
    if ($1 ne $currname) {
        $currname = $1;
        print "$i $currname\n";
        $i++;
        push @radii, $currname;
        $circlegraph_obj->draw_circle($radius + ($i*15)-3);
        $j=0;
    }
    my @temp = ();
    push @temp, $line;
    print "\t".($j%10)." $line\n";
    draw_region (\@temp,$circlegraph_obj,$i,$radius + ($i*15)+($j%10),$circle_size);
    $j++;
}

for ($i=0; $i<@radii; $i++) {
    $circlegraph_obj->append_to_legend(@radii[$i],$i);
}

$circlegraph_obj->set_font("Times-Roman", 10, "black");
$circlegraph_obj->draw_legend_text;

open OUT, ">", $outfile.".ps" or die "couldn't make output file $outfile.ps";
print OUT $circlegraph_obj->output_ps . "\n";
close OUT;

if (!$keepfiles) {
    foreach my $file (@files) {
        unlink $file or warn "could not delete $file";
    }
}

sub draw_region {
    my $inputptr = shift;
    my $circlegraph_obj = shift;
    my $color = shift;
    my $radius = shift;
    my $circle_size = shift;

    my @inputs = @$inputptr;

    unless ($circlegraph_obj) {
        $circlegraph_obj = new CircleGraph();
    }

    if ($color) {
        $circlegraph_obj->set_color($color);
    } else {
        $circlegraph_obj->set_color([0, 0, 0.8]);
    }

    unless ($radius) {
        $radius = $circlegraph_obj->inner_radius;
    }

    my @labels = ();
    for (my $i = 0; $i < @inputs; $i++) {
        my $line = @inputs[$i];
        my ($label, $start, $stop, $value) = split /\t/, $line;
        $value =~ s/\n//;
        if ($value eq "") {
            $value = 0;
        }

        my $start_angle = ($start/$circle_size) * 360;
        my $stop_angle = ($stop/$circle_size) * 360;

        $circlegraph_obj->draw_arc ($radius, $start_angle, $stop_angle, {color => $color, width => 3});
    }

    return $circlegraph_obj;
}

__END__

=head1 NAME

map_regions

=head1 SYNOPSIS

map_regions -input -output [-genbank] [-labels] [--keepfiles]
GetOptions ('input=s' => \$infile,
            'outputfile=s' => \$outfile,
            'genbank|gb:s' => \$gb_file,
            'labels:s' => \$labelfile,
            'keepfiles!' => \$keepfiles,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

=head1 OPTIONS
  -input:           list of chromosomes and regions to be mapped along the plastome
  -outputfile:      prefix of output files
  -genbank|gb:      optional: genbank file to generate a map along the graph
  -labels:          optional: tab-delimited list of samples and their corresponding labels
  --keepfiles:      optional: keep temporary files (default is no)

=head1 DESCRIPTION

Given a list of regions to map along the plastome, maps them according to chromosome group.

=cut

