require "subfuncs.pl";

use CircleGraph;

my $usage  = "graph_by_genes.pl input_file output_file\n";
$usage .= qq{columns should be "gene	start	stop	value"\n};
$usage .= qq{last row should have the stop value of the total size of the circle, in bp.\n};

my $datafile = shift or die $usage;
my $out_file = shift or die $usage;

open INPUTFILE, "<$datafile" or die "$datafile failed to open\n";
my @inputs = <INPUTFILE>;
close INPUTFILE;

while (@inputs[0] !~ /\t/) { #there's some sort of header
    shift @inputs;
    if (@inputs == 0) {
        die "no data in input file.";
    }
}

my @sorted = sort (@differences);
my $diff_len = @sorted;
my $max_diffs = @sorted[@sorted-1];
(undef, undef, my $circle_size, undef) = split /\t/, pop @inputs;
$circle_size =~ s/\n//;
my $x = new CircleGraph();

$x->draw_circle($x->outer_radius);

for (my $i = 0; $i < @inputs; $i++) {
	my $line = @inputs[$i];
	my ($label, $start, $stop, $value) = split /\t/, $line;
	$value =~ s/\n//;
	if ($value eq "") {
	    $value = 0;
	}

	my $start_angle = ($start/$circle_size) * 360;
	my $stop_angle = ($stop/$circle_size) * 360;
 	my $radius = $x->inner_radius + 20;

 	$x->set_color_by_percent((1-$value)*100);
	$x->draw_filled_arc ($radius, $start_angle, $stop_angle);

	# label this element
	my $center_angle = ($start_angle + $stop_angle) / 2;
    $x->set_font("Helvetica", 6, "black");
 	$x->circle_label($center_angle, $x->inner_radius + 22, $label);

}

$x->draw_circle($x->inner_radius, {filled => 1, color => "white"});
$x->draw_circle($x->inner_radius);

$x->output_ps();
open OUT, ">", $out_file;
print OUT $x->output_ps . "\n";
close OUT;
