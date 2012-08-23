require "subfuncs.pl";
use PostScript::Simple;
use constant CENTER_X => 600; # X coordinate of circle center
use constant CENTER_Y => 600; # Y coordinate of circle center
use constant PS_X_SIZE => 1200; # X size of the PostScript object
use constant PS_Y_SIZE => 1200; # Y size of the PostScript object

my $usage  = "graph_by_genes.pl input_file output_file\n";
$usage .= qq{columns should be "gene	start	stop	value"\n};
$usage .= qq{last row should have the stop value of the total size of the circle, in bp.\n};

my $datafile = shift or die $usage;
my $out_file = shift or die $usage;
my $p = new PostScript::Simple(	xsize => PS_X_SIZE,
								ysize => PS_Y_SIZE,
                                colour => 1,
                                eps => 0,
                                units => "bp");

$p->newpage;
my $OUTER_RADIUS = 387.5;
my $INNER_RADIUS = 337.5;

open INPUTFILE, "<$datafile" or die "$datafile failed to open\n";
my @labels, @starts, @stops, @values;
my $max_diffs = 0;

my @inputs = <INPUTFILE>;
my $total_elems = scalar @inputs;
close INPUTFILE;

my @sorted = sort (@differences);
my $diff_len = @sorted;
$max_diffs = @sorted[@sorted-1];
(undef, undef, my $circle_size, undef) = split /\t/, @inputs[$total_elems-1];

$p->setcolour(black);
$p->setlinewidth(1);
$p->circle(CENTER_X,CENTER_Y, $OUTER_RADIUS);
$p->setlinewidth(2);

for (my $i = 0; $i < $total_elems; $i++) {
	my $line = @inputs[$i];
	my ($label, $start, $stop, $value) = split /\t/, $line;
	chomp $value;

	my $start_angle = ($start/$circle_size) * 360;
	my $stop_angle = ($stop/$circle_size) * 360;
 	#my $radius = $INNER_RADIUS + (($OUTER_RADIUS-$INNER_RADIUS)*(1-$value));
 	my $radius = $OUTER_RADIUS - 5;
 	my @start_coords = coords_on_circle($start_angle, $radius);
 	my @stop_coords = coords_on_circle($stop_angle, $radius);

 	set_percent_red($p, (1-$value)*100);
	draw_filled_arc ($radius, $start_angle, $stop_angle);

	# label this element
	my $center_angle = ($start_angle + $stop_angle) / 2;
 	my @new_coords = coords_on_circle($center_angle,$OUTER_RADIUS + 10);

	$p->setfont("Helvetica", 6);
	$p->setcolour(black);

 	$p->text( {rotate => $center_angle}, @new_coords[0], @new_coords[1], "$label");
}

$p->setcolour("white");
$p->circle({filled => 1}, CENTER_X,CENTER_Y, $INNER_RADIUS);
$p->setcolour(black);
$p->setlinewidth(1);
$p->circle(CENTER_X,CENTER_Y, $INNER_RADIUS);

$p->output("$out_file.ps");

sub draw_filled_arc {
	my $radius = shift;
	my $start_angle = shift;
	my $stop_angle = shift;
	my $center_angle = ($start_angle + $stop_angle) / 2;

 	my @start_coords = coords_on_circle($start_angle, $radius);
 	my @stop_coords = coords_on_circle($stop_angle, $radius);
 	my @center_coords = coords_on_circle($center_angle, $radius);

	$p->arc({filled => 1}, CENTER_X,CENTER_Y,$radius, $start_angle, $stop_angle);
 	$p->polygon({filled => 1}, CENTER_X,CENTER_Y, @start_coords[0], @start_coords[1], @center_coords[0], @center_coords[1], @stop_coords[0], @stop_coords[1]);
}

sub set_percent_red {
	my $p = shift;
	my $percent_red = shift;
	my $scaling = ($percent_red/100);

	my @zero_red = (255,240,240);
	my @full_red = (204,0,0);
	my $r = int(((@full_red[0]-@zero_red[0])*$scaling) + @zero_red[0]);
	my $g = int(((@full_red[1]-@zero_red[1])*$scaling) + @zero_red[1]);
	my $b = int(((@full_red[2]-@zero_red[2])*$scaling) + @zero_red[2]);
	$p->setcolour($r, $g, $b);
}
