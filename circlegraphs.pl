use PostScript::Simple;
use CircleGraph;
use Bio::AlignIO;
use Bio::SeqIO;

use constant PI 	=> 3.1415926535897932384626433832795;
use constant CENTER_X => 600; # X coordinate of circle center500
use constant CENTER_Y => 600; # Y coordinate of circle center
use constant PS_X_SIZE => 1200; # X size of the PostScript object1000
use constant PS_Y_SIZE => 1200; # Y size of the PostScript object

#
# graphics default values
#
my $MAIN_RADIUS				= 387.5;
my @colours = (	[255,0,0],		# red
				[0,255,0],		# green
				[0,0,255],		# blue
				[200,200,0],	# yellow-ish
				[255,0,255],	# magenta
				[0,255,255],	# cyan
				[255,100,0],	# orange
				[200,50,200],	# fuchsia
				[70,150,0],		# dark green
				[100,0,255]		# violet
				);

sub coords_on_circle {
	my $angle = shift;
	my $radius = shift;

	return ( CENTER_X + ($radius*cos(($angle * PI)/180)), CENTER_Y + ($radius*sin(($angle * PI)/180)));
}

sub draw_circle_graph_from_file {
	my $datafile = shift;
	my $graphfile = shift;
	my $p = new PostScript::Simple( colour => 1, eps => 0, units => "bp", xsize => PS_X_SIZE, ysize => PS_Y_SIZE );
	my $OUTER_RADIUS = 387.5;
	my $INNER_RADIUS = 337.5;

	open my $F, "<$datafile" or die "$datafile failed to open\n";
	my $line = readline $F;

	my @labels = split /\t/, $line;
	my $num_graphs = @labels-1;
	print "drawing $num_graphs graphs\n";

	# print legend
	$p->setfont("Helvetica", 12);
	if (@labels[1] =~ m/.*[:alpha:].*/) {
		for (my $i = 1; $i <= $num_graphs; $i++) {
			print ">@labels[$i]";
			$p->setcolour($colours[$i-1][0],$colours[$i-1][1],$colours[$i-1][2]);
			my $max_height = ($num_graphs * 15) + 60;
			$p->text(10,($max_height - (15*$i)), "@labels[$i]");
		}
		$line = readline $F;
	}

	my $max_diffs = 0;
	my $total_elems;
	my @positions, @differences;
	while ($line ne "") {
		my @items = split ('\t', $line);
		my $pos = shift @items;
		print scalar @items . "\n";
		$total_elems = push (@positions, $pos);
		push (@differences, @items);
		$line = readline $F;
	}

	my @sorted = sort (@differences);
	my $diff_len = @sorted;
	$max_diffs = @sorted[@sorted-1];
	my $window_size = @positions[1]-@positions[0];
	my $circle_size = @positions[$total_elems-1];

	$p->setlinewidth(2);
	for (my $j = 0; $j < $num_graphs; $j++) {
		my @coords = coords_on_circle(0,$INNER_RADIUS);
		my ($last_x, $last_y, $this_x, $this_y);
		$last_x = @coords[0];
		$last_y = @coords[1];
 		$p->setcolour($colours[$j][0],$colours[$j][1],$colours[$j][2]);
		$p->{pspages} .= "@coords[0] @coords[1] newpath moveto\n";
		for (my $i = 0; $i < $total_elems; $i++) {
			my $angle = (@positions[$i]/$circle_size) * 360;
			my $radius = $INNER_RADIUS + (($OUTER_RADIUS-$INNER_RADIUS)*(@differences[($i*$num_graphs)+$j]/$max_diffs));
			my @new_coords = coords_on_circle($angle,$radius);
			$this_x = @new_coords[0];
			$this_y = @new_coords[1];
			$p->{pspages} .= "$this_x $this_y lineto\n";
			$last_x = $this_x;
			$last_y = $this_y;
		}
		$p->{pspages} .= "@coords[0] @coords[1] lineto\nclosepath\nstroke\n";
	}
}

# must return 1 for the file overall.
1;
