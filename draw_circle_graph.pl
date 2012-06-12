require "subfuncs.pl";
use PostScript::Simple;
use constant CENTER_X => 600; # X coordinate of circle center
use constant CENTER_Y => 600; # Y coordinate of circle center
use constant PS_X_SIZE => 1200; # X size of the PostScript object
use constant PS_Y_SIZE => 1200; # Y size of the PostScript object

my $usage  = "draw_circl_graph.pl input_file output_file\n";

my $in_file = shift or die $usage;
my $out_file = shift or die $usage;
my $p = new PostScript::Simple(	xsize => PS_X_SIZE,
								ysize => PS_Y_SIZE,
                                colour => 1,
                                eps => 0,
                                units => "bp");

$p->newpage;
draw_circle_graph_from_file($in_file, $p);
$p->output("$out_file.ps");

