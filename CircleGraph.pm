#! /usr/bin/perl -w

package CircleGraph;

use strict;
use Carp;
use PostScript::Simple;

use constant PI 	=> 3.1415926535897932384626433832795;
use constant CENTER_X => 500; # X coordinate of circle center
use constant CENTER_Y => 500; # Y coordinate of circle center
use constant PS_X_SIZE => 1000; # X size of the PostScript object
use constant PS_Y_SIZE => 1000; # Y size of the PostScript object

#
# graphics default values
#

my $MAIN_RADIUS	 = 387.5;
my $OUTER_RADIUS = 387.5;
my $INNER_RADIUS = 337.5;


sub new {
    my ($class, %data) = @_;
    my $self = {
        ps_object => undef,
        legend => "",
    };

    foreach (keys %data) {
        $self->{$_} = $data{$_};
    }

    bless $self, $class;
    $self->init();

    return $self;
}

sub init {
    my $self = shift;
    $self->{ps_object} = new PostScript::Simple( colour => 1, eps => 0, units => "bp", xsize => PS_X_SIZE, ysize => PS_Y_SIZE );
    $self->{ps_object}->newpage();
}

sub ps_string {
    my $self = shift;
    return $self->{ps_object}->{pspages};
}

sub ps_object {
    my $self = shift;
    return $self->{ps_object};
}

sub inner_radius {
    my $self = shift;
    my $arg = shift;
    if ($arg) {
        $INNER_RADIUS = $arg;
    }
    return $INNER_RADIUS;
}

sub main_radius {
    my $self = shift;
    my $arg = shift;
    if ($arg) {
        $MAIN_RADIUS = $arg;
    }
    return $MAIN_RADIUS;
}

sub outer_radius {
    my $self = shift;
    my $arg = shift;
    if ($arg) {
        $OUTER_RADIUS = $arg;
    }
    return $OUTER_RADIUS;
}

sub legend {
    my $self = shift;
    my $arg = shift;
    if ($arg) {
        $self->{legend} = "<item color=black>$arg</item>";
    }
    return $self->{legend};
}

sub append_to_legend {
    my $self = shift;
    my $text = shift;
    my $color = shift;
    if (defined($color)) {
        $self->{legend} .= "<item color=$color>$text</item>";
    } else {
        $self->{legend} .= "<item color=black>$text</item>";
    }
}

sub draw_legend_text {
    my $self = shift;
    my $params = shift;

    my $font = "Helvetica";
    my $size = 12;
    my $height = 20;
    my $margin = 10;

    if (ref($params) eq "HASH") {
        if ($params->{"font"}) {
            $font = $params->{"font"};
        }
        if ($params->{"size"}) {
            $size = $params->{"size"};
        }
        if ($params->{"height"}) {
            $height = $params->{"height"};
        }
        if ($params->{"margin"}) {
            $margin = $params->{"margin"};
        }
    }

    $self->set_font($font, $size);

    my $legend_text = $self->legend;
    my $num_lines = split(/<item/,$legend_text) - 1;
    my $top_text = $margin + ($num_lines * $height);
    my $i = 0;

    while ($legend_text =~ s/<item color=(.+?)>(.+?)<\/item>//) {
        my $color = $1;
        my $text = $2;
        $self->set_color($color);
        $self->{ps_object}->text($margin, $top_text - ($i*$height), $text);
        $i++;
    }
}

sub output_ps {
    my $self = shift;
    my $output_string = "";
    my $ps_array = $self->{ps_object}->_builddocument;
    foreach my $i (@$ps_array) {
        if (ref($i) eq "SCALAR") {
            $output_string .= $$i;
        } else {
            $output_string .= $i;
        }
    }
    return $output_string;
}

sub set_color {
    my $self = shift;
    my $arg = shift;
    my @color = (0,0,0);
    my @colors_by_name = qw(orange dark_green blue yellow red green magenta cyan fuchsia violet grey brown slate pink gold tardis teal brick);

    my %colors = (
        blue => [0.04,0.33,0.63],
        orange => [1,0.57,0.12],
        dark_green => [0,0.5,0.08],
        red => [1,0,0],
        green => [0,1,0],
        yellow => [0.8,0.8,0],
        magenta => [1,0,1],
        cyan => [0,1,1],
        fuchsia => [0.8,0.2,0.8],
        violet => [0.4,0,1],
        white => [1,1,1],
        black => [0,0,0],
        grey => [0.5,0.5,0.5],
        light_grey => [0.8,0.8,0.8],
        dark_grey => [0.3,0.3,0.3],
        brown => [0.65,0.35,0.2],
        slate => [0.35,0.45,0.6],
        pink => [1,0.8,0.8],
        gold => [1,0.8,0],
        tardis => [0,0,0.65],
        teal => [0.2,0.75,0.75],
        brick => [0.6,0.05,0.05],
    );

    if (ref($arg) =~ /ARRAY/) {
        @color = @$arg;
        if ((@color[0]>1) && (@color[1]>1) && (@color[2]>1)) {
            @color[0] = @color[0]/255;
            @color[1] = @color[1]/255;
            @color[2] = @color[2]/255;
        }
    } elsif (ref($arg) eq "") {
        if (exists $colors{$arg}) {
            @color = @{$colors{$arg}};
        } elsif ($arg =~ /\d+/) {
            $arg = $arg % @colors_by_name;
            @color = @{$colors{@colors_by_name[$arg]}};
        }
    }
    $self->{ps_object}->{pspages} .= "@color[0] @color[1] @color[2] setrgbcolor\n";
}

sub set_font {
    my $self = shift;
    my $font = shift;
    my $size = shift;
    my $color = shift;

    $self->ps_object->setfont($font, $size);
    if ($color) { $self->set_color($color); }
}

sub draw_circle {
    my $self = shift;
    my $r = shift;
    my $params = shift;
    my $filled = 0;
    my $linewidth = 1;
    my $p = $self->{ps_object};

    if (ref($params) eq "HASH") {
        if ($params->{filled}) {
            $filled = $params->{filled};
        }
        if ($params->{color}) {
            $self->set_color($params->{color});
        }
        if ($params->{width}) {
            $linewidth = $params->{width};
        }
    } else {
        $self->set_color("black");
    }
    $p->{pspages} .= "$linewidth u setlinewidth\n";
    $p->{pspages} .= "newpath\n";
    $p->{pspages} .= CENTER_X . " ux " . CENTER_Y . " uy $r 0 360 arc closepath\n";
    if ($filled) {
        $p->{pspages} .= "fill\n";
    } else {
        $p->{pspages} .= "stroke\n";
    }
}

sub plot_line {
    my $self = shift;
    my $arg1 = shift;
    my $arg2 = shift;
	my $params = shift;
    my @x_vals = @$arg1;
    my @y_vals = @$arg2;

    my $p = $self->{ps_object};
	my $color;
	my $width = 1;

    if (ref($params) eq "HASH") {
        if (exists $params->{"color"}) {
            $color = $params->{"color"};
        }
        if (exists $params->{"width"}) {
            $width = $params->{"width"};
        }
    }


	my $window_size = @x_vals[1]-@x_vals[0];
	my $circle_size = @x_vals[(scalar @x_vals)-1] + $window_size;

    my @coords = $self->coords_on_circle(0,$INNER_RADIUS + (($OUTER_RADIUS-$INNER_RADIUS)*(@y_vals[0])));
    my ($last_x, $last_y, $last_angle, $this_x, $this_y, $this_angle, $this_radius);
    my @blank_sectors;
    $last_x = @coords[0];
    $last_y = @coords[1];
    $last_angle = 0;

    if (defined ($color)) {
    	$self->set_color("$color");
    }
    $p->{pspages} .= "$width u setlinewidth\n";
    $p->{pspages} .= "@coords[0] @coords[1] newpath moveto\n";
    for (my $i = 0; $i < (scalar @x_vals); $i++) {
        $this_angle = (@x_vals[$i]/$circle_size) * 360;
        $this_radius = $INNER_RADIUS + (($OUTER_RADIUS-$INNER_RADIUS)*(@y_vals[$i]));

        # if the value to map is negative, we don't actually want to map it as-is:
        # map radius 0 instead, then cover it up with a sector.
        if (@y_vals[$i] < 0) {
            $this_radius = $self->inner_radius;
            my @angle_to_blank = ($last_angle, $this_angle);
            push @blank_sectors, \@angle_to_blank;
        }

        my @new_coords = $self->coords_on_circle($this_angle,$this_radius);
        $this_x = @new_coords[0];
        $this_y = @new_coords[1];
        $p->{pspages} .= "$this_x $this_y lineto\n";


        $last_x = $this_x;
        $last_y = $this_y;
        $last_angle = $this_angle;
    }
    $p->{pspages} .= "closepath\nstroke\n";

    if (@blank_sectors > 0) {
        foreach my $sector (@blank_sectors) {
            my $start_angle = @$sector[0];
            my $end_angle = @$sector[1];

            $self->draw_filled_arc ($self->outer_radius, $start_angle, $end_angle, $color);
        }
    }

}

sub plot_points {
    my $self = shift;
    my $arg1 = shift;
    my $arg2 = shift;
	my $params = shift;
    my @x_vals = @$arg1;
    my @y_vals = @$arg2;

    my $p = $self->{ps_object};
	my $color;
	my $width = 1;
	my $radius = $INNER_RADIUS;
	my $angle_size = 0.1;

    if (ref($params) eq "HASH") {
        if (exists $params->{"color"}) {
            $color = $params->{"color"};
        }
        if (exists $params->{"width"}) {
            $width = $params->{"width"};
        }
        if (exists $params->{"radius"}) {
            $radius += $params->{"radius"};
        }
        if (exists $params->{"angle"}) {
            $angle_size = $params->{"angle"};
        }
    }


	my $window_size = @x_vals[1]-@x_vals[0];
	my $circle_size = @x_vals[(scalar @x_vals)-1] + $window_size;

    my @coords = $self->coords_on_circle(0,$radius);
    my ($last_x, $last_y, $last_angle, $this_x, $this_y, $this_angle);
    my @blank_sectors;
    $last_x = @coords[0];
    $last_y = @coords[1];
    $last_angle = 0;

    if (defined ($color)) {
    	$self->set_color("$color");
    }
    $p->{pspages} .= "$width u setlinewidth\n";
    for (my $i = 0; $i < (scalar @x_vals); $i++) {
    	if (@y_vals[$i]) {
			$this_angle = (@x_vals[$i]/$circle_size) * 360;
			my $increment_angle = $this_angle + $angle_size;
			my @new_coords = $self->coords_on_circle($this_angle,$radius);
			$this_x = @new_coords[0];
			$this_y = @new_coords[1];
			$p->{pspages} .= "$this_x $this_y newpath moveto\n";
			@new_coords = $self->coords_on_circle($increment_angle,$radius);
			$this_x = @new_coords[0];
			$this_y = @new_coords[1];
			$p->{pspages} .= "$this_x $this_y lineto\n";
			$p->{pspages} .= "closepath\nstroke\n";
		}
    }
}


sub circle_label {
    my $self = shift;
    my $angle = shift;
    my $radius = shift;
    my $text = shift;
    my $align = shift;

    if ($text eq "") {
    	return;
    }

    my @new_coords = $self->coords_on_circle($angle,$radius);
    if ($align) {
        $self->ps_object->text( {rotate => $angle, align => $align}, @new_coords[0], @new_coords[1], "$text");
    } else {
        $self->ps_object->text( {rotate => $angle}, @new_coords[0], @new_coords[1], "$text");
    }
}

sub coords_on_circle {
    my $self = shift;
	my $angle = shift;
	my $radius = shift;

	return ( CENTER_X + ($radius*cos(($angle * PI)/180)), CENTER_Y + ($radius*sin(($angle * PI)/180)));
}

sub draw_filled_arc {
    my $self = shift;
	my $radius = shift;
	my $start_angle = shift;
	my $stop_angle = shift;
	my $params = shift;

	my $color = "black";

    if (ref($params) eq "HASH") {
        if (exists $params->{"color"}) {
            $color = $params->{"color"};
        }
    } else { # it's a deprecated color specification
    	$color = $params;
    }

    $self->set_color("$color");
    my $p = $self->{ps_object};

    $p->{pspages} .= "newpath\n";
    $p->{pspages} .= CENTER_X . " ux " . CENTER_Y . " uy $radius $start_angle $stop_angle arc\n";
    $p->{pspages} .= CENTER_X . " ux " . CENTER_Y . " uy lineto\n";
    $p->{pspages} .= "closepath fill\n";

}

sub draw_arc {
    my $self = shift;
	my $radius = shift;
	my $start_angle = shift;
	my $stop_angle = shift;
	my $params = shift;

	my $color = "black";
	my $width = 5;

    if (ref($params) eq "HASH") {
        if (exists $params->{"color"}) {
            $color = $params->{"color"};
        }
        if (exists $params->{"width"}) {
            $width = $params->{"width"};
        }
    }

    $self->set_color("$color");
    my $p = $self->{ps_object};

    $p->{pspages} .= "newpath\n";
    $p->{pspages} .= "$width u setlinewidth\n";
    $p->{pspages} .= CENTER_X . " ux " . CENTER_Y . " uy $radius $start_angle $stop_angle arc\n";
    $p->{pspages} .= "stroke\n";

}

sub set_color_by_percent {
    my $self = shift;
	my $percent_red = shift;
	my $zero_color = shift;
	my $full_color = shift;
	my $scaling = ($percent_red/100);

    # red is the default color
	my @zero_red = (255,240,240);
	my @full_red = (204,0,0);

    if (ref($zero_color) eq "ARRAY") {
        @zero_red = @$zero_color;
    }
    if (ref($full_color) eq "ARRAY") {
        @full_red = @$full_color;
    }

	my $r = int(((@full_red[0]-@zero_red[0])*$scaling) + @zero_red[0]);
	my $g = int(((@full_red[1]-@zero_red[1])*$scaling) + @zero_red[1]);
	my $b = int(((@full_red[2]-@zero_red[2])*$scaling) + @zero_red[2]);
	$self->{ps_object}->setcolour($r, $g, $b);
}


1;
