#! /usr/bin/perl -w

package CircleGraph;

use strict;
use Carp;
use PostScript::Simple;

use constant PI 	=> 3.1415926535897932384626433832795;
use constant CENTER_X => 600; # X coordinate of circle center500
use constant CENTER_Y => 600; # Y coordinate of circle center
use constant PS_X_SIZE => 1200; # X size of the PostScript object1000
use constant PS_Y_SIZE => 1200; # Y size of the PostScript object

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

    my $legend_text = $self->legend;
    my $num_lines = split(/<item/,$legend_text) - 1;
    my $top_text = 10 + ($num_lines * 20);
    my $i = 0;

    while ($legend_text =~ s/<item color=(.+?)>(.+?)<\/item>//) {
        my $color = $1;
        my $text = $2;
        $self->set_color($color);
        $self->{ps_object}->text(10, $top_text - ($i*20), $text);
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
    my @colors_by_name = qw(red green blue yellow magenta cyan orange fuchsia dark_green violet);

    my %colors = (
        red => [255,0,0],
        green => [0,255,0],
        blue => [0,0,255],
        yellow => [200,200,0],
        magenta => [255,0,255],
        cyan => [0,255,255],
        orange => [255,100,0],
        fuchsia => [200,50,200],
        dark_green => [70,150,0],
        violet => [100,0,255],
        white => [255,255,255],
        black => [0,0,0]
    );

    if (ref($arg) =~ /ARRAY/) {
        print "array\n";
        @color = @$arg;
    } elsif (ref($arg) eq "") {
        if (exists $colors{$arg}) {
            @color = @{$colors{$arg}};
        } elsif ($arg =~ /\d+/) {
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

    if (ref($params)) {
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
    my @x_vals = @$arg1;
    my @y_vals = @$arg2;

    my $p = $self->{ps_object};

	my $window_size = @x_vals[1]-@x_vals[0];
	my $circle_size = @x_vals[(scalar @x_vals)-1] + $window_size;

    my @coords = $self->coords_on_circle(0,$INNER_RADIUS + (($OUTER_RADIUS-$INNER_RADIUS)*(@y_vals[0])));
    my ($last_x, $last_y, $this_x, $this_y);
    $last_x = @coords[0];
    $last_y = @coords[1];
    $p->{pspages} .= "@coords[0] @coords[1] newpath moveto\n";
    for (my $i = 0; $i < (scalar @x_vals); $i++) {
        my $angle = (@x_vals[$i]/$circle_size) * 360;
        my $radius = $INNER_RADIUS + (($OUTER_RADIUS-$INNER_RADIUS)*(@y_vals[$i]));
        my @new_coords = $self->coords_on_circle($angle,$radius);
        $this_x = @new_coords[0];
        $this_y = @new_coords[1];
        $p->{pspages} .= "$this_x $this_y lineto\n";
        $last_x = $this_x;
        $last_y = $this_y;
    }
    $p->{pspages} .= "closepath\nstroke\n";

}

sub circle_label {
    my $self = shift;
    my $angle = shift;
    my $radius = shift;
    my $text = shift;
    my $align = shift;
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

    my $p = $self->{ps_object};

    $p->{pspages} .= "newpath\n";
    $p->{pspages} .= CENTER_X . " ux " . CENTER_Y . " uy $radius $start_angle $stop_angle arc\n";
    $p->{pspages} .= CENTER_X . " ux " . CENTER_Y . " uy lineto\n";
    $p->{pspages} .= "closepath fill\n";

}

sub set_percent_red {
    my $self = shift;
	my $percent_red = shift;
	my $scaling = ($percent_red/100);

    my $p = $self->{ps_object};
	my @zero_red = (255,240,240);
	my @full_red = (204,0,0);
	my $r = int(((@full_red[0]-@zero_red[0])*$scaling) + @zero_red[0]);
	my $g = int(((@full_red[1]-@zero_red[1])*$scaling) + @zero_red[1]);
	my $b = int(((@full_red[2]-@zero_red[2])*$scaling) + @zero_red[2]);
	$p->setcolour($r, $g, $b);
}


1;
