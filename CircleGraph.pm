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

my $MAIN_RADIUS				= 387.5;
my $OUTER_RADIUS = 387.5;
my $INNER_RADIUS = 337.5;

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


sub new {
    my ($class, %data) = @_;
    my $self = {
        ps_object => undef,
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

sub draw_inner_circle {
    my $self = shift;
    my $p = $self->{ps_object};
    $p->setcolour(0,0,0);
	$p->setlinewidth(1);
	$p->circle(CENTER_X,CENTER_Y, $INNER_RADIUS);
}

sub draw_outer_circle {
    my $self = shift;
    my $p = $self->{ps_object};
    $p->setcolour(0,0,0);
	$p->setlinewidth(1);
	$p->circle(CENTER_X,CENTER_Y, $OUTER_RADIUS);
}

sub plot_line {
    my $self = shift;
    my $arg1 = shift;
    my $arg2 = shift;
    my $arg3 = shift;
    my @x_vals = @$arg1;
    my @y_vals = @$arg2;
    my @color;
    if (ref($arg3) =~ /ARRAY/) {
        @color = @$arg3;
    } elsif (ref($arg3) eq "") {
        @color = @{@colours[$arg3]};
        print "@color[0]\n";
    } else {
        @color = (0,0,0);
    }
    my $p = $self->{ps_object};

	my $window_size = @x_vals[1]-@x_vals[0];
	my $circle_size = @x_vals[(scalar @x_vals)-1] + $window_size;

    my @coords = coords_on_circle(0,$INNER_RADIUS + (($OUTER_RADIUS-$INNER_RADIUS)*(@y_vals[0])));
    my ($last_x, $last_y, $this_x, $this_y);
    $last_x = @coords[0];
    $last_y = @coords[1];
    $p->{pspages} .= "@color[0] @color[1] @color[2] setrgbcolor\n";
    $p->{pspages} .= "@coords[0] @coords[1] newpath moveto\n";
    for (my $i = 0; $i < (scalar @x_vals); $i++) {
        my $angle = (@x_vals[$i]/$circle_size) * 360;
        my $radius = $INNER_RADIUS + (($OUTER_RADIUS-$INNER_RADIUS)*(@y_vals[$i]));
        my @new_coords = coords_on_circle($angle,$radius);
        $this_x = @new_coords[0];
        $this_y = @new_coords[1];
        $p->{pspages} .= "$this_x $this_y lineto\n";
        $last_x = $this_x;
        $last_y = $this_y;
    }
    $p->{pspages} .= "closepath\nstroke\n";

}

sub coords_on_circle {
	my $angle = shift;
	my $radius = shift;

	return ( CENTER_X + ($radius*cos(($angle * PI)/180)), CENTER_Y + ($radius*sin(($angle * PI)/180)));
}


1;
