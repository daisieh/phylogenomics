require "subfuncs.pl";

use CircleGraph;
use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use Getopt::Long;
use Pod::Usage;
require "circlegraphs.pl";

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my $runline = "running " . basename($0) . " " . join (" ", @ARGV) . "\n";

my ($gb_file, $out_file) = 0;
my $datafile = 0;
my $plastid_name = "";
my $help = 0;
my $circle_size = 0;
my $ira = "";
my $irb = "";


my $result = GetOptions ('fasta:s' => \$datafile,
            'outputfile=s' => \$out_file,
            'circlesize|size=i' => \$circle_size,
            'ira=s' => \$ira,
            'irb=s' => \$irb,
            'name:s' => \$plastid_name,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

print $runline;

my @inputs = ();
my @fwd_strand = ();
my @rev_strand = ();

open FH, "<", $datafile;
my $line = readline FH;

while ($line) {
	$line =~ m/>(.*?): (.*)$/;
	my $label = $1;
	my $location = $2;
	if ($location =~ /,/) { #we have exons.
		my $start, $stop = 0;
		my @exon_block = ();
		while ($location =~ /.*?(\d+)-(\d+)(.*)/) {
			if ($start == 0) {
				$start = $1;
			}
			push @exon_block, "$label\t$1\t$2";
			$location = $3;
			$stop = $2;
		}
		push @inputs, "exons $label\t$start\t$stop";
		push @inputs, @exon_block;

	} else {
		$location =~ /.*?(\d+)-(\d+).*/;
		$location = "$1\t$2";
		push @inputs, "$label\t$location";
	}
	$line = readline FH;
	while ($line !~ m/>(.*)$/) {
		chomp $line;
		$line = readline FH;
		unless ($line) { last; }
	}
}
close FH;

foreach my $line (@inputs) {
	my ($label, $start, $stop) = split /\t/, $line;
	chomp $stop;
	if ($label =~ /(.*) \(reverse strand\)/) {
		$label = $1;
		push @rev_strand,  "$label\t$start\t$stop\n";
	} else {
		push @fwd_strand, "$label\t$start\t$stop\n";
	}
}

my $x = new CircleGraph();
my $radius = $x->inner_radius + 20;
my @fwd_strand_labels = ();

#draw forward strand along the outside
while (@fwd_strand > 0) {
	$line = shift @fwd_strand;
	my ($label, $start, $stop) = split /\t/, $line;
	chomp $stop;

	if ($label =~ /exons (.+)$/) {
		$label = $1;
		my $start_angle = ($start/$circle_size) * 360;
		my $stop_angle = ($stop/$circle_size) * 360;

		$x->set_color_by_percent (30);
		$x->draw_filled_arc ($radius, $start_angle, $stop_angle, {width=>5});
		# label this element
		my $center_angle = ($start_angle + $stop_angle) / 2;
		push @fwd_strand_labels, "$label\t$center_angle";
		while ($line =~ /$label/) {
			$line = shift @fwd_strand;
			my (undef, $start, $stop) = split /\t/, $line;
			my $start_angle = ($start/$circle_size) * 360;
			my $stop_angle = ($stop/$circle_size) * 360;
			$x->draw_filled_arc ($radius, $start_angle, $stop_angle, {color=>"red", width=>5});
		}
		unshift @fwd_strand, $line;
	} else {
		my $start_angle = ($start/$circle_size) * 360;
		my $stop_angle = ($stop/$circle_size) * 360;
		$x->draw_filled_arc ($radius, $start_angle, $stop_angle, {color=>"red", width=>5});

		# label this element
		my $center_angle = ($start_angle + $stop_angle) / 2;
		push @fwd_strand_labels, "$label\t$center_angle";
	}
}

#label fwd strand along the inside
$x->set_font("Helvetica", 10, "black");
for (my $i = 0; $i < @fwd_strand_labels; $i++) {
	my $line = @fwd_strand_labels[$i];
	my ($label, $center_angle) = split /\t/, $line;
	chomp $center_angle;

 	$x->circle_label($center_angle, $x->inner_radius + 22, $label);

}


$x->draw_circle($x->inner_radius, {filled => 1, color => "white"});

my @rev_strand_labels = ();
#draw reverse strand along the inside
$radius = $x->inner_radius;

while (@rev_strand > 0) {
	$line = shift @rev_strand;
	my ($label, $start, $stop) = split /\t/, $line;
	chomp $stop;

	if ($label =~ /exons (.+)$/) {
		$label = $1;
		my $start_angle = ($start/$circle_size) * 360;
		my $stop_angle = ($stop/$circle_size) * 360;

		$x->set_color_by_percent (30);
		$x->draw_filled_arc ($radius, $start_angle, $stop_angle, {width=>5});
		# label this element
		my $center_angle = ($start_angle + $stop_angle) / 2;
		push @rev_strand_labels, "$label\t$center_angle";
		while ($line =~ /$label/) {
			$line = shift @rev_strand;
			my (undef, $start, $stop) = split /\t/, $line;
			my $start_angle = ($start/$circle_size) * 360;
			my $stop_angle = ($stop/$circle_size) * 360;
			$x->draw_filled_arc ($radius, $start_angle, $stop_angle, {color=>"red", width=>5});
		}
		unshift @rev_strand, $line;
	} else {
		my $start_angle = ($start/$circle_size) * 360;
		my $stop_angle = ($stop/$circle_size) * 360;
		$x->draw_filled_arc ($radius, $start_angle, $stop_angle, {color=>"red", width=>5});

		# label this element
		my $center_angle = ($start_angle + $stop_angle) / 2;
		push @rev_strand_labels, "$label\t$center_angle";
	}
}

$x->draw_circle($x->inner_radius - 20, {filled => 1, color => "white"});

#label reverse strand along the inside
$x->set_font("Helvetica", 10, "black");
for (my $i = 0; $i < @rev_strand_labels; $i++) {
	my $line = @rev_strand_labels[$i];
	my ($label, $center_angle) = split /\t/, $line;
	chomp $center_angle;

 	$x->circle_label($center_angle, $x->inner_radius - 22, $label, "right");

}

$x->draw_circle($x->inner_radius);

$ira =~ /(\d+)-(\d+)/;
my $start_angle = ($1/$circle_size) * 360;
my $stop_angle = ($2/$circle_size) * 360;
$x->draw_arc($x->inner_radius, $start_angle, $stop_angle, {width=>5, color=>"black"});
$irb =~ /(\d+)-(\d+)/;
my $start_angle = ($1/$circle_size) * 360;
my $stop_angle = ($2/$circle_size) * 360;
$x->draw_arc($x->inner_radius, $start_angle, $stop_angle, {width=>5, color=>"black"});

$x->draw_center_text("$plastid_name\n$circle_size base pairs");

$x->output_ps();
open OUT, ">", "$out_file.ps";
print OUT $x->output_ps . "\n";
close OUT;


__END__

=head1 NAME

draw_gene_map

=head1 SYNOPSIS

draw_gene_map [-locations locationfile |-gb genbankfile] -output

=head1 OPTIONS

  -locations|genbank:   either a genbank file with genes to map or a tab-delimited file of genes and locations.
  -output:              output file

=head1 DESCRIPTION

draws a plastome graph from a genbank file or from a tab-delimited list of genes and locations.

=cut

