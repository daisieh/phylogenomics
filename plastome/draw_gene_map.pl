#!/usr/bin/env perl
use FindBin;
use lib "$FindBin::Bin/../lib";
use CircleGraph;
use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use Getopt::Long;
use Pod::Usage;

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my $runline = "running " . basename($0) . " " . join (" ", @ARGV) . "\n";

my $out_file = 0;
my $datafile = 0;
my $plastid_name = "";
my $help = 0;
my $circle_size = 0;
my $ira = "";
my $irb = "";
my $color = "red";
my @special = ();

my $result = GetOptions ('fasta:s' => \$datafile,
            'outputfile=s' => \$out_file,
            'circlesize|size=i' => \$circle_size,
            'ira=s' => \$ira,
            'irb=s' => \$irb,
            'name=s' => \$plastid_name,
            'color=s' => \$color,
            'special=s{1,}' => \@special,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 2);
}
if (($out_file eq "") or ($datafile eq "") or ($circle_size == 0)) {
    pod2usage(-msg => "Need to specify output file ($out_file), input file ($datafile), and plastome size ($circle_size).", -exitval => 2);
}

print $runline;

my @inputs = ();
my @fwd_strand = ();
my @rev_strand = ();

open FH, "<:crlf", $datafile;
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
		push @inputs, "$label exons\t$start\t$stop";
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

my $circlegraph_obj = new CircleGraph();

my %special_genes = ();
foreach my $line (@special) {
	$line =~ /(.*)=(.*)/;
	chomp $2;
	if ($circlegraph_obj->get_color($2)) {
		$special_genes{$1} = $2;
	}
}

my $radius = $circlegraph_obj->inner_radius + 20;
my @fwd_strand_labels = ();

#draw forward strand along the outside
while (@fwd_strand > 0) {
	$line = shift @fwd_strand;
	my ($label, $start, $stop) = split /\t/, $line;
	chomp $stop;
	my $currcolor = $color;
	foreach my $code (keys %special_genes) {
		if ($label =~ /$code/) {
			$currcolor = $special_genes{$code};
			last;
		}
	}
	if ($label =~ /(.+) exons$/) {
		$label = $1;
		my $start_angle = ($start/$circle_size) * 360;
		my $stop_angle = ($stop/$circle_size) * 360;

		$circlegraph_obj->set_color_by_percent (30, {zero_color=>"white",full_color=>$currcolor});
		$circlegraph_obj->draw_filled_arc ($radius, $start_angle, $stop_angle, {width=>5});
		# label this element
		my $center_angle = ($start_angle + $stop_angle) / 2;
		push @fwd_strand_labels, "$label\t$center_angle";
		while ($line =~ /$label/) {
			$line = shift @fwd_strand;
			my (undef, $start, $stop) = split /\t/, $line;
			my $start_angle = ($start/$circle_size) * 360;
			my $stop_angle = ($stop/$circle_size) * 360;
			$circlegraph_obj->draw_filled_arc ($radius, $start_angle, $stop_angle, {color=>$currcolor, width=>5});
		}
		unshift @fwd_strand, $line;
	} else {
		my $start_angle = ($start/$circle_size) * 360;
		my $stop_angle = ($stop/$circle_size) * 360;
		$circlegraph_obj->draw_filled_arc ($radius, $start_angle, $stop_angle, {color=>$currcolor, width=>5});

		# label this element
		my $center_angle = ($start_angle + $stop_angle) / 2;
		push @fwd_strand_labels, "$label\t$center_angle";
	}
}

$circlegraph_obj->draw_circle($circlegraph_obj->inner_radius, {filled => 1, color => "white"});

my @rev_strand_labels = ();
#draw reverse strand along the inside
$radius = $circlegraph_obj->inner_radius;

while (@rev_strand > 0) {
	$line = shift @rev_strand;
	my ($label, $start, $stop) = split /\t/, $line;
	chomp $stop;

	my $currcolor = $color;
	foreach my $code (keys %special_genes) {
		if ($label =~ /$code/) {
			$currcolor = $special_genes{$code};
			last;
		}
	}
	if ($label =~ /exons (.+)$/) {
		$label = $1;
		my $start_angle = ($start/$circle_size) * 360;
		my $stop_angle = ($stop/$circle_size) * 360;

		$circlegraph_obj->set_color_by_percent (30, {zero_color=>"white",full_color=>$color});
		$circlegraph_obj->draw_filled_arc ($radius, $start_angle, $stop_angle, {width=>5});
		# label this element
		my $center_angle = ($start_angle + $stop_angle) / 2;
		push @rev_strand_labels, "$label\t$center_angle";
		while ($line =~ /$label/) {
			$line = shift @rev_strand;
			my (undef, $start, $stop) = split /\t/, $line;
			my $start_angle = ($start/$circle_size) * 360;
			my $stop_angle = ($stop/$circle_size) * 360;
			$circlegraph_obj->draw_filled_arc ($radius, $start_angle, $stop_angle, {color=>$color, width=>5});
		}
		unshift @rev_strand, $line;
	} else {
		my $start_angle = ($start/$circle_size) * 360;
		my $stop_angle = ($stop/$circle_size) * 360;
		$circlegraph_obj->draw_filled_arc ($radius, $start_angle, $stop_angle, {color=>$color, width=>5});

		# label this element
		my $center_angle = ($start_angle + $stop_angle) / 2;
		push @rev_strand_labels, "$label\t$center_angle";
	}
}

$circlegraph_obj->draw_circle($circlegraph_obj->inner_radius - 20, {filled => 1, color => "white"});

#label reverse strand along the inside
$circlegraph_obj->set_font("Helvetica", 16, "black");

for (my $i = 0; $i < @rev_strand_labels; $i++) {
	my $line = @rev_strand_labels[$i];
	my ($label, $center_angle) = split /\t/, $line;
	chomp $center_angle;

 	$circlegraph_obj->circle_label($center_angle, $circlegraph_obj->inner_radius - 22, $label, "right");

}

#label fwd strand along the outside
for (my $i = 0; $i < @fwd_strand_labels; $i++) {
	my $line = @fwd_strand_labels[$i];
	my ($label, $center_angle) = split /\t/, $line;
	chomp $center_angle;

 	$circlegraph_obj->circle_label($center_angle, $circlegraph_obj->inner_radius + 22, $label);

}

$circlegraph_obj->draw_circle($circlegraph_obj->inner_radius);

if ($ira ne "") {
	$ira =~ /(\d+)-(\d+)/;
	my $start_angle = ($1/$circle_size) * 360;
	my $stop_angle = ($2/$circle_size) * 360;
	$circlegraph_obj->draw_arc($circlegraph_obj->inner_radius, $start_angle, $stop_angle, {width=>5, color=>"black"});
}

if ($irb ne "") {
	$irb =~ /(\d+)-(\d+)/;
	my $start_angle = ($1/$circle_size) * 360;
	my $stop_angle = ($2/$circle_size) * 360;
	$circlegraph_obj->draw_arc($circlegraph_obj->inner_radius, $start_angle, $stop_angle, {width=>5, color=>"black"});
}

$circlegraph_obj->draw_center_text("$plastid_name\n$circle_size base pairs");
$circlegraph_obj->output_ps();
open OUT, ">", "$out_file.ps";
print OUT $circlegraph_obj->output_ps . "\n";
close OUT;


__END__

=head1 NAME

draw_gene_map

=head1 SYNOPSIS

draw_gene_map -fasta -output -circlesize -ira -irb -name

=head1 OPTIONS

  -fasta:       fasta file as output by DOGMA.
  -output:      output file
  -circlesize:  size of plastome
  -ira:         coordinates of inverted repeat A, written as xxxx-xxxx
  -irb:         coordinates of inverted repeat B, written as xxxx-xxxx
  -name:        name of plastome
  -color:       main color of gene map
  -special:     special colorization of genes, in the format ndh=green, psb=red, etc.

=head1 DESCRIPTION

Draws a plastome graph from a DOGMA-generated list of genes.

=cut

