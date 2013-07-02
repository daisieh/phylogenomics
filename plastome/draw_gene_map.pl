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
my $min_coverage = 0;
my $circle_size = 0;


GetOptions ('locations:s' => \$datafile,
            'outputfile=s' => \$out_file,
            'genbank|gb:s' => \$gb_file,
            'name' => \$plastid_name,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

print $runline;
unless ($datafile) {
	# process the gb file
	if ($gb_file) {
		if ($gb_file =~ /\.gb$/) {
			my ($fh, $filename) = tempfile(UNLINK => 1);
			my $gb_locs = get_locations_from_genbank_file($gb_file);
			print $fh $gb_locs;
			close $fh;
			$datafile = "$filename";

			$gb_locs =~ /(.+)\n(.+?)$/;
			my $source = $2;
			$source =~ /.+?\t.+?\t(.+)$/;
			$circle_size = $1;
		}
	} else {
		pod2usage (-msg => "Need to specify either a locations file or a genbank file.", -exitval => 2);
	}
}

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

my @fwd_strand = ();
my @rev_strand = ();
for (my $i = 0; $i < @inputs; $i++) {
	my $line = @inputs[$i];
	my ($label, $start, $stop) = split /\t/, $line;
	chomp $stop;
	if ($start > $stop) {
		push @rev_strand,  "$label\t$stop\t$start\n";
	} else {
		push @fwd_strand, "$label\t$start\t$stop\n";
	}
}

my $x = new CircleGraph();

#draw forward strand along the outside
for (my $i = 0; $i < @fwd_strand; $i++) {
	my $line = @fwd_strand[$i];
	my ($label, $start, $stop) = split /\t/, $line;
	chomp $stop;

	my $start_angle = ($start/$circle_size) * 360;
	my $stop_angle = ($stop/$circle_size) * 360;
 	my $radius = $x->inner_radius + 20;

	$x->draw_filled_arc ($radius, $start_angle, $stop_angle, {color=>"red", width=>5});
	print "$label\t$start_angle\t$stop_angle\n";

	# label this element
	my $center_angle = ($start_angle + $stop_angle) / 2;
    $x->set_font("Helvetica", 10, "black");
 	$x->circle_label($center_angle, $x->inner_radius + 22, $label);

}

$x->draw_circle($x->inner_radius, {filled => 1, color => "white"});
# $x->draw_circle($x->inner_radius);

#draw reverse strand along the inside
for (my $i = 0; $i < @rev_strand; $i++) {
	my $line = @rev_strand[$i];
	my ($label, $start, $stop) = split /\t/, $line;
	chomp $stop;

	my $start_angle = ($start/$circle_size) * 360;
	my $stop_angle = ($stop/$circle_size) * 360;
 	my $radius = $x->inner_radius;

	$x->draw_filled_arc ($radius, $start_angle, $stop_angle, {color=>"red", width=>5});
	print "$label\t$start_angle\t$stop_angle\n";

}
$x->draw_circle($x->inner_radius - 20, {filled => 1, color => "white"});

#label reverse strand along the inside
for (my $i = 0; $i < @rev_strand; $i++) {
	my $line = @rev_strand[$i];
	my ($label, $start, $stop) = split /\t/, $line;
	chomp $stop;

	my $start_angle = ($start/$circle_size) * 360;
	my $stop_angle = ($stop/$circle_size) * 360;

	# label this element
	my $center_angle = ($start_angle + $stop_angle) / 2;
    $x->set_font("Helvetica", 10, "black");
 	$x->circle_label($center_angle, $x->inner_radius - 22, $label, "right");

}

$x->draw_circle($x->inner_radius);


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

