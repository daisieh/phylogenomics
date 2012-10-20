# use PostScript::Simple;
use CircleGraph;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
require "subfuncs.pl";

print "running " . basename($0) . " " . join (" ", @ARGV) . "\n";

my ($fastafile, $datafile, $out_file, $gb_file, $window_size) = 0;
GetOptions ('fasta:s' => \$fastafile,
            'data:s' => \$datafile,
            'outputfile:s' => \$out_file,
            'genbank|gb_file:s' => \$gb_file,
            'window:i' => \$window_size) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

my $circle_size;

if ($gb_file) {
    open FH, ">", "$out_file.genes";
    print FH parse_genbank_file($gb_file);
    close FH;
}

if ($fastafile) {   # if we were given a fasta file, we should create the diffs file.
    unless ($window_size > 0) { pod2usage(-msg => "no window size specified.", -exitval => 2); }

    my $start_pos = 1;
    my $stop_pos = $window_size;
    my $len;
    my $flag = 1;

    my $curr_aln = make_aln_from_fasta_file ($fastafile);
    $datafile = "$out_file.diffs";

    open FH, ">", "$datafile" ;
    truncate $diffs_file, 0;
    while ($flag >= 0) {
        my $gene_name = "$start_pos";
        $flag = perc_diff_partition ($curr_aln, $start_pos, $stop_pos);
        if ($flag >= 0) {
            my $val = $flag;
            print FH "$start_pos\t$val\n";
        }
        $start_pos = $stop_pos;
        $stop_pos = $start_pos + $window_size;
    }
    $circle_size = $curr_aln->length();
    print FH "$circle_size\t-1\n";
    close FH;
}

open my $F, "<$datafile" or die "$datafile failed to open\n";
my $line = readline $F;
my $max_diffs = 0;
my @positions, @differences;
while ($line ne "") {
    my @items = split ('\t', $line);
    my $pos = shift @items;
    push (@positions, $pos);
    push (@differences, @items);
    $line = readline $F;
}
my $circle_size = pop @positions;
pop @differences;

my @sorted = sort (@differences);
my $diff_len = @sorted;
$max_diffs = @sorted[@sorted-1];
$max_diffs =~ s/\n//;

for(my $i=0; $i<@differences; $i++) {
    @differences[$i] = (@differences[$i]/$max_diffs);
}

my $x = new CircleGraph();
$x->draw_circle($x->inner_radius);
$x->draw_circle($x->outer_radius);
$x->set_color(0);
$x->plot_line(\@positions, \@differences);

# labeling positions along the outer circle
$x->set_font("Helvetica", 6, "black");
my $radius = $x->outer_radius;
$x->circle_label(0, $radius, "  1");

for (my $i = 0; $i < @positions; $i++) {
    my $angle = (@positions[$i]/$circle_size) * 360;
    my $label = "  @positions[$i]";
    if ((($label % 1000) != 0) && ($label != 0)) {
        $label = "-";
    }
    $x->circle_label($angle, $radius, "$label");
}

$x->set_font("Helvetica", 12, "black");
$x->legend($x->legend . "Maximum percent difference ($max_diffs) is scaled to 1\n");
if ($window_size > 0) {
    $x->legend($x->legend . "Sliding window of $window_size bp\n");
}
$x->draw_legend_text;

if ($gb_file) {
    open INPUTFILE, "<$out_file.genes" or die "$out_file.genes failed to open\n";
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

    my @labels = ();
    for (my $i = 0; $i < @inputs; $i++) {
        my $line = @inputs[$i];
        my ($label, $start, $stop, $value) = split /\t/, $line;
        $value =~ s/\n//;
        if ($value eq "") {
            $value = 0;
        }

        my $start_angle = ($start/$circle_size) * 360;
        my $stop_angle = ($stop/$circle_size) * 360;
        my $radius = $x->inner_radius;

        $x->set_percent_red((1-$value)*100);
        $x->draw_filled_arc ($radius, $start_angle, $stop_angle);

        # label this element
        my $center_angle = ($start_angle + $stop_angle) / 2;
        push @labels, "$label\t$center_angle";
    }

    $x->draw_circle($x->inner_radius - 5, {filled => 1, color => "white"});
    $x->draw_circle($x->inner_radius);
    $x->set_font("Helvetica", 6, "black");
    foreach my $line (@labels) {
        $line =~ /(.+?)\t(.+?)$/;
        $x->circle_label($2, $x->inner_radius - 5, $1, "right");
    }
}

$x->output_ps();
open OUT, ">", "$out_file.ps" or die "couldn't make output file $out_file";
print OUT $x->output_ps . "\n";
close OUT;

__END__

=head1 NAME

tree_omega

=head1 SYNOPSIS

slice_fasta_file [options]

=head1 OPTIONS

    -fasta:     fasta file of aligned sequences
	-outputfile:    name of output file or directory for output files (if using -genbank)
	-genbank|gb_file:	genbank file specifying genes
	-start:	start position of single slice
	-end:   end position of single slice

=head1 DESCRIPTION

=cut
