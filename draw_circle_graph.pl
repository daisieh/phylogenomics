use PostScript::Simple;
use CircleGraph;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
require "subfuncs.pl";

print "running " . basename($0) . " " . join (" ", @ARGV) . "\n";
print "This should ony be for testing!!!\n";

my ($fastafile, $datafile, $out_file, $gb_file, $window_size) = 0;
GetOptions ('fasta:s' => \$fastafile,
            'data:s' => \$datafile,
            'outputfile:s' => \$out_file,
            'genbank|gb_file:s' => \$gb_file,
            'window:i' => \$window_size) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($fastafile) {
    unless ($window_size > 0) { pod2usage(-msg => "no window size specified.", -exitval => 2); }

    my $start_pos = 1;
    my $stop_pos = $window_size;
    my $len;
    my $flag = 1;

    my $curr_aln = make_aln_from_fasta_file ($fastafile);
    $datafile = "$out_file.txt";

    open FH, ">", "$datafile";
    truncate $diffs_file, 0;

    while ($flag >= 0) {
        my $gene_name = "$start_pos";
        $flag = perc_diff_partition ($curr_aln, $start_pos, $stop_pos);
        if ($flag >= 0) {
            my $val = $flag;
            print FH "$start_pos\t$val\n";
        }
        $start_pos += $window_size;
        $stop_pos += $window_size;
    }

    close FH;
}

open my $F, "<$datafile" or die "$datafile failed to open\n";
my $line = readline $F;
my $max_diffs = 0;
my $total_elems;
my @positions, @differences;

my @labels = split /\t/, $line;
my $num_graphs = @labels-1;
print "drawing $num_graphs graphs\n";
my %graphs = ();
for (my $i=0; $i<$num_graphs; $i++) {
    my @curr_array = ();
    $graphs{$i} = \@curr_array;
}

while ($line ne "") {
    my @items = split ('\t', $line);
    my $pos = shift @items;
    $total_elems = push (@positions, $pos);
    for (my $i=0; $i<@items; $i++) {
        @items[$i] =~ s/\n//;
        push (@{$graphs{$i}}, @items[$i]);
    }
    push (@differences, @items);
    $line = readline $F;
}

my @sorted = sort (@differences);
my $diff_len = @sorted;
$max_diffs = @sorted[@sorted-1];
$max_diffs =~ s/\n//;
print "max diffs is $max_diffs\n";

for(my $i=0; $i<scalar @differences; $i++) {
    @differences[$i] = (@differences[$i]/$max_diffs);
}

my $x = new CircleGraph();

# draw background circles
$x->draw_circle($x->inner_radius);
$x->draw_circle($x->outer_radius);

##todo: make this work for more than one column.
for (my $j=0; $j< $num_graphs; $j++) {
    $x->set_color($j);
    my @plot_me = @{$graphs{$j}};
    print "###" . ref ($graphs);
    $x->plot_line(\@positions, \@plot_me);
}
##end todo


# draw labels around the edge
$x->set_font("Helvetica", 6, "black");
my $circle_size = @positions[@positions-1];

for (my $i = 0; $i < $total_elems; $i++) {
    my $angle = (@positions[$i]/$circle_size) * 360;
    my $radius = $x->outer_radius + 10;
    $x->circle_label($angle, $radius, "@positions[$i]");
}

# draw legend in lower corner
$x->set_font("Helvetica", 12, "black");
$x->legend($x->legend . "Maximum percent difference ($max_diffs) is scaled to 1\n");
if ($window_size > 0) {
    $x->legend($x->legend . "Sliding window of $window_size bp\n");
}
$x->draw_legend_text;

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
