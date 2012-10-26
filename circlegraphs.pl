use CircleGraph;
use Bio::AlignIO;
use Bio::SeqIO;

=pod

CircleGraph $circlegraph_obj draw_circle_graph ( String $datafile, CircleGraph $circlegraph_obj )

$datafile is a tab-delimited file:
    -   the first column contains positions around a circle,
        scaled to the size of the last row's value
    -   each subsequent column contains the values to graph, one graph per column
    -   the first row contains the names of the columns
$circlegraph_obj is an optional parameter for an existing CircleGraph.

The function returns a CircleGraph object with the values plotted,
positions mapped around the edge, and a legend describing the graph.

=cut

sub draw_circle_graph {
    my $datafile = shift;
    my $circlegraph_obj = shift;

    open DATAFH, "<$datafile" or die "$datafile failed to open\n";
    my $max_diffs = 0;
    my $total_elems;
    my @positions, @differences;

    my $line = readline DATAFH;
    $line =~ s/\n//;
    my @labels = split /\t/, $line;
    shift @labels;
    my @graphs = ();
    for (my $i=0; $i<@labels; $i++) {
        my @curr_array = ();
        push (@graphs, \@curr_array);
    }
    $line = readline DATAFH;

    while ($line ne "") {
        $line =~ s/\n//;
        my @items = split ('\t', $line);
        my $pos = shift @items;
        $total_elems = push (@positions, $pos);
        for (my $i=0; $i<@items; $i++) {
            @items[$i] =~ s/\n//;
            push (@{$graphs[$i]}, @items[$i]);
        }
        push (@differences, @items);
        $line = readline DATAFH;
    }

    my @sorted = sort (@differences);
    my $diff_len = @sorted;
    $max_diffs = @sorted[@sorted-1];
    $max_diffs =~ s/\n//;

    for(my $i=0; $i<@graphs; $i++) {
        for (my $j=0; $j<@{$graphs[$i]}; $j++) {
            @{$graphs[$i]}[$j] = (@{$graphs[$i]}[$j]/$max_diffs);
        }
    }

    unless ($circlegraph_obj) {
        $circlegraph_obj = new CircleGraph();
    }

    # draw background circles
    $circlegraph_obj->draw_circle($circlegraph_obj->inner_radius);
    $circlegraph_obj->draw_circle($circlegraph_obj->outer_radius);

    for (my $j=0; $j<@graphs; $j++) {
        $circlegraph_obj->set_color($j);
        $circlegraph_obj->plot_line(\@positions, $graphs[$j]);
        $circlegraph_obj->append_to_legend("@labels[$j]", "$j");
    }

    # draw labels around the edge
    $circlegraph_obj->set_font("Helvetica", 6, "black");
    my $circle_size = @positions[@positions-1];

    for (my $i = 0; $i < $total_elems; $i++) {
        my $angle = (@positions[$i]/$circle_size) * 360;
        my $radius = $circlegraph_obj->outer_radius + 10;
        $circlegraph_obj->circle_label($angle, $radius, "@positions[$i]");
    }

    $circlegraph_obj->append_to_legend("Maximum percent difference ($max_diffs) is scaled to 1");
    return $circlegraph_obj;
}

# must return 1 for the file overall.
1;
