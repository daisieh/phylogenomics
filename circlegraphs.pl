use CircleGraph;
use Bio::AlignIO;
use Bio::SeqIO;

sub draw_circle_graph {
    my $datafile = shift;
    my $x = shift;

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

    unless ($x) {
        $x = new CircleGraph();
    }

    # draw background circles
    $x->draw_circle($x->inner_radius);
    $x->draw_circle($x->outer_radius);

    for (my $j=0; $j<@graphs; $j++) {
        $x->set_color($j);
        $x->plot_line(\@positions, $graphs[$j]);
        $x->append_to_legend("@labels[$j]", "$j");
    }

    # draw labels around the edge
    $x->set_font("Helvetica", 6, "black");
    my $circle_size = @positions[@positions-1];

    for (my $i = 0; $i < $total_elems; $i++) {
        my $angle = (@positions[$i]/$circle_size) * 360;
        my $radius = $x->outer_radius + 10;
        $x->circle_label($angle, $radius, "@positions[$i]");
    }

    $x->append_to_legend("Maximum percent difference ($max_diffs) is scaled to 1");
    return $x;
}

# must return 1 for the file overall.
1;
