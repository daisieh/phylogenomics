use PostScript::Simple;
use Bio::AlignIO;
use Bio::SeqIO;

use constant PI 	=> 3.1415926535897932384626433832795;
use constant CENTER_X => 600; # X coordinate of circle center500
use constant CENTER_Y => 600; # Y coordinate of circle center
use constant PS_X_SIZE => 1200; # X size of the PostScript object1000
use constant PS_Y_SIZE => 1200; # Y size of the PostScript object

#
# graphics default values
#
my $MAIN_RADIUS				= 387.5;
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

sub convert_aln_to_nexus {
	my $aln = shift;
	my $blocksize=2000;
	my $nexblock = "";
	my $result = "";
	my $ntax = 0;
	my $nchar = 1;
	my $i = 1;
	my $flag = 1;
	my $len;
	while ($flag) {
		$ntax=0;
		foreach my $seq ( $aln->each_seq()) {
			$len = $seq->length;
			if ($i > $len - $blocksize) { $flag = 0; $blocksize = $len - $i + 1;}
			$nexblock .= "" . $seq->display_name . "\t";
			$nchar = length($seq->seq());
			my $seq_str = $seq->subseq($i, $i+$blocksize-1);
			$nexblock .= "$seq_str\n";
			$ntax++;
		}
		$nexblock .= "\n";
		$i += $blocksize;
	}
	$result .= "#NEXUS\n\nBegin DATA;\nDimensions ntax=$ntax nchar=$nchar;\n";
	$result .= "Format datatype=dna gap=- interleave=yes;\n";
 	$result .= "Matrix\n$nexblock\n;\nEnd;\n";

	return $result;
}

sub perc_diff_partition {
	my $newaln = shift;
	my $start_pos = shift;
	my $stop_pos = shift;

	if ($stop_pos < $newaln->length()) {
		my $aln_slice = $newaln->slice($start_pos, $stop_pos);
		my $p = $aln_slice->percentage_identity();
		return 100-$p;
	} else {
		return -1;
	}
}

sub make_aln_from_fasta_file {
	my $fa_file = shift;
	my $min_length = 0;

	my $inseq = Bio::SeqIO->new(-file => "<$fa_file", -format => "fasta");
	my $newaln = Bio::SimpleAlign->new();

	my $seq;
	eval {$seq = $inseq->next_seq;} or die "file not in fasta format.\n";
	$min_length = $seq->length();
	while ($seq ne "") {
		my $name = $seq->display_name;
		#remove weird chars, file suffixes
		$name =~ s/(.*?)\..*/$1/;
		#$name =~ s/[\Q !@#$%^&*.-?<>,|\/\E]//g;
		#shorten name if it's too long
		if (length($name) > 12) {
			$name =~ /(.{12})/;
			$name = $1;
		}
 		my $outseq = Bio::LocatableSeq->new(-seq => $seq->seq(), -id => $name);
 		$newaln->add_seq ($outseq);
 		if ($min_length > $seq->length() ) {
 			$min_length = $seq->length();
 		}
 		$seq = $inseq->next_seq;
	}

	my $flush_aln = $newaln->slice(1,$min_length);

	return $flush_aln;
}

sub test_ps { # create a new PostScript object
	my $outfile = shift;
	my $p = new PostScript::Simple( colour => 1, eps => 0, units => "bp", xsize => PS_X_SIZE, ysize => PS_Y_SIZE );
	$p->newpage;

	# draw some lines and other shapes
	$p->circle(CENTER_X,CENTER_Y, $MAIN_RADIUS);
	$p->circle(CENTER_X,CENTER_Y, $MAIN_RADIUS-50);

	my @inner_coords = coords_on_circle(45,$MAIN_RADIUS-50);
	my @outer_coords = coords_on_circle(45,$MAIN_RADIUS);
	$p->line ( @inner_coords[0], @inner_coords[1], @outer_coords[0], @outer_coords[1]);

	@inner_coords = coords_on_circle(22.5,$MAIN_RADIUS-50);
	$p->line ( @inner_coords[0], @inner_coords[1], @outer_coords[0], @outer_coords[1]);

	@outer_coords = coords_on_circle(0,$MAIN_RADIUS);
	$p->line ( @inner_coords[0], @inner_coords[1], @outer_coords[0], @outer_coords[1]);

	# add some text in red
	$p->setcolour("red");
	$p->setfont("Times-Roman", 20);
	$p->text({align => 'centre'}, CENTER_X,CENTER_Y, "Hello");

	# write the output to a file
	$p->output("file.ps");
}

sub coords_on_circle {
	my $angle = shift;
	my $radius = shift;

	return ( CENTER_X + ($radius*cos(($angle * PI)/180)), CENTER_Y + ($radius*sin(($angle * PI)/180)));
}

sub draw_circle_graph_from_file {
	my $datafile = shift;
	my $p = shift;
	#my $p = new PostScript::Simple( colour => 1, eps => 0, units => "bp", xsize => PS_X_SIZE, ysize => PS_Y_SIZE );
	my $OUTER_RADIUS = 387.5;
	my $INNER_RADIUS = 337.5;

	open my $F, "<$datafile" or die "$datafile failed to open\n";
	my $line = readline $F;
	my @positions, @differences;
	my $total_elems;
	my $max_diffs = 0;
	my @labels = split /\t/, $line;
	my $num_graphs = @labels-1;
	print "$num_graphs\n";

	# print legend
	$p->setfont("Helvetica", 12);
	if (@labels[1] =~ m/.*[:alpha:].*/) {
		for (my $i = 1; $i <= $num_graphs; $i++) {
			print "@labels[$i]\n";
			$p->setcolour($colours[$i-1][0],$colours[$i-1][1],$colours[$i-1][2]);
			my $max_height = ($num_graphs * 15) + 60;
			$p->text(10,($max_height - (15*$i)), "@labels[$i]");
		}
		$line = readline $F;
	}

	while ($line ne "") {
		my @items = split /\t/, $line;
		#$line =~ /(.+?)\t(.+?)$/;
		my $pos = shift @items;
		$total_elems = push (@positions, $pos);
		push (@differences, @items);
		$line = readline $F;
	}

	my @sorted = sort (@differences);
	my $diff_len = @sorted;
	$max_diffs = @sorted[@sorted-1];
	my $window_size = @positions[1]-@positions[0];
	my $circle_size = @positions[$total_elems-1];

	$p->setcolour(black);
	$p->setlinewidth(1);
	$p->circle(CENTER_X,CENTER_Y, $OUTER_RADIUS);
	$p->circle(CENTER_X,CENTER_Y, $INNER_RADIUS);
	$p->setlinewidth(2);

	for (my $j = 0; $j < $num_graphs; $j++) {
		my @coords = coords_on_circle(0,$INNER_RADIUS);
		my ($last_x, $last_y, $this_x, $this_y);
		$last_x = @coords[0];
		$last_y = @coords[1];

		$p->setcolour($colours[$j][0],$colours[$j][1],$colours[$j][2]);
		for (my $i = 0; $i < $total_elems; $i++) {
			my $angle = (@positions[$i]/$circle_size) * 360;
			#my $radius = $INNER_RADIUS + (($OUTER_RADIUS-$INNER_RADIUS)*(@differences[$i]/$max_diffs));
			my $radius = $INNER_RADIUS + (($OUTER_RADIUS-$INNER_RADIUS)*(@differences[($i*$num_graphs)+$j]/$max_diffs));
			my @new_coords = coords_on_circle($angle,$radius);
			$this_x = @new_coords[0];
			$this_y = @new_coords[1];
			$p->line($last_x, $last_y, $this_x, $this_y);
			$last_x = $this_x;
			$last_y = $this_y;
			print @positions[$i] . ": $this_x, $this_y\n";
		}
		$p->line($last_x, $last_y, @coords[0], @coords[1]);
	}

	$p->setfont("Helvetica", 6);
	$p->setcolour(black);

	for (my $i = 0; $i < $total_elems; $i++) {
		my $angle = (@positions[$i]/$circle_size) * 360;
		my $radius = $OUTER_RADIUS + 10;
		my @new_coords = coords_on_circle($angle,$radius);
		$p->text( {rotate => $angle}, @new_coords[0], @new_coords[1], "@positions[$i]");
	}

	$p->setfont("Helvetica", 12);
	$p->setcolour(black);
	$p->text(10, 10, "Maximum percent difference ($max_diffs) is scaled to 1");
	$p->text(10, 30, "Sliding window size of $window_size bp");
#	$p->output("$outfile.ps");

}

sub main_name_for_gb_feature {
	my $feat = shift;
	my @names;
	my $curr_name = "n/a";
	eval {@names = $feat->get_tag_values('locus_tag');} or @names="";
	if (@names) {
		$curr_name = @names[0];
	}
	eval {@names = $feat->get_tag_values('gene');} or return "$curr_name";
	return @names[0];
}

sub timestamp {
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
    $mon++;
    $mon = sprintf("%02d", $mon);
    $min = sprintf("%02d", $min);
    $sec = sprintf("%02d", $sec);
    $hour = sprintf("%02d", $hour);
    $mday = sprintf("%02d", $mday);

    $year -= 100;
    my $time = "$hour:$min:$sec";
    my $date = "$year$mon$mday";
    return ($time, $date);
}

# must return 1 for the file overall.
1;
