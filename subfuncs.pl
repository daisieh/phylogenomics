use PostScript::Simple;
use Bio::AlignIO;
use Bio::SeqIO;

use constant PI 	=> 3.1415926535897932384626433832795;
use constant MAX_LABEL_LENGTH => 16; # if the gb file contains labels longer than this labelsOK will return FALSE
use constant MIN_INNER => 1.4; # minimal angle required for non-overlapping inner labels
use constant MIN_OUTER => 1.1; # minimal angle required for non-overlapping outer labels
use constant MIN_OPT_OUTER => 1.2; # minimal angle required for outwards-shifted labels
use constant MIN_OPT_INNER => 1.8; # minimal angle required for inwards-shifted labels
use constant MAX_OPT_RUNS => 3; # emergency brake for optimization loop
use constant LABEL_SHIFT => 3.9; # shift for centering of labels above features
use constant MIN_FRAMED_FEATURE_WIDTH => 0.3; # minimal width for a feature to be drawn with frame
use constant CENTER_X => 600; # X coordinate of circle center500
use constant CENTER_Y => 600; # Y coordinate of circle center
use constant PS_X_SIZE => 1200; # X size of the PostScript object1000
use constant PS_Y_SIZE => 1200; # Y size of the PostScript object
use constant CHAR_WIDTH			=> 10; # median char width - - NOT CORRECT
use constant CHAR_HEIGHT		=> 10; # Character height - - has to be checked again

#
# graphics default values
#
my $MAIN_RADIUS				= 387.5;
my $RE_LABEL_RADIUS 		= $MAIN_RADIUS + 127.5; # Radius for RE labels
my $OPT_OUTER_RADIUS 		= $MAIN_RADIUS + 67.5; # radius for outwards-shifted labels
my $OPT_INNER_RADIUS 		= $MAIN_RADIUS - 67.5; # radius for inwards-shifted labels
my $INSIDE_FEATURE_RADIUS 	= $MAIN_RADIUS - 12.15;
my $OUTSIDE_FEATURE_RADIUS 	= $MAIN_RADIUS + 12.15;
my $INSIDE_LABEL_RADIUS		= $MAIN_RADIUS - 27.5;
my $OUTSIDE_LABEL_RADIUS	= $MAIN_RADIUS + 27.5;

sub convert_aln_to_nexus {
	my $aln = shift;

	my $nexblock = "";
	my $result = "";
	my $ntax = 0;
	my $nchar = 0;
	foreach my $seq ( $aln->each_seq()) {
		my $len = $seq->length;
 		$nexblock .= "'" . $seq->display_name . "'\t";
 		my $seq_str = $seq->seq();
 		$nexblock .= "$seq_str\n";
 		$nchar = length($seq_str);
 		$ntax++;
	}
	$result .= "#NEXUS\n\nBegin DATA;\nDimensions ntax=$ntax nchar=$nchar;\n";
	$result .= "Format datatype=dna gap=-;\n";
 	$result .= "Matrix\n$nexblock\n;\nEnd;\n";

	return $result;
}

sub return_partition_2 {
	my $inseq = shift;
	my $start_pos = shift;
	my $stop_pos = shift;
	my $newaln = shift;

	my $seq;
	eval {$seq = $inseq->next_seq;} or die "3 not fasta\n";
	my $nexblock = "";
	my $result = "";
	my $ntax = 0;
	my $nchar = 0;
	while ($seq ne "") {
		my $len = $seq->length;
		if ($stop_pos > $len) {
			return 0;
		}
 		my $seq_seg = $seq->subseq($start_pos, $stop_pos);
 		my $outseq = Bio::LocatableSeq->new(-seq => $seq_seg, -id => $seq->display_name);
 		$newaln->add_seq ($outseq);

 		$seq = $inseq->next_seq;
	}

	my $p = $newaln->percentage_identity();
	print "$start_pos\t$stop_pos\t$p\n";

	return 1;
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

	my $inseq = Bio::SeqIO->new(-file => "<$fa_file", -format => "fasta");
	my $newaln = Bio::SimpleAlign->new();

	my $seq;
	eval {$seq = $inseq->next_seq;} or die "not fasta\n";
	while ($seq ne "") {
 		my $outseq = Bio::LocatableSeq->new(-seq => $seq->seq(), -id => $seq->display_name);
 		$newaln->add_seq ($outseq);
 		$seq = $inseq->next_seq;
	}
	return $newaln;
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
	my $OUTER_RADIUS = 387.5;
	my $INNER_RADIUS = 337.5;

	open my $F, "<$datafile" or die "$datafile failed to open\n";
	my $line = readline $F;
	my @positions, @differences;
	my $total_elems;
	my $max_diffs = 0;
	while ($line ne "") {
		$line =~ /(.+?)\t(.+?)$/;
		$total_elems = push (@positions, $1);
		my $curr_diffs = $2;
		push (@differences, $curr_diffs);
		if ($max_diffs < $curr_diffs) {
			$max_diffs = $curr_diffs;
		}
		$line = readline $F;
	}

	my $circle_size = @positions[$total_elems-1];

	my $p = new PostScript::Simple( colour => 1, eps => 0, units => "bp", xsize => PS_X_SIZE, ysize => PS_Y_SIZE );
	$p->newpage;

	$p->circle(CENTER_X,CENTER_Y, $OUTER_RADIUS);
	$p->circle(CENTER_X,CENTER_Y, $INNER_RADIUS);

	my @coords = coords_on_circle(0,$INNER_RADIUS);
	my ($last_x, $last_y, $this_x, $this_y);
	$last_x = @coords[0];
	$last_y = @coords[1];

	$p->setcolour(red);
	for (my $i = 0; $i < $total_elems; $i++) {
		my $angle = (@positions[$i]/$circle_size) * 360;
		my $radius = $INNER_RADIUS + (($OUTER_RADIUS-$INNER_RADIUS)*(@differences[$i]/$max_diffs));
		print "$angle, $radius\n";
		my @new_coords = coords_on_circle($angle,$radius);
		$this_x = @new_coords[0];
		$this_y = @new_coords[1];
		$p->line($last_x, $last_y, $this_x, $this_y);
		$last_x = $this_x;
		$last_y = $this_y;
	}
	$p->line($last_x, $last_y, @coords[0], @coords[1]);
	$p->setfont("Helvetica", 12);
	$p->text(10, 10, "Maximum percent difference ($max_diffs) is scaled to 1");

	# write the output to a file
	$p->output("file.ps");
}


1;
