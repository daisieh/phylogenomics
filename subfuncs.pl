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

sub return_partition {
	my $inseq = shift;
	my $start_pos = shift;
	my $stop_pos = shift;
	my $newaln = shift;

	my $seq;
	eval {$seq = $inseq->next_seq;} or die "not fasta\n";
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
		return $p;
	} else {
		return 0;
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

1;
