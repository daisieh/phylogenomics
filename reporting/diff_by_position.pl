my $inputfile = shift;

open(fileIN, "<", "$inputfile")  or die "no file named $inputfile";

my $input = readline fileIN;
my $length = 0;
my $curr_seq_name = "";
my @seq_names = ();
my @seqs = ();
while ($input ne "") {
	if ($input =~ /^>(.+)\s*/) {
		if ($length > 0) {
			# we are at the next taxon; push the last one onto the taxon array.
			push @seq_names, $curr_seq_name;
			push @seqs, $sequence;
			$length = 0;
			$sequence = "";
		}
		$curr_seq_name = $1;
	} else {
		$input =~ /^\s*(.+)\s*/;
		$sequence .= $1;
		$length += length($1);
	}
	$input = readline fileIN;
}

if ($length > 0) {
	push @seq_names, $curr_seq_name;
	push @seqs, $sequence;
}

print "pushed " . @seq_names . " sequences\n";
for (my $i=0; $i<(@seq_names);) {
	my $diffprot1 = uc(@seqs[$i]);
	my $diffseq1 = uc(@seqs[$i+1]);
	my $diffprot2 = uc(@seqs[$i+2]);
	my $diffseq2 = uc(@seqs[$i+3]);

	@seq_names[$i] =~ /(.+?)\s+(.*)/;
	my $seq1 = $1;
	my $seq_type = $2;
	@seq_names[$i+2] =~ /(.+?)\s+(.*)/;
	my $seq2 = $1;

	my $is_prot = 0;

	print "comparing $seq1 to $seq2, $seq_type\n";
	if ($seq_type =~ /prot/) {
		$is_prot = 1;
	}

	for (my $j=0; $j<length($diffprot1); $j++) {
		my $diff_type = "";

		$diffprot1 =~ /(.)(.*)/;
		my $aa1 = $1;
		$diffprot1 = $2;

		$diffprot2 =~ /(.)(.*)/;
		my $aa2 = $1;
		$diffprot2 = $2;

		my $start = ($j*3)+1;
		my $end = $start + 2;

		$diffseq1 =~ /(.{3})(.*)/;
		my $codon1 = $1;
		$diffseq1 = $2;

		$diffseq2 =~ /(.{3})(.*)/;
		my $codon2 = $1;
		$diffseq2 = $2;

		if ($codon1 ne $codon2) {
			if ($aa1 eq $aa2) {
				$diff_type = "S";
			} else {
				$diff_type = "NS";
			}
			print "$diff_type\t$start-$end\t$aa1\t$aa2\t$codon1\t$codon2\n";
		}
	}
	$i = $i + 4;
}

close (fileIN);
