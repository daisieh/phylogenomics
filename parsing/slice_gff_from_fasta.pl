use Data::Dumper;
use lib "$FindBin::Bin/../";
use Subfunctions qw (split_seq disambiguate_str get_iupac_code);

my $gff_file = shift;
my $gene = shift;
my $fastafile = shift;
my $outfile = shift;

open FH, "<", $gff_file;
my $gff_block = "";
my $line = readline FH;
while (defined $line) {
	# scaffold_99	phytozome9_0	gene	16787	19271	.	+	.	ID=Potri.T085300;Name=Potri.T085300
	if ($line =~ /gene.*$gene/) {
		$gff_block .= $line;
	} elsif ($gff_block ne "") {
		if ($line =~ /gene/) {
			last;
		}
		$gff_block .= $line;
	}
	$line = readline FH;
}
close FH;

if ($gff_block eq "") {
	print "No gene named $gene found.\n";
	exit;
}
my $gff_hash = parse_gff_block ($gff_block);
$gff_hash->{"sequence"} = subseq_from_fasta ($fastafile, $gff_hash->{"start"}, $gff_hash->{"end"});

open OUT_FH, ">", $outfile or die "couldn't create $outfile";
print OUT_FH ">$gff_hash->{Name}.gene\n$gff_hash->{sequence}\n";

my $params = { 'padded' => 0, 'separate' => 1 };
for (my $mRNA_num = 1; $mRNA_num <= (keys $gff_hash->{"mRNA"}); $mRNA_num++) {
	print "mRNA $mRNA_num:\n";
	my $seq = feature_to_seq ($gff_hash->{"sequence"}, $gff_hash->{"mRNA"}->{$mRNA_num}->{"five_prime_UTR"}, $params);
	print @$seq . " 5' UTR\n";
	if ((ref $seq) =~ /ARRAY/ ) {
		for (my $i=1; $i<=@$seq; $i++) {
			print "$i @$seq[$i-1]\n";
			print OUT_FH ">$gff_hash->{Name}.$mRNA_num.5UTR$i\n@$seq[$i-1]\n";
		}
	}
	$seq = feature_to_seq ($gff_hash->{"sequence"}, $gff_hash->{"mRNA"}->{$mRNA_num}->{"exon"}, $params);
	print @$seq . " exons\n";
	if ((ref $seq) =~ /ARRAY/ ) {
		for (my $i=1; $i<=@$seq; $i++) {
			print OUT_FH ">$gff_hash->{Name}.$mRNA_num.exon$i\n@$seq[$i-1]\n";
		}
	}
	$seq = feature_to_seq ($gff_hash->{"sequence"}, $gff_hash->{"mRNA"}->{$mRNA_num}->{"three_prime_UTR"}, $params);
	print @$seq . " 3' UTR\n";
	if ((ref $seq) =~ /ARRAY/ ) {
		for (my $i=1; $i<=@$seq; $i++) {
			print OUT_FH ">$gff_hash->{Name}.$mRNA_num.3UTR$i\n@$seq[$i-1]\n";
		}
	}
	$seq = feature_to_seq ($gff_hash->{"sequence"}, $gff_hash->{"mRNA"}->{$mRNA_num}->{"CDS"}, $params);
	print @$seq . " CDS\n";
	if ((ref $seq) =~ /ARRAY/ ) {
		for (my $i=1; $i<=@$seq; $i++) {
			print "CDS $i\n";
			print OUT_FH ">$gff_hash->{Name}.$mRNA_num.CDS.$i\n@$seq[$i-1]\n";
		}
	}
}

sub feature_to_seq {
	my $sequence = shift;
	my $feature = shift;
	my $params = shift;

	my ($padded, $separate) = 0;
	if ((ref $params) =~ /HASH/) {
		if (defined $params->{"padded"}) {
			$padded = $params->{"padded"};
		}
		if (defined $params->{"separate"}) {
			$separate = $params->{"separate"};
		}
		if (defined $params->{"offset"}) {

		}

	}

	my $finalseq = "";
	my @seqarray = ();
	if (! (defined $feature)) {
		return "";
	}
	for (my $i = 1; $i <= keys $feature; $i++) {
		my $feat = $feature->{$i};
		my ($startseq, $seq, $endseq) = split_seq ($sequence, $feat->{"start"}, $feat->{"end"});
		if ($padded == 1) {
			push @seqarray, "x" x length($startseq) . $seq . "x" x length($endseq);
		} else {
			push @seqarray, $seq;
		}
	}

	if ($separate == 0) {
		while ($seqarray[0] ne "") {
			my $currchars = "";
			foreach my $seq (@seqarray) {
				if ($seq =~ m/(.)(.*)$/) {
					$currchars .= $1;
					$seq = $2;
				}
			}
			$currchars =~ s/x//g;
			if ($currchars eq "") {
				$finalseq .= "N";
			} else {
				$finalseq .= get_iupac_code ($currchars);
			}
		}
		return \($finalseq);
	} else {
		return \@seqarray;
	}


}

sub subseq_from_fasta {
	my $fastafile = shift;
	my $start = shift;
	my $end = shift;

	my $sequence = "";
	my $pos = 0;
	my $newstart = 0;
	my $length = $end - $start;
	my $newend = 0;
	open FH, "<", $fastafile or die "couldn't open $fastafile";

	my $line = readline FH; # first line is the name
	$line = readline FH;
	while (defined $line) {
		$line =~ s/\s//g;
		my $linelen = length ($line);
		if (($pos + $linelen) >= $start) {
			if ($newstart == 0) {
				$newstart = $start - $pos;
				$newend = $newstart + $length;
			}
			$sequence .= "$line";
			if ($pos >= $end) {
				last;
			}
		}
		$pos += $linelen;

		$line = readline FH;
	}
	close FH;
	my (undef, $finalseq, undef) = split_seq ($sequence, $newstart, $newend);

	return $finalseq;
}


sub parse_gff_block {
	# must pass in a block corresponding to a single gene.
	# should create a dictionary where the "gene" key contains the actual full location info
	# and all subsequent type entries contain offset information.
	my $gff_block = shift;

	my @lines = split(/\n/, $gff_block);
	my ($seqid, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = 0;
	my $offset = 0;
	my $gff_hash = {};

	my $gene_line = shift @lines;
	# this first line is the gene line. We are going to need to use this to get the whole sequence and the offset info.
	# scaffold_99	phytozome9_0	gene	16787	19271	.	+	.	ID=Potri.T085300;Name=Potri.T085300
	($seqid, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = split(/\t/, $gene_line);
	if ($type =~ /gene/) {
		my $attr_hash = parse_attributes($attributes);
		$offset = $start - 1;
		$gff_hash->{"start"} = $start;
		$gff_hash->{"end"} = $end;
		$gff_hash->{"seqid"} = $seqid;
		$gff_hash->{"ID"} = delete $attr_hash->{"ID"};
		$gff_hash->{"Name"} = delete $attr_hash->{"Name"};
		$gff_hash->{"attributes"} = $attr_hash;
	}
	my $gene_name = $gff_hash->{"Name"};

	# read in the block and hash it:
	foreach my $line (@lines) {
		if ($line =~ /^\s*$/) {
			next;
		}
		($seqid, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = split(/\t/, $line);
		# the first two are the same as its parent gene, so can be discarded.
		# type is the key under which this should be filed, with a subkey for the unique id of this feature.
		#	the id is named after the parent: #PAC:27046906.exon.3
		# the start and end are offsets from the parent gene.
		# the next three are irrelevant for our purposes, so can be discarded.
		# keep the attributes as is.

		$start = $start - $offset;
		$end = $end - $offset;

		my $attr_hash = parse_attributes($attributes);
		my $id = delete $attr_hash->{"ID"};
		my $name = delete $attr_hash->{"Name"};

		my $hash_ptr = {};
		$hash_ptr->{"start"} = $start;
		$hash_ptr->{"end"} = $end;
		$hash_ptr->{"attributes"} = $attr_hash;


		if ($type eq "mRNA") {
			if ($name =~ /$gene_name\.(\d+)/) {
				$gff_hash->{"mRNA"}->{$id} = $hash_ptr;
# 				$gff_hash->{$type} = $hash_ptr;
				$hash_ptr->{"ID"} = $id;
				$hash_ptr->{"Name"} = $name;
			}
		} else {
			if ($id =~ /(.*?)\.$type\.(\d+)/) {
				$gff_hash->{"mRNA"}->{$1}->{$type}->{$2} = $hash_ptr;
			}
		}
	}
	# now that we've finished hashing, we can rename the mRNAs with the numerical index.
	my $mRNA_hash = {};
	foreach my $k (keys $gff_hash->{"mRNA"}) {
		my $name = $gff_hash->{"mRNA"}->{$k}->{"Name"};
		my $this_hash = delete $gff_hash->{"mRNA"}->{$k};
		$name =~ /$gene_name\.(\d+)/;
		$mRNA_hash->{$1} = $this_hash;
	}
	$gff_hash->{"mRNA"} = $mRNA_hash;

	return $gff_hash;
}

sub parse_attributes {
	my $attributes = shift;
	# ID=Potri.019G067600;Name=Potri.019G067600
	# ID=PAC:27027285;Name=Potri.T085300.1;pacid=27027285;longest=1;Parent=Potri.T085300
	my $attr_hash = {};
	my @attrs = split (/;/, $attributes);
	foreach my $attr (@attrs) {
		if ($attr =~ /(.*?)=(.*)/) {
			$attr_hash->{$1} = $2;
		}
	}
	return $attr_hash;
}
