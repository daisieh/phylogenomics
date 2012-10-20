use Bio::AlignIO;
use Bio::SeqIO;
use Bio::Align::Utilities qw(cat);



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

sub parse_genbank_file {
    my $gb_file = shift;
    my @gene_alns;

    my $seqio_object = Bio::SeqIO->new(-file => $gb_file);
    my $seq_object = $seqio_object->next_seq;
    my $source_start, $source_end;

    my $result_str = "";
    while ($seq_object) {
        for my $feat_object ($seq_object->get_SeqFeatures) {
            if ($feat_object->primary_tag eq "source") {
                my @locations = $feat_object->location->each_Location;
                $source_start = @locations[0]->start;
                $source_end = @locations[0]->end;
            }
            if ($feat_object->primary_tag eq "CDS") {
                my $name = main_name_for_gb_feature($feat_object);
                my @locations = $feat_object->location->each_Location;
                my $exon = 1;
                foreach my $loc (@locations) {
                    my $exon_name = $name . "_" . $exon++;
                    my $start = $loc->start;
                    my $end = $loc->end;
                    if (@locations == 1) { $exon_name = $name; }
                    $result_str .= "$exon_name\t$start\t$end\n";
                }
            }
        }
        $seq_object = $seqio_object->next_seq;
    }
    $result_str .= "source\t$source_start\t$source_end\n";
    return "$result_str";
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
