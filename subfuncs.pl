use Bio::AlignIO;
use Bio::SeqIO;
use Bio::Align::Utilities qw(cat);

=head1

B<Arbitrary collection of helper subfunctions that don't have another home.>

=cut


=head1

B<String $nexus_str convert_aln_to_nexus ( SimpleAlign $aln )>

Takes a SimpleAlign object and returns a NEXUS-formatted string representing the alignment.

=cut

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


=head1

B<String $genes_str get_locations_from_genbank_file ( String $genbank_file, <optional> String $type )>

Takes a genbank-formatted file and returns a tab-delimited list of genes and positions.
The final entry in the list will be the size of the source genbank file.

$genbank_file is the name of the genbank file to parse
$type is the tag of interest (will be "gene" if this is not specified)

=cut

sub get_locations_from_genbank_file {
    my $gb_file = shift;
    my $type = shift;
    my @gene_alns;

    # default type is "gene"
    unless ($type) { $type = "gene"; }

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
            if ($feat_object->primary_tag eq $type) {
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

=head1

B<SimpleAlign @gene_alns parse_aln_into_genes ( SimpleAlign $aln, String $genbank_file, <optional> String $type )>

Parses a SimpleAlign into an array of SimpleAligns, using a genbank-formatted file to determine gene positions.

$genbank_file is the name of the genbank file to parse
$type is the tag of interest (will be "gene" if this is not specified)

=cut

sub parse_aln_into_genes {
	my $whole_aln = shift;
    my $gb_file = shift;
    my $type = shift;
    my @gene_alns;

    # default type is "gene"
    unless ($type) { $type = "gene"; }

    my $gb_seqio = Bio::SeqIO->new(-file => $gb_file);
    my $start, $end;

	while (my $seq_object = $gb_seqio->next_seq) {
		for my $feat_object ($seq_object->get_SeqFeatures) {
			if ($feat_object->primary_tag eq $type) {
				my $name = main_name_for_gb_feature($feat_object);
				my @locations = $feat_object->location->each_Location;
				my $cat_aln = 0;
				my $strand = 0;
				foreach $loc (@locations) {
					$strand = $loc->strand;
					my $start = $loc->start;
					my $end = $loc->end;
					my $curr_slice = $whole_aln->slice($start, $end,1);
					if ($cat_aln == 0) {
						$cat_aln = $curr_slice;
					} else {
						$cat_aln = cat($cat_aln, $curr_slice);
					}
				}
				if ($strand < 0) {
					# must flip each seq in the curr_slice
					my $flipped_aln = Bio::SimpleAlign->new();
					foreach $seq ( $cat_aln->each_seq() ) {
						$seq = $seq->revcom();
						$flipped_aln->add_seq($seq);
					}
					$cat_aln = $flipped_aln;
				}

				$cat_aln = $cat_aln->slice(1, $cat_aln->length()-3,1);
				$cat_aln->description($name);
				push @gene_alns, $cat_aln;
			}
		}
	}
    return \@gene_alns;
}



=head1

B<Int $percent_diff perc_diff_partition ( SimpleAlign $aln, Int $start, Int $end )>

Finds the overall percentage difference for the specified region of the SimpleAlign object.

=cut

sub perc_diff_partition {
	my $aln = shift;
	my $start_pos = shift;
	my $stop_pos = shift;
	my $check_for_num_seqs = shift;
	my $delete_gaps = shift;

	if ($stop_pos < $aln->length()) {
		my $aln_slice = $aln->slice($start_pos, $stop_pos);

		# remove any sequences with gaps
		if ($delete_gaps) {
			foreach my $seq ($aln_slice->each_seq()) {
				my $seqstr = $seq->seq();
				if ($seqstr =~ /-|n|N|\?/) {
					$aln_slice->remove_seq($seq);
				}
			}
		}

		# check to see if one of the sequences was excluded.
		if ($aln->num_sequences != $aln_slice->num_sequences) {
            # if there were only two seqs in the first place, there won't be a valid number here.
            if ($aln->num_sequences == 2) {
                return -1;
            } else {
                my $p = $aln_slice->percentage_identity();
                if ($check_for_num_seqs) {
                    # return the negative of the perc identity to denote that a seq was excluded
#                     return -(100-$p);
		return -20000;
                } else {
                    return 100-$p;
                }
            }
        } else {
            my $p = $aln_slice->percentage_identity();
            return 100-$p;
        }
	} else {
		return -200;
	}
}


=head1

B<SimpleAlign $fa_aln make_aln_from_fasta_file ( String $fasta_file )>

Takes a fasta-formatted file and returns a SimpleAlign object.
The sequences in the fasta file must already be aligned.

=cut

sub make_aln_from_fasta_file {
	my $fa_file = shift;
	my $flush = shift;
	my $min_length = 0;

	if ($flush eq "") { # unless otherwise specified, return a flush alignment
		$flush = 1;
	}

	my $inseq = Bio::SeqIO->new(-file => "<$fa_file", -format => "fasta");
	my $newaln = Bio::SimpleAlign->new();

	my $seq;
	eval {$seq = $inseq->next_seq;} or die "file not in fasta format.\n";
	$min_length = $seq->length();
	while ($seq ne "") {
		my $name = $seq->display_name;
		#remove weird chars, file suffixes
		$name =~ s/(.*?)\..*/$1/;
 		my $outseq = Bio::LocatableSeq->new(-seq => $seq->seq(), -id => $name);
 		$newaln->add_seq ($outseq);
 		if ($min_length > $seq->length() ) {
 			$min_length = $seq->length();
 		}
 		$seq = $inseq->next_seq;
	}

	if ($flush) {
		my $flush_aln = $newaln->slice(1,$min_length);
		return $flush_aln;
	}
	return $newaln;
}


=head1

B<String $name main_name_for_gb_feature ( SeqFeatureI $feature )>

Finds the main name for the given SeqFeatureI from a BioPerl-parsed Genbank object.

=cut

sub main_name_for_gb_feature {
	my $feat = shift;
	my @names;
	my $curr_name = "n/a";

	if ($feat->has_tag('locus_tag')) {
        @names = $feat->get_tag_values('locus_tag');
		$curr_name = @names[0];
	}

	if ($feat->has_tag('gene')) {
        @names = $feat->get_tag_values('gene');
		$curr_name = @names[0];
	}

	return $curr_name;
}


=head1

B<(String $time, String $date) timestamp ()>

Convenience function to get the current time and date as formatted strings.

=cut

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

=head1

B<String $result_str combine_files ( Array \@files, Boolean $has_names, Boolean $has_header)>

Takes any number of input files containing tab-delimited lists of the same length
and creates a tab-delimited string with each list as a column.

@files:     a pointer to a list of filenames
$has_names: 1 if files have the same column of names (first column of each file), otherwise 0.
$has_header:1 if files have column labels in first row, otherwise 0.

Default is no column names, no common row names.

=cut

sub combine_files {
    my $fileptr = shift;
    my $has_names = shift;
    my $has_header = shift;
    my $out_file = shift;

    my @files = @$fileptr;

    if (@files < 1) { die "no files provided."; }

    my @inputs;
    for (my $i=0; $i<@files; $i++) {
        open FH, "<", @files[$i] or die "can't open @files[$i]\n";
        my @data = <FH>;
        close FH;
        push @inputs, \@data;
    }

    my $result = "";
    my @labels = ();
    my $num_entries = scalar @{@inputs[0]};
    for (my $i = 0; $i < @files; $i++) {
        if (scalar @{@inputs[$i]} != $num_entries) {
            die "Error: files have different numbers of inputs." . scalar @{@inputs[$i]};
        }
        if ($has_header) {
            my @heads = split /\t/, (@{@inputs[$i]}[0]);
#             foreach $head (@heads) {
#                 $head = @files[$i] . "|" . $head;
#             }
            @{@inputs[$i]}[0] = join ("\t", @heads);
        }
    }

    #if different files have same column names, must disambiguate column names.
#     my @sortedheads = sort(@heads);
#     my $flag = 0;
#     my $common_name = @sortedheads[0];
#     for(my $j=1; $j<@sortedheads; $j++) {
#         if ($common_name eq @sortedheads[$j]) {
#
#         }
#     }
#

    for (my $j = 0; $j < $num_entries; $j++) {
        if ($has_names) {
            my $entry = @{@inputs[0]}[$j];
            $entry =~ /(.+?)\t(.*)/;
            $result .= "$1\t";
        }
        for (my $i = 0; $i < @files; $i++) {
            my $entry = @{@inputs[$i]}[$j];
            $entry =~ s/\n//g;
            if ($has_names) {
                $entry =~ /(.+?)\t(.*)/;
                $entry = $2;
            }
            $result .= $entry . "\t";
        }
        if ($has_header && $has_names && ($j==0)) {
            #clean up the header row
            $result =~ s/^.*?\|//;
        }
        $result .= "\n";
    }
    return $result;
}

=head1

B<slice_fasta_to_exons ( String $fa_file, String $gb_file, String $out_file )>

Given a fasta file of aligned sequences and a corresponding genbank file
with the CDS coordinates, will create a fasta file with each
CDS corresponding to a separate sequence block.

$fa_file:   fasta file of aligned sequences
$gb_file:   genbank file with CDS coordinates
$out_file:  prefix of output files

=cut

sub slice_fasta_to_exons {
    my $fa_file = shift;
    my $gb_file = shift;
    my $out_file = shift;

    my $whole_aln = make_aln_from_fasta_file ($fa_file);
    my @gene_alns;

    my $result_str = get_locations_from_genbank_file ($gb_file, "CDS");
    my @exons = split (/\n/,$result_str);

    # we don't need the last line; it just says the size of the gb source sequence
    while ({pop @exons} =~ /source/) {
    }

    foreach my $exon (@exons) {
        $exon =~ /(.+?)\t(.+?)\t(.+)\s*.*/;
        my $name = $1;
        my $start = $2;
        my $end = $3;
        my $curr_slice = $whole_aln->slice($start, $end);
        $curr_slice->description($name);
        push @gene_alns, $curr_slice;
    }

    open my $gene_file, ">$out_file.fasta";

    foreach my $aln (@gene_alns) {
        my $gene_name = $aln->description();
        foreach my $seq ($aln->each_seq()) {
            my $name = $seq->id() . "_$gene_name: " . $seq->length();
            print $gene_file ">$name\n";
            print $gene_file $seq->seq() . "\n";
        }
    }
    close $gene_file;
}

=head1

B<Hashref make_label_lookup ( String $labelfile )>

Given a tab-delimited file of sample ids and human-readable labels, returns
a hash ref for quick lookup.

$labelfile:   tab-delimited file of sample ids and human-readable labels

=cut

sub make_label_lookup {
    my $labelfile = shift;
    my %labels;
    if ($labelfile) {
        open FH, "<", "$labelfile" or die "make_label_lookup died: couldn't open $labelfile\n";
        my @items = <FH>;
        close FH;
        foreach my $line (@items) {
            (my $name, my $label) = split (/\t/,$line);
            $label =~ s/\r|\n//;
            $labels{$name} = $label;
        }
    }
    return \%labels;
}

# must return 1 for the file overall.
1;
