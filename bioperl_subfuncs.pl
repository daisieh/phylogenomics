use Bio::AlignIO;
use Bio::SeqIO;
use Bio::Align::Utilities qw(cat);
use FindBin;
use lib "$FindBin::Bin/..";
use Nexus qw(write_nexus_character_block write_nexus_taxa_block);

=head1

B<Arbitrary collection of helper subfunctions that don't have another home.>

=cut


=head1

B<String $nexus_str convert_aln_to_nexus ( SimpleAlign $aln )>

Takes a SimpleAlign object and returns a NEXUS-formatted string representing the alignment.

=cut

sub convert_aln_to_nexus {
	my $aln = shift;

	my $taxa_hash = {};
	my @taxa_names = ();

	foreach my $seq ( $aln->each_seq()) {
		my $name = $seq->display_name;
		$name =~ s/-/_/g;

		$taxa_hash->{$name} = $seq->seq();
		push @taxa_names, $name;
	}
	my $result = "#NEXUS\n\n";

	$result .= write_nexus_taxa_block (\@taxa_names);
	$result .= write_nexus_character_block ($taxa_hash, \@taxa_names);
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
                my $name = main_name_for_gb_feature($feat_object, 0);
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
The exons will be rev-comped if necessary.

$genbank_file is the name of the genbank file to parse
$type is the tag of interest (will be "gene" if this is not specified)

=cut

sub parse_aln_into_genes {
	my $whole_aln = shift;
    my $gb_file = shift;
    my $remove_stop = shift;

    my @gene_alns;

    my $type = "CDS";

    my $gb_seqio = Bio::SeqIO->new(-file => $gb_file);
    my $start, $end;

	while (my $seq_object = $gb_seqio->next_seq) {
		for my $feat_object ($seq_object->get_SeqFeatures) {
			if ($feat_object->primary_tag eq $type) {
				my @loc_list = ();
				my $name = main_name_for_gb_feature($feat_object, 1);
				my @locations = $feat_object->location->each_Location;
				my $cat_aln = 0;
				my $strand = 0;
				foreach $loc (@locations) {
					$strand = $loc->strand;
					my $start = $loc->start;
					my $end = $loc->end;
					push @loc_list, "$start..$end";
					my $curr_slice = $whole_aln->slice($start, $end,1);
					print "slice\n";
					if ($cat_aln == 0) {
						$cat_aln = $curr_slice;
					} else {
					print "cat\n";
						$cat_aln = cat($cat_aln, $curr_slice);
						print "cat2\n";
					}
					print "CDS\t$name\t$start\t$end\n";
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
				if ($remove_stop) {
					$cat_aln = $cat_aln->slice(1, $cat_aln->length()-3,1);
				}
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
		my $aln_slice = $aln->slice($start_pos, $stop_pos, 1);

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
		my $flush_aln = $newaln->slice(1,$min_length,1);
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
	my $search_others = shift;
	my @names;
	my $curr_name = "";

	if ($feat->has_tag('gene')) {
        @names = $feat->get_tag_values('gene');
		$curr_name = @names[0];
	}
	if (($curr_name eq "") && ($search_others)) {
		if ($feat->has_tag('locus_tag')) {
			@names = $feat->get_tag_values('locus_tag');
			$curr_name = @names[0];
		}
	}
	return $curr_name;
}


=head1

B<@SimpleAlign slice_fasta_to_exons ( String $fa_file, String $gb_file )>

Given a fasta file of aligned sequences and a corresponding genbank file
with the CDS coordinates, returns an array of SimpleAligns corresponding to each CDS.

$fa_file:   fasta file of aligned sequences
$gb_file:   genbank file with CDS coordinates

=cut

sub slice_fasta_to_exons {
    my $fa_file = shift;
    my $gb_file = shift;
    return slice_fasta_from_genbank_file ($fa_file, $gb_file, "CDS");
}

=head1

B<@SimpleAlign slice_fasta_from_genbank_file ( String $fa_file, String $gb_file, String $type )>

Given a fasta file of aligned sequences and a corresponding genbank file
with locus coordinates, returns an array of SimpleAligns corresponding to each locus.

$fa_file:   fasta file of aligned sequences
$gb_file:   genbank file with locus coordinates
$type:  locus type from genbank file (default is "gene")

=cut

sub slice_fasta_from_genbank_file {
    my $fa_file = shift;
    my $gb_file = shift;
    my $type = shift;

    unless ($type) { $type = "gene"; }

    my $whole_aln = make_aln_from_fasta_file ($fa_file);
    my @gene_alns;

    my $result_str = get_locations_from_genbank_file ($gb_file, $type);
    my @exons = split (/\n/,$result_str);

    # we don't need the last line; it just says the size of the gb source sequence
    while ({pop @exons} =~ /source/) {
    }

    foreach my $exon (@exons) {
        $exon =~ /(.+?)\t(.+?)\t(.+)\s*.*/;
        my $name = $1;
        my $start = $2;
        my $end = $3;
        my $curr_slice = $whole_aln->slice($start, $end, 1);
        $curr_slice->description($name);
        push @gene_alns, $curr_slice;
    }

	return \@gene_alns;
}

# must return 1 for the file overall.
1;

