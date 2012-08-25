use Bio::SeqIO;
use Bio::Align::Utilities qw(cat);

require "subfuncs.pl";

# Bring in the file and format, or die with a nice
# usage statement if one or both arguments are missing.
my $usage  = "test_biopl.pl genbank_file fasta_file result_dir\n";
my $gb_file   = shift or die $usage;
my $fa_file = shift or die $usage;
my $result_dir = shift or die $usage;

my $whole_aln = make_aln_from_fasta_file ($fa_file);
my @gene_alns;

my $seqio_object = Bio::SeqIO->new(-file => $gb_file);
my $seq_object = $seqio_object->next_seq;
print $seq_object;

my $result_str = "";
while ($seq_object) {
	for my $feat_object ($seq_object->get_SeqFeatures) {
		print $feat_object->primary_tag;
		if ($feat_object->primary_tag eq "CDS") {
			my $name = main_name_for_gb_feature($feat_object);
			my @locations = $feat_object->location->each_Location;
			my $cat_aln = 0;
			my $strand = 0;
			my $last_end = 0;
			foreach $loc (@locations) {
				$strand = $loc->strand;
				my $start = $loc->start;
				my $end = $loc->end;
				$last_end = $end;
				my $curr_slice = $whole_aln->slice($start, $end);
				if ($cat_aln == 0) {
					$cat_aln = $curr_slice;
				} else {
					$cat_aln = cat($cat_aln, $curr_slice);
				}
				if ($result_str eq "") { $result_str = "$name\t$start"; }
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

			$cat_aln = $cat_aln->slice(1, $cat_aln->length()-3);
			$cat_aln->description($name);
			push @gene_alns, $cat_aln;
			print "$result_str\t$last_end\n";
			$result_str = "";
		}
	}
	$seq_object = $seqio_object->next_seq;
}

#open my $gene_file, ">", "$result_dir\/$gb_file.fa";
foreach my $aln (@gene_alns) {
	foreach my $seq ( $aln->each_seq()) {
		print ">" . $aln->description . "\n" . $seq->seq() . "\n";
	}
}

# foreach my $aln (@gene_alns) {
# 	my $gene_name = $aln->description();
# 	#print "writing $gene_name...\n";
# #	my $outfile = "$result_dir\/$gene_name.nex";
# #	my $result = convert_aln_to_nexus ($aln);
# #	open my $gene_file, ">$outfile";
# #	truncate $gene_file, 0;
# 	$result =
# 	print $gene_file $result;
# #	close $gene_file;
# }
