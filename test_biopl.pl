require "subfuncs.pl";
use Bio::SeqIO;
use PostScript::Simple;
use Bio::Align::Utilities qw(cat);
use Bio::Tools::Run::Phylo::PAML::Yn00;


use constant CENTER_X => 600; # X coordinate of circle center
use constant CENTER_Y => 600; # Y coordinate of circle center
use constant PS_X_SIZE => 1200; # X size of the PostScript object
use constant PS_Y_SIZE => 1200; # Y size of the PostScript object

my $usage  = "test_biopl.pl \n";

my $gb_file = shift or die $usage;
my $fa_file = shift or die $usage;

my $whole_aln = make_aln_from_fasta_file ($fa_file);
my @gene_alns;

my $seqio_object = Bio::SeqIO->new(-file => $gb_file);
my $seq_object = $seqio_object->next_seq;

while ($seq_object) {
	for my $feat_object ($seq_object->get_SeqFeatures) {
		if ($feat_object->primary_tag eq "CDS") {
			my $name = main_name_for_gb_feature($feat_object);
			my @locations = $feat_object->location->each_Location;
			my $cat_aln = 0;
			my $strand = 0;
			foreach $loc (@locations) {
				$strand = $loc->strand;
				my $start = $loc->start;
				my $end = $loc->end;
				my $curr_slice = $whole_aln->slice($start, $end);
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

			$cat_aln = $cat_aln->slice(1, $cat_aln->length()-3);
			$cat_aln->description($name);
			push @gene_alns, $cat_aln;
		}
	}
	$seq_object = $seqio_object->next_seq;
}

my $x = @gene_alns;
print "gene\tseq1\tseq2\tka\tks\tka/ks\n";

foreach my $aln (@gene_alns) {
	my $yn = Bio::Tools::Run::Phylo::PAML::Yn00->new();
	my $perc_id = $aln->percentage_identity();
	my $name = $aln->description();
	$yn->alignment($aln);
	my ($rc,$parser) = $yn->run();
		foreach my $seq ($aln->each_seq()) {
			print $seq->display_name(), "\t", $seq->seq(), "\n";
		}
	if ($rc == 0) {
		my $t = $yn->error_string();
		print "problem in $name: $t\n";
		foreach my $seq ($aln->each_seq()) {
			print $seq->display_name(), "\t", $seq->seq(), "\n";
		}
	} else {
		while( my $result = $parser->next_result ) {
			my @otus = $result->get_seqs();
			my $MLmatrix = $result->get_MLmatrix();

#  			for (my $i=1; $i < scalar @otus; $i++) {
# 					my $seq1 = @otus[$i]->display_name();
# 					print "$seq1...\n";
# 			}
			for (my $i=1; $i < scalar @otus; $i++) {
				for (my $j=1; $j<scalar @otus; $j++) {
# 				for (my $j=(scalar @otus) - 1; $j > $i; $j--) {
					my $dN = $MLmatrix->[$i]->[$j]->{dN};
					my $dS = $MLmatrix->[$i]->[$j]->{dS};
					my $kaks =$MLmatrix->[$i]->[$j]->{omega};
					my $seq1 = @otus[$i]->display_name();
					my $seq2 = @otus[$j]->display_name();
					print "$name\t$seq1\t$seq2\t$dN\t$dS\t$kaks\n";
				}
			}
		}
	}
}

