#!/usr/bin/env perl

use FindBin;
use lib "$FindBin::Bin/../..";
use Subfunctions;
use Bio::SeqIO;
use PostScript::Simple;
use Bio::Align::Utilities qw(cat);
use Bio::Tools::Run::Phylo::PAML::Yn00;
use Bio::Tools::Run::Phylo::PAML::Codeml;

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

my $paml_exec = Bio::Tools::Run::Phylo::PAML::Codeml->new
			   ( -params => { 'runmode' => -2,
							  'seqtype' => 1,
							} );
#my $paml_exec = Bio::Tools::Run::Phylo::PAML::Yn00->new();

#print "gene\tseq1\tseq2\tdN/dS\tdN\tdS\n";


foreach my $aln (@gene_alns) {
	my $name = $aln->description();
	my $resultstr = $name;
	$paml_exec->alignment($aln);
	my ($rc,$parser) = $paml_exec->run();
	if ($rc == 0) {
		my $t = $paml_exec->error_string();
# 		print "problem in $name: $t\n";
# 		foreach my $seq ($aln->each_seq()) {
# 			print $seq->display_name(), "\t", $seq->seq(), "\n";
# 		}
	} else {
		while( my $result = $parser->next_result ) {
			my @otus = $result->get_seqs();
			my $MLmatrix = $result->get_MLmatrix();
			#	this loop set is excessively complicated because I am trying to get the output to correspond to yn00's output block.
			my $j = 1;
			for (my $i=1;$i<=$j+1;$i++) {
				if ($i == scalar @otus) { last; }
				for ($j=0;$j<$i;$j++) {
					my $dN = $MLmatrix->[$j]->[$i]->{dN};
					my $dS = $MLmatrix->[$j]->[$i]->{dS};
					my $kaks =$MLmatrix->[$j]->[$i]->{omega};
					my $seq1 = @otus[$i]->display_name();
					my $seq2 = @otus[$j]->display_name();
					$resultstr .= "\t$seq1#$seq2#" . $kaks;
					#print "$name\t$seq1\t$seq2\t$kaks\t$dN\t$dS\n";
				}
			}
		}
	}
	print $resultstr, "\n";
}


