require "subfuncs.pl";
use Bio::SeqIO;
use PostScript::Simple;
use Bio::Align::Utilities qw(cat);
use Bio::Tools::Run::Phylo::PAML::Codeml;


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
				print ("$strand\n");
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
print "$x genes in alignment\n";

#### need to read PAML documentation to understand inputs and outputs: http://abacus.gene.ucl.ac.uk/software/pamlDOC.pdf

my $kaks_factory = Bio::Tools::Run::Phylo::PAML::Codeml->new
	(-verbose => $verbose,
	 -params => { 'runmode' => -2,
		      'seqtype' => 1,
		  }
	 );

foreach my $dna_aln (@gene_alns) {
	my $name = $dna_aln->description();
	my $temp = $dna_aln->get_seq_by_pos(1)->seq();
# 	print "$name, $temp\n";
	$kaks_factory->alignment($dna_aln);
	my ($rc,$parser) = $kaks_factory->run();
	foreach my $result ($parser->next_result()) {
		my $MLmatrix = $result->get_MLmatrix();
		my @otus = $result->get_seqs();
		my $size = @otus;
		my $temp = scalar ($MLmatrix->[0]);
		print "$temp\n";

		print join("\t", qw(Ka Ks Ka/Ks)), "\n";
		for( my $i = 0; $i < ($size-1) ; $i++) {
			for( my $j = $i+1; $j < ($size); $j++ ) {
				#my $sub_aa_aln = $aa_aln->select_noncont($pos[$i],$pos[$j]);
				#my $sub_dna_aln = $dna_aln->select_noncont($pos[$i],$pos[$j]);
				print "$i $j\n";
				print join("\t",
					   #$otus[$i]->display_id,
					   #$otus[$j]->display_id,
					   $MLmatrix->[$i]->[$j]->{'dN'},
					   $MLmatrix->[$i]->[$j]->{'dS'},
					   $MLmatrix->[$i]->[$j]->{'omega'},
					   #sprintf("%.2f",$sub_aa_aln->percentage_identity),
					   #sprintf("%.2f",$sub_dna_aln->percentage_identity),
					   ), "\n";
			}
		}
	}
}
