require "subfuncs.pl";
use Bio::SeqIO;
use Bio::TreeIO;
use Bio::Root::Test;
use PostScript::Simple;
use Bio::Align::Utilities qw(cat);
use Bio::Tools::Run::Phylo::PAML::Codeml;
use File::Basename;
use Getopt::Long;

my $usage = "perl " . basename($0);
$usage .= " gb_file fa_file tree_file output_name\n";

my $gb_file = shift or die $usage;
my $fa_file = shift or die $usage;
my $tree_file = shift or die $usage;
my $output_name = shift or die $usage;
#GetOptions ('query=s' => \$queryfile, 'subject=s' => \$subjectfile);


my $whole_aln = make_aln_from_fasta_file ($fa_file);
my @gene_alns;

my $seqio_object = Bio::SeqIO->new(-file => $gb_file);
my $seq_object = $seqio_object->next_seq;

while ($seq_object) {
	for my $feat_object ($seq_object->get_SeqFeatures) {
		if ($feat_object->primary_tag eq "CDS") {
			my $name = main_name_for_gb_feature($feat_object);
			# no point dealing with this feature if it doesn't have a name...
			if ($name eq "") { next; }
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

my $paml_exec =	Bio::Tools::Run::Phylo::PAML::Codeml->new
#				PAML with fixed omega value (model=0 designates single omega)
#				( -params => { 'runmode' => 0, 'seqtype' => 1, 'model' => 0, 'fix_omega' => 1, 'omega' => 1 } );

#				PAML with free omegas (model=1)
				( -params => { 'runmode' => 0, 'seqtype' => 1, 'model' => 1, 'fix_blength' => 1 });
#				( -params => { 'runmode' => 0, 'seqtype' => 1, 'model' => 1 }, -branchlengths => 1);

#				PAML with single omega value (model=0 designates single omega)
#				( -params => { 'runmode' => 0, 'seqtype' => 1, 'model' => 0 } );

my $treeio = Bio::TreeIO->new(-format => "nexus", -file => "$tree_file");
#read in all of the trees
my %trees = ();
my $tree = $treeio->next_tree;
$paml_exec->tree($tree, {'branchLengths' => 0 }); # initialize with this first tree

while ($tree) {
	$trees{$tree->id()} = $tree;
	$tree = $treeio->next_tree;
}

foreach my $aln (@gene_alns) {
	my $name = $aln->description();
	#print "$name\n";
	my $resultstr = $name;
 	$paml_exec->alignment($aln);
 	if (keys(%trees) != 1) {
		$paml_exec->tree($trees{$name}, {'branchLengths' => 1 });
		print "using " . $trees{$name}->id() . " for tree\n";
	}
	my %params = $paml_exec->get_parameters();
	foreach my $k (keys %params) {
		print $k . "=>" . $params{$k} . "\n";
	}
 	$paml_exec->outfile_name("$output_name"."_$name.mlc");
 	my ($rc,$parser) = $paml_exec->run();
	if ($rc == 0) {
		my $t = $paml_exec->error_string();
		print $t . "\n";
	} else {
		while( my $result = $parser->next_result ) {
			my @otus = $result->get_seqs();
			my $MLmatrix = $result->get_MLmatrix();
			print $result->model() . "\n";
		}
	}
}
