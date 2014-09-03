#!/usr/bin/env perl

use FindBin;
use lib "$FindBin::Bin/../..";
use Subfunctions;
use Bio::SeqIO;
use Bio::TreeIO;
use Bio::Root::Test;
use PostScript::Simple;
use Bio::Align::Utilities qw(cat);
use Bio::Tools::Run::Phylo::Hyphy::REL;
use Bio::Tools::Run::Phylo::Hyphy::BatchFile;
use Bio::Tools::Run::Phylo::Hyphy::SLAC;
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


my $treeio = Bio::TreeIO->new(-format => "nexus", -file => "$tree_file");
#read in all of the trees
my %trees = ();
my $tree = $treeio->next_tree;
my $firsttree = $tree;

while ($tree) {
	$trees{$tree->id()} = $tree;
	$tree = $treeio->next_tree;
}

foreach my $aln (@gene_alns) {
	my $name = $aln->description();
#          1 {'tempalnfile' => undef }, # aln file goes here
#          2 {'temptreefile' => undef }, # tree file goes here
#          3 {'Number of Rate Classes' => [ '4' ] },
#          4 {'Model Selection Method' => [ 'Both',
#                                         'Hierarchical Test',
#                                         'AIC Test'] },
#          5 {'Model rejection level' => '0.05' },
#          6 {'outfile' => undef },
#          7 {'aicoutfile' => undef }
	my $bf_exec = Bio::Tools::Run::Phylo::Hyphy::BatchFile->new(-params => {'bf' => "ModelTest.bf", 'order' => [$aln, $firsttree, '4', 'AIC Test',  "$output_name"."_$name.aic"]});
# 	my $bf_exec = Bio::Tools::Run::Phylo::Hyphy::BatchFile->new(-params => {'bf' => "/Users/daisie/Documents/Work/Sandbox/hyphy/examples/LRT.bf", 'order' => ["Universal", "Custom", $aln, "001001", $firsttree]});
	my $resultstr = $name;
 	$bf_exec->alignment($aln);
# 	$bf_exec->set_parameter('3', "012012");
 	if ($trees{$name} == undef) {
 		print "skipping $name because tree is not available\n";
 		next;
 	}
 	if (keys(%trees) != 1) {
		$bf_exec->tree($trees{$name}, {'branchLengths' => 1 });
		print "using " . $trees{$name}->id() . " for tree\n";
	}
 	$bf_exec->outfile_name("$output_name"."_$name.bfout");
 	my ($rc,$parser) = $bf_exec->run();
	if ($rc == 0) {
		my $t = $bf_exec->error_string();
		print ">>" . $t . "\n";
	}
	open FH, "<", $bf_exec->outfile_name();
	my @output_fh = <FH>;
	close FH;

	my $output = join("\n", @output_fh);
	$output =~ m/Model String:(\d+)/g;
	my $model = $1;
	print "running LRT on $name...\n";
	$bf_exec = Bio::Tools::Run::Phylo::Hyphy::BatchFile->new(-params => {'bf' => "/Users/daisie/Documents/Work/Sandbox/hyphy/examples/LRT.bf", 'order' => ["Universal", "Custom", $bf_exec->alignment, $model, $bf_exec->tree]});
	$bf_exec->alignment($aln);
	$bf_exec->tree($trees{$name}, {'branchLengths' => 1 });
	$bf_exec->outfile_name("$output_name"."_$name.bfout");
 	my ($rc,$parser) = $bf_exec->run();
	if ($rc == 0) {
		my $t = $bf_exec->error_string();
		print ">>" . $t . "\n";
	}

}
