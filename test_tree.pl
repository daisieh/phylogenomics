use Bio::SeqIO;
use Bio::TreeIO;
use Bio::Align::Utilities qw(cat);
require "subfuncs.pl";

my $treefile = shift;
my $treeio = Bio::TreeIO->new(-format => 'nexus',
                     -file   => $treefile);

my $tree = $treeio->next_tree;

my @nodes = $tree->get_nodes();
foreach my $node (@nodes) {
	$node->remove_all_tags();
	print $node->id() . ": " . $node->branch_length() . "\n";
}

print $tree->as_text('newick') . "\n";
