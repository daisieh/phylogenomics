use Bio::SeqIO;
use Bio::TreeIO;
use Bio::Align::Utilities qw(cat);
use Bio::Tools::Run::Phylo::Hyphy::BatchFile;
use Bio::Tools::Run::Phylo::Hyphy::SLAC;

my $alignio = Bio::AlignIO->new(-format => 'phylip',
                     -file   => '/Users/daisie/Desktop/rbcL.phy');

my $treeio = Bio::TreeIO->new(-format => 'newick',
                     -file   => '/Users/daisie/Desktop/newick.tre');

my $aln = $alignio->next_aln;
my $tree = $treeio->next_tree;
my ($rc,$results);
print ref ($aln) . "\n";

# my $slac = Bio::Tools::Run::Phylo::Hyphy::SLAC->new();
# $slac->alignment($aln);
# $slac->tree($tree);
# $slac->save_tempfiles(1);
# ($rc,$results) = $slac->run();

my $bf_exec = Bio::Tools::Run::Phylo::Hyphy::BatchFile->new(-params => {'bf' => "ModelTest.bf", 'order' => [$aln, $tree, '4', 'AIC Test', ""]});
$bf_exec->alignment($aln);
$bf_exec->tree($tree);
$bf_exec->save_tempfiles(1);
print $bf_exec->tempdir()."\n";
$bf_exec->set_parameter(5, $bf_exec->tempdir() . "/output");
# $bf_exec->save_tempfiles(1);
($rc,$results) = $bf_exec->run();
