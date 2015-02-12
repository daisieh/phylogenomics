#!/usr/bin/env perl

use FindBin;
use lib "$FindBin::Bin/../../lib";
use Subfunctions;
use Bio::SeqIO;
use Bio::TreeIO;
use Bio::Align::Utilities qw(cat);
use Bio::Tools::Run::Phylo::Hyphy::BatchFile;
use File::Basename;
use Cwd 'abs_path';
use Getopt::Long;

my $usage = "perl " . basename($0);
$usage .= " gb_file fa_file tree_file output_name [num_threads]\n";

my $gb_file = shift or die $usage;
my $fa_file = shift or die $usage;
my $tree_file = shift or die $usage;
my $output_name_as_entered = shift or die $usage;
my $num_threads = shift;
unless (defined $num_threads) {
    $num_threads = 1;
}

my $output_name = abs_path( $output_name_as_entered );

my $whole_aln = make_aln_from_fasta_file ($fa_file);

my @gene_alns = @{parse_aln_into_genes($whole_aln, $gb_file, 1)};

my $treeio = Bio::TreeIO->new(-format => "nexus", -file => "$tree_file");
#read in all of the trees
my %trees = ();
my $tree = $treeio->next_tree;
my $firsttree = $tree;
my %pids = ();

while ($tree) {
	$trees{$tree->id()} = $tree;
	$tree = $treeio->next_tree;
}
open OUT_FH, ">", "$output_name.lrt";
print OUT_FH "gene\tHa=single omega\tHa=local omegas\tglobal omega\n";
my @pids;
my $last_aln = @gene_alns[(scalar @gene_alns) - 1];
foreach my $aln (@gene_alns) {
	my $name = $aln->description();
    if ($trees{$name} == undef) {
        print "skipping $name because tree is not available\n";
        next;
    }
    if (keys(%trees) != 1) { # if there's only one tree, use that tree for every alignment.
        $tree = $trees{$name};
    } else {
        $tree = $firsttree;
    }
    my $hyphy_pid = fork_hyphy($aln, $tree);
    push @pids, $hyphy_pid;
    if ("$last_aln" eq "$aln") {
        wait_for_pids(@pids);
        last;
    }
    if ((scalar @pids) >= $num_threads) {
        wait_for_pids(@pids);
        @pids = ();
    }
}

print "Done with HYPHY. Processing outputs...\n";

foreach my $aln (@gene_alns) {
	my $name = $aln->description();
	open FH, "<:crlf", "$output_name"."_$name.bfout" or next;
	my @output_fh = <FH>;
	close FH;
    print "$name...";
	my $output = join("\n", @output_fh);
	my ($omega, $p_value1, $pvalue2) = "";
	$output =~ m/Global omega calculated to be (.+?)\n/g;
	$omega = $1;
	$output =~ m/LRT for single omega across the tree: p-value = (.+?),/g;
	$p_value1 = $1;
	$output =~ m/LRT for variable omega across the tree: p-value = (.+?),/g;
	$p_value2 = $1;
	print OUT_FH "$name\t$p_value1\t$p_value2\t$omega\n";
}

print "\nResults in $output_name.lrt\n";

sub wait_for_pids {
    @pid_in = @_;
    if ((scalar @pid_in) > 1) {
        print "waiting on pids " . join (", ",@pids) . "...\n";
    }
    foreach my $item (@pids) {
        print "waiting for pid $item to return...\n";
        waitpid $item, 0;
    }
    return;
}

sub fork_hyphy {
    $aln = shift;
    $trees = shift;
    my $name = $aln->description();
    my $child_pid = fork();
    unless ($child_pid) { #child process
        my $bf_exec = Bio::Tools::Run::Phylo::Hyphy::BatchFile->new(-params => {'bf' => "ModelTest.bf", 'order' => [$aln, $firsttree, '4', 'AIC Test',  "$output_name"."_$name.aic"]});
        my $resultstr = $name;
        $bf_exec->alignment($aln);
        if ($trees{$name} == undef) {
            print "skipping $name because tree is not available\n";
            next;
        }
        if (keys(%trees) != 1) {
            $bf_exec->tree($trees{$name}, {'branchLengths' => 1 });
        }
        $bf_exec->outfile_name("$output_name"."_$name.bfout");
        my ($rc,$parser) = $bf_exec->run();
        if ($rc == 0) {
            my $t = $bf_exec->error_string();
            print $t . "\n";
        }
        open FH, "<:crlf", $bf_exec->outfile_name();
        my @output_fh = <FH>;
        close FH;

        my $output = join("\n", @output_fh);
        $output =~ m/Model String:(\d+)/g;
        my $model = $1;
        print "Model string $model chosen for $name\n";
        $bf_exec = Bio::Tools::Run::Phylo::Hyphy::BatchFile->new(-params => {'bf' => "", 'order' => ["Universal", $bf_exec->alignment, $model, $bf_exec->tree]});
        $bf_exec->alignment($aln);
        $bf_exec->tree($trees{$name}, {'branchLengths' => 1 });
        $bf_exec->outfile_name("$output_name"."_$name.bfout");
        my $bf = $bf_exec->make_batchfile_with_contents(batchfile_text());
        my ($rc,$parser) = $bf_exec->run();
        if ($rc == 0) {
            my $t = $bf_exec->error_string();
            print "There was an error: " . $t . "\n";
        }
        exit;
    } else { #parent process
        return $child_pid;
    }
}

sub batchfile_text {
    return qq{
RequireVersion ("0.9920060830");
VERBOSITY_LEVEL = -1;

ExecuteAFile (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"Utility"+DIRECTORY_SEPARATOR+"DescriptiveStatistics.bf");
ExecuteAFile (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"TemplateModels"+DIRECTORY_SEPARATOR+"chooseGeneticCode.def");

ModelMatrixDimension = 64;
for (k=0; k<64; k=k+1)
{
	if (_Genetic_Code[k] == 10)
	{
		ModelMatrixDimension = ModelMatrixDimension -1;
	}
}

ExecuteAFile (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"2RatesAnalyses"+DIRECTORY_SEPARATOR+"MG94xREV.mdl");

SetDialogPrompt     ("Choose a nucleotide alignment");
DataSet ds        = ReadDataFile (PROMPT_FOR_FILE);

DataSetFilter	  	filteredData = CreateFilter (ds,3,"","",GeneticCodeExclusions);

SKIP_MODEL_PARAMETER_LIST = 1;
done 					  = 0;

while (!done)
{
	fprintf (stdout,"\nPlease enter a 6 character model designation (e.g:010010 defines HKY85):");
	fscanf  (stdin,"String", modelDesc);
	if (Abs(modelDesc)==6)
	{
		done = 1;
	}
}
modelType 				  = 0;
ExecuteAFile (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"TemplateModels"+DIRECTORY_SEPARATOR+"MG94custom.mdl");
SKIP_MODEL_PARAMETER_LIST = 0;

ExecuteAFile 		(HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"queryTree.bf");

brNames				= BranchName (givenTree,-1);
COVARIANCE_PARAMETER 				= {};
global global_OMEGA = 1;
for (k=0; k < Columns (brNames)-1; k=k+1)
{
	ExecuteCommands ("givenTree."+brNames[k]+".nonSynRate:=givenTree."+brNames[k]+".omega*givenTree."+brNames[k]+".synRate;");
	COVARIANCE_PARAMETER["givenTree."+brNames[k]+".omega"] = 1;
}

LikelihoodFunction  theLnLik = (filteredData, givenTree);


for (k=0; k < Columns (brNames)-1; k=k+1)
{
//  set all of the branches to have omega constrained to the same global_OMEGA
	ExecuteCommands ("givenTree."+brNames[k]+".omega:=global_OMEGA;");
}

fprintf 					   (stdout, "\nFitting the global model to the data...\n");
Optimize 					   (res_global, theLnLik);
fprintf						   (stdout, theLnLik,"\n\n");
global omega_MLE = global_OMEGA;

for (k=0; k < Columns (brNames)-1; k=k+1)
{
//  set each branch to have unconstrained omega
	ExecuteCommands ("givenTree."+brNames[k]+".omega=global_OMEGA;");
}

fprintf 					   (stdout, "\nFitting the local model to the data...\n");
Optimize 					   (res_local, theLnLik);
fprintf						   (stdout, theLnLik,"\n\n");

for (k=0; k < Columns (brNames)-1; k=k+1)
{
//  set each branch to have omega = 1
	ExecuteCommands ("givenTree."+brNames[k]+".omega:=1;");
}

fprintf 					   (stdout, "\nFitting the neutral model to the data...\n");
Optimize 					   (res_neutral, theLnLik);
fprintf						   (stdout, theLnLik,"\n\n");

fprintf (stdout, "\nGlobal omega calculated to be ", omega_MLE, "\n");

LR = 2(res_global[1][0]-res_neutral[1][0]);
DF = res_global[1][1]-res_neutral[1][1];

fprintf (stdout, "\nLRT for single omega across the tree: p-value = ", 1-CChi2(LR,DF), ", LR = ", LR, ", Constraints = ", DF, "\n\n");

LR = 2(res_local[1][0]-res_global[1][0]);
DF = res_local[1][1]-res_global[1][1];

fprintf (stdout, "\nLRT for variable omega across the tree: p-value = ", 1-CChi2(LR,DF), ", LR = ", LR, ", Constraints = ", DF, "\n\n");

COVARIANCE_PRECISION = 0.95;
CovarianceMatrix (covMx, theLnLik);
//
VERBOSITY_LEVEL = 0;
//
for (k=0; k < Columns (brNames)-1; k=k+1)
{
	fprintf (stdout, "Branch :", brNames[k], "\n\tomega MLE = ", Format (covMx[k][1],6,3), "\n\t95% CI = (",Format (covMx[k][0],6,3), ",", Format (covMx[k][2],6,3), ")\n");
}

    };
}
