#!/usr/bin/env perl

use FindBin;
use lib "$FindBin::Bin/../../lib";
use Subfunctions;
use Bio::SeqIO;
use Bio::TreeIO;
use Bio::Align::Utilities qw(cat);
use Bio::Tools::Run::Phylo::PAML::Codeml;
use File::Basename;
use Getopt::Long;
use Pod::Usage;

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my ($gb_file, $fa_file, $tree_file, $output_name, $analysis, $cleandata) = 0;
GetOptions ('genbank|gb:s' => \$gb_file,
            'fasta=s' => \$fa_file,
            'treefile=s' => \$tree_file,
            'output:s' => \$output_name,
            'cleandata' => \$cleandata,
            'analysis|model=i' => \$analysis) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

unless ($fa_file && $tree_file) {
    my $msg = qq{Error: an option was mis-specified:
    genbank = $gb_file
    fasta = $fa_file
    treefile = $tree_file
    output = $output_name
    analysis = $analysis
};
    pod2usage(-msg => $msg, -exitval => 2);
}


my $whole_aln = make_aln_from_fasta_file ($fa_file);
my @gene_alns;
if ($gb_file) {
	print "parsing whole fasta file into genes...\n";
	@gene_alns = @{parse_aln_into_genes($whole_aln, $gb_file, 1)};
} else {
	@gene_alns = ();
	my $temp_aln = make_aln_from_fasta_file($fa_file);
	$fa_file =~ /(.*?)\.fa.*/;
	$temp_aln->description("$1");
	push @gene_alns, $temp_aln;
}

my $paml_exec;

my %paramset = ();
$paramset{'runmode'} = 0;
$paramset{'seqtype'} = 1;

if ($cleandata) {
	$paramset{'cleandata'} = 1;
}

my $analysis_name = "";

if ($analysis == 0) {
#   PAML with single omega value
#     $paml_exec =	Bio::Tools::Run::Phylo::PAML::Codeml->new(-params => { 'runmode' => 0, 'seqtype' => 1, 'model' => 0 });
	$analysis_name = "single_omega";
	$paramset{'model'} = 0;
} elsif ($analysis == 1) {
#	PAML with fixed omega value of 1
#     $paml_exec =	Bio::Tools::Run::Phylo::PAML::Codeml->new(-params => { 'runmode' => 0, 'seqtype' => 1, 'model' => 0, 'fix_omega' => 1, 'omega' => 1 });
	$analysis_name = "fixed_omega";
	$paramset{'model'} = 0;
	$paramset{'fix_omega'} = 1;
	$paramset{'omega'} = 1;
} elsif ($analysis == 2) {
#	PAML with free omegas (model=1)
#     $paml_exec =	Bio::Tools::Run::Phylo::PAML::Codeml->new(-params => { 'runmode' => 0, 'seqtype' => 1, 'model' => 1, 'fix_blength' => 1 });
	$analysis_name = "free_omega";
	$paramset{'model'} = 1;
	$paramset{'fix_blength'} = 1;
} elsif ($analysis == 3) {
#	PAML with branch-site model A (model=2 NSSites=2)
#     $paml_exec =	Bio::Tools::Run::Phylo::PAML::Codeml->new(-params => { 'runmode' => 0, 'seqtype' => 1, 'model' => 2, 'NSsites' => 2});
	$analysis_name = "branch_site_A";
	$paramset{'model'} = 2;
	$paramset{'NSsites'} = 2;
} elsif ($analysis == 4) {
#	PAML with null model A1 (model = 2, NSsites = 2, fix_omega = 1, omega = 1)
#     $paml_exec =	Bio::Tools::Run::Phylo::PAML::Codeml->new(-params => { 'runmode' => 0, 'seqtype' => 1, 'model' => 2, 'NSsites' => 2, 'fix_omega' => 1, 'omega' => 1 });
	$analysis_name = "branch_site_A1";
	$paramset{'model'} = 2;
	$paramset{'fix_omega'} = 1;
	$paramset{'omega'} = 1;
} else { pod2usage(1); }

if ($output_name == 0) {
	$output_name = "$analysis_name";
}

$paml_exec =	Bio::Tools::Run::Phylo::PAML::Codeml->new(-params => \%paramset);
#     $paml_exec =	Bio::Tools::Run::Phylo::PAML::Codeml->new(-params => { 'runmode' => 0, 'seqtype' => 1, 'model' => 1, 'fix_blength' => 1 , 'cleandata' => 1});

$paml_exec->save_tempfiles(1);
print $paml_exec->tempdir() . "\n";

my $treeio = Bio::TreeIO->new(-format => "nexus", -file => "$tree_file");
#read in all of the trees
my %trees = ();
my $tree = $treeio->next_tree;
my $firsttree = $tree;

while ($tree) {
	print "tree\n";
	$trees{$tree->id()} = $tree;
	$tree = $treeio->next_tree;
}

foreach my $aln (@gene_alns) {
	my $name = $aln->description();
	print "$name\n";
	my $resultstr = $name;
 	$paml_exec->alignment($aln);
 	if (keys(%trees) != 1) {
        if ($trees{$name} == undef) {
            print "skipping $name because tree is not available\n";
            next;
        }
		$paml_exec->tree($trees{$name}, {'branchLengths' => 1 });
		print "using " . $trees{$name}->id() . " for tree\n";
	} else {
        $tree = $firsttree;
        $paml_exec->tree($tree, {'branchLengths' => 0 }); # initialize with this first tree
    }
	my %params = $paml_exec->get_parameters();
 	$paml_exec->outfile_name("$name.mlc");
 	my ($rc,$parser) = $paml_exec->run();
	if ($rc == 0) {
		my $t = $paml_exec->error_string();
		print "Error: " . $t . "\n";
	} else {
		while( my $result = $parser->next_result ) {
			my @otus = $result->get_seqs();
			my $MLmatrix = $result->get_MLmatrix();
			print $result->model() . "\n";
		}
		system ("mkdir $output_name");
		system ("cp " . $paml_exec->tempdir() . "/$name.mlc" . " $output_name");
	}
}

__END__

=head1 NAME

tree_omega

=head1 SYNOPSIS

tree_omega [options]

=head1 OPTIONS

    -genbank|gb     genbank file listing the genes to analyze
    -fasta          fasta file of whole aligned sequences
    -treefile       treefile for use by PAML
    -output         output file name
    -model|analysis CodeML analysis to perform:
                    0:  PAML with single omega
                    1:  PAML with fixed omega = 1 (neutral)
                    2:  PAML with free omegas (corresponds to model=1)
                    3:	PAML with branch-site model A (model=2 NSSites=2)
                    4:	PAML with null model A1 (model = 2, NSsites = 2, fix_omega = 1, omega = 1)



=head1 DESCRIPTION

=cut
