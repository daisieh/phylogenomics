require "subfuncs.pl";
use Bio::SeqIO;
use Bio::TreeIO;
use Bio::Align::Utilities qw(cat);
use Bio::Tools::Run::Phylo::PAML::Codeml;
use File::Basename;
use Getopt::Long;
use Pod::Usage;

my $usage = "perl " . basename($0);
$usage .= " gb_file fa_file tree_file output_name analysis\n";

my ($gb_file, $fa_file, $tree_file, $output_name, $analysis) = 0;
GetOptions ('genbank|gb=s' => \$gb_file,
            'fasta=s' => \$fa_file,
            'treefile=s' => \$tree_file,
            'output=s' => \$output_name,
            'analysis|model=i' => \$analysis) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

unless ($gb_file && $fa_file && $tree_file && $output_name) {
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

my $paml_exec;

if ($analysis == 0) {
#   PAML with single omega value (model=0 designates single omega)
    $paml_exec =	Bio::Tools::Run::Phylo::PAML::Codeml->new(-params => { 'runmode' => 0, 'seqtype' => 1, 'model' => 0 });
} elsif ($analysis == 1) {
#	PAML with fixed omega value (model=0 designates single omega)
    $paml_exec =	Bio::Tools::Run::Phylo::PAML::Codeml->new(-params => { 'runmode' => 0, 'seqtype' => 1, 'model' => 0, 'fix_omega' => 1, 'omega' => 1 });
} elsif ($analysis == 2) {
#	PAML with free omegas (model=1)
    $paml_exec =	Bio::Tools::Run::Phylo::PAML::Codeml->new(-params => { 'runmode' => 0, 'seqtype' => 1, 'model' => 1, 'fix_blength' => 1 });
#				( -params => { 'runmode' => 0, 'seqtype' => 1, 'model' => 1 }, -branchlengths => 1);
} else { pod2usage(1); }

$paml_exec->save_tempfiles(1);
print $paml_exec->tempdir() . "\n";

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
 	$paml_exec->outfile_name("$output_name"."_$name.mlc");
# 	foreach my $k (keys %params) {
# 		print $k . "=>" . $params{$k} . "\n";
# 	}
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
                    1:  PAML with free omegas (corresponds to model=1)
                    2:  PAML with fixed omega = 1 (neutral)

=head1 DESCRIPTION

=cut
