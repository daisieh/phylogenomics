#!/usr/bin/env perl

use FindBin;
use lib "$FindBin::Bin/../../lib";
use Subfunctions;
use Bio::SeqIO;
use Bio::TreeIO;
use Bio::Root::Test;
use PostScript::Simple;
use Bio::Align::Utilities qw(cat);
use Bio::Tools::Run::Phylo::PAML::Codeml;

my $usage  = "process_tree_omega.pl gb_file result_dir result_prefix\n";

my $gb_file = shift or die $usage;
my $result_dir = shift or die $usage;
my $result_prefix = shift;

my $output_name = "$result_dir/";
if ($result_prefix) {
	$output_name .= "$result_prefix"."_";
}

my @genes;

if ($gb_file =~ /\.gb/) {
	my $seqio_object = Bio::SeqIO->new(-file => $gb_file);
	my $seq_object = $seqio_object->next_seq;

	while ($seq_object) {
		for my $feat_object ($seq_object->get_SeqFeatures) {
			if ($feat_object->primary_tag eq "CDS") {
				my $name = main_name_for_gb_feature($feat_object);
				if ($name eq "") { next; }
				my @locations = $feat_object->location->each_Location;
				my $strand = 0;
				push @genes, $name;
			}
		}
		$seq_object = $seqio_object->next_seq;
	}
} else {
	open FH, "<", $gb_file;
	my $line = readline FH;
	while ($line) {
		chomp $line;
		push @genes, $line;
		$line = readline FH;
	}
	close FH;
}

my $outfilename = "$output_name" . "summary.txt";
open my $summaryfile, ">$outfilename";
print $summaryfile "gene\tnumparams\tlog-likelihood\n";

my $treestr;
open my $temptreefile, ">", \$treestr;
my $treeio = Bio::TreeIO->new(-format =>"newick", -fh => $temptreefile);

my $omega_tree = 0;
my %branch_omegas = ();
my @actual_genes;
foreach my $name (@genes) {
	my $resultstr = $name;
	my $filename = "$output_name"."$name.mlc";
	my $x = (-e $filename); # check to see if the file exists; if not, skip to next.
	if ($x != 1) {
		print "missing $filename\n";
		next;
	}
	print "processing $filename...\n";
	my $parser = Bio::Tools::Phylo::PAML->new (-file => $filename);
	my $result;
	eval {$result = $parser->next_result} or next;
	if ($result) {
		my @otus = $result->get_seqs();
		while (	my $result_tree = $result->next_tree() ) {
			my $likelihood = $result_tree->score();
			my $numparams = $result_tree->id();
			$treeio->write_tree($result_tree);
			$numparams =~ s/num_param://;
			print $summaryfile "$name\t$numparams\t$likelihood\n";
			for my $x ($result_tree->get_nodes()) {
				if ($x->ancestor()) {
					if ((exists $branch_omegas{$x->id()}) == 0) {
						$branch_omegas{$x->id()} = $x;
					}
					if ($x->has_tag("dN/dS")) {
						my @tags = $x->get_tag_values("dN/dS");
						$branch_omegas{$x->id()}->add_tag_value($name, $tags[0]);
						@tags = $x->get_tag_values("t");
						$branch_omegas{$x->id()}->add_tag_value("$name-t", $tags[0]);
					}
					my $newname;
					$newname = $x->ancestor()->id() . "-" . $x->id();
					$branch_omegas{$x->id()}->add_tag_value("label", $newname);
				}
			}
		}
	}
	push @actual_genes, $name;
}
close $summaryfile;

#print output
$outfilename = "$output_name" . "omegas.txt";
open my $outfile, ">$outfilename";
$outfilename = "$output_name" . "lengths.txt";
open my $outfile2, ">$outfilename";

print $outfile "node\t";
print $outfile2 "node\t";

foreach my $name (@actual_genes) {
	print $outfile "$name\t";
	print $outfile2 "$name\t";
}
print $outfile "\n";
print $outfile2 "\n";

for my $k (keys %branch_omegas) {
	my $i = $branch_omegas{$k};
	my @label = $i->get_tag_values("label");
	print $outfile "$label[0]\t";
	print $outfile2 "$label[0]\t";
	foreach my $name (@actual_genes) {
		my @values = $i->get_tag_values($name);
		print $outfile "$values[0]\t";
		@values = $i->get_tag_values("$name-t");
		print $outfile2 "$values[0]\t";
	}
	print $outfile "\n";
	print $outfile2 "\n";
}

close $outfile;
close $outfile2;

#rewrite tree file into nexus format:
$outfilename = "$output_name" . "trees.tre";
open my $outtreefile, ">$outfilename";

my @trees = split (/;/,$treestr);
print $outtreefile "#NEXUS\n\nbegin trees;\n";
for (my $i=0; $i < scalar @actual_genes; $i++) {
	#strip out the internal node labels, as this confuses the parser
	$trees[$i] =~ s/\)\d+:/\):/g;
	$trees[$i] =~ s/\)\d+$/\)/;
	print $outtreefile "tree $actual_genes[$i] = $trees[$i];\n";
}
print $outtreefile "end;\n";
close $outtreefile;
