#!/usr/bin/perl

use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/../lib";
use Genbank;
use Subfunctions qw(parse_fasta);

my $regionfile = shift;
my $fastafile = shift;
my $featurefile = shift;
my $outfile = shift;

my ($gene_array, $gene_index_array) = parse_regionfile($regionfile);
my ($destination_gene_array, $destination_gene_index_array) = parse_featurefile($featurefile);

my ($fastahash, undef) = parse_fasta($fastafile);
my $name = "";
my $seqlen = 0;
foreach my $k (keys $fastahash) {
	# there should be only one key, so just one name.
	$name = $k;
	$seqlen = length ($fastahash->{$k});
	set_sequence($fastahash->{$k});
}

my $gene_hash = {};
foreach my $id (@$gene_index_array) {
	my $gene = shift $gene_array;
	$gene_hash->{$id} = $gene;
}

my $dest_gene_hash = {};
foreach my $id (@$destination_gene_index_array) {
	my $gene = shift $destination_gene_array;
	$dest_gene_hash->{$id} = $gene;
}

# fill in the genes from the regionfile with the info from the destination gene array
my @final_gene_array = ();
foreach my $id (@$gene_index_array) {
	my $gene = $gene_hash->{$id};
	my $dest_gene = $dest_gene_hash->{$id};
	print "gene id $id $gene\n";

	foreach my $q (keys %{$dest_gene->{"qualifiers"}}) {
		if (!(exists $gene->{"qualifiers"}->{$q})) {
			$gene->{"qualifiers"}->{$q} = $dest_gene->{"qualifiers"}->{$q}
		}
	}
	my @new_contains = ();
	$gene->{"id"} = $id;
	print Dumper($gene);
	foreach my $destcontains (@{$dest_gene->{"contains"}}) {
		my $genecontains = shift $gene->{"contains"};
		$destcontains->{"region"} = $genecontains->{"region"};
		push @new_contains, $destcontains;
	}
	$gene->{"contains"} = \@new_contains;
	push @final_gene_array, $gene;
}

# print Dumper(@final_gene_array);
open FH, ">", $outfile;

# print header
print FH ">Features\t$name\n";

open FASTA_FH, ">", "$outfile.fasta";
# start printing genes
foreach my $gene (@final_gene_array) {
	# first, print overall gene information
	my $genename = $gene->{"qualifiers"}->{"gene"};
	my $gene_id = $gene->{"id"};
	my $feat_id = 0;
	foreach my $r (@{$gene->{'region'}}) {
		$r =~ /(\d+)\.\.(\d+)/;
		print FH "$1\t$2\tgene\n";
# 		my $regseq = sequence_for_interval ($r);
# 		print FASTA_FH ">$gene_id"."_$feat_id"."_$genename"."_$featname($strand)\t$start\t$end\n$regseq\n";
	}
	foreach my $q (keys %{$gene->{'qualifiers'}}) {
		print FH "\t\t\t$q\t$gene->{qualifiers}->{$q}\n";
	}
# 	print FH "\t\t\t$q\t$gene->{qualifiers}->{$q}\n";

	# then, print each feature contained.
	foreach my $feat (@{$gene->{'contains'}}) {
		foreach my $reg (@{$feat->{"region"}}) {
			my $strand = "+";
			my ($start, $end) = split (/\.\./, $reg);
			if ($end < $start) {
				$strand = "-";
				my $oldend = $start;
				$start = $end;
				$end = $oldend;
			}
			my $regseq = sequence_for_interval ($reg);
			my $featname = $feat->{"type"};
			print FASTA_FH ">$gene_id"."_$feat_id"."_$genename"."_$featname($strand)\t$start\t$end#$regseq\n";
			$feat_id++;
		}
		print FH sequin_feature ($feat->{'region'}, $feat);
	}
}

close FH;

foreach my $gene (@final_gene_array) {
}
close FASTA_FH;

# >Features Populus trichocarpa chloroplast, complete genome.
# 1	157033	source
# 			mol_type	genomic DNA
# 			organism	Populus trichocarpa
# 77	4	gene
# 			gene	trnH
# 77	4	tRNA
# 			gene	trnH
# 			product	tRNA-His

