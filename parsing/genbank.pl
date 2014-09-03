#!/usr/bin/perl

use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/";
use Genbank;

my $gbfile = shift;

my $gene_array = parse_genbank($gbfile);


open FASTA_FH, ">", "$gbfile.fasta";
open FEATURES_FH, ">", "$gbfile.features";
my $gene_id = 0;
foreach my $gene (@$gene_array) {
	if ($gene->{"type"} eq "gene") {
		my $interval_str = flatten_interval ($gene->{"region"});
		my $geneseq = sequence_for_interval ($interval_str);
		my $genename = $gene->{"qualifiers"}->{"gene"};
		foreach my $feat (@{$gene->{"contains"}}) {
			my $feat_id = 0;
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
				print FASTA_FH ">$gene_id"."_$feat_id"."_$genename"."_$featname($strand)\t$start\t$end\n$regseq\n";
				$feat_id++;
			}
		}
		my $interval_str = flatten_interval ($gene->{"region"});
		my @gene_interval = ($interval_str);
		print FEATURES_FH "$gene_id\t" . stringify_feature(\@gene_interval, $gene);
		foreach my $feat (@{$gene->{"contains"}}) {
			print FEATURES_FH "$gene_id\t" . stringify_feature ($feat->{"region"}, $feat);
		}
		$gene_id++;
	}
}
close FASTA_FH;
close FEATURES_FH;

# open FEATURES_FH, ">", "$gbfile.tbl";
# print FEATURES_FH ">Features $name\n";
# foreach my $gene (@$gene_array) {
# 	my $interval_str = flatten_interval ($gene->{"region"});
# 	my @gene_interval = ($interval_str);
# 	print FEATURES_FH sequin_feature(\@gene_interval, $gene);
# 	foreach my $feat (@{$gene->{"contains"}}) {
# 		my @regions = ();
# 		foreach my $r (@{$feat->{"region"}}) { push @regions, $r; }
# 		print FEATURES_FH sequin_feature (\@regions, $feat);
# 	}
# }
# close FEATURES_FH;

open DUMP_FH, ">", "$gbfile.dump";
print DUMP_FH Dumper($gene_array);
close DUMP_FH;
