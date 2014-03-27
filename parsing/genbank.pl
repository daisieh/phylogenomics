#!/usr/bin/perl

use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/..";
use Subfunctions qw(split_seq reverse_complement);

my $gbfile = shift;

open FH, "<", $gbfile;

my @gene_array = ();
my $line = "";
my $in_features = 0;
my $in_sequence = 0;
my $feat_desc_string = "";
my $sequence = "";
my $curr_gene = 0;

my $name = "";
$line = readline FH;
while (defined $line) {
	if ($line =~ /^\s+$/) {
		# if the line is blank, skip to end
	} elsif ($line =~ /^\/\//) { # we've finished the file
		print "DONE!\n";
		last;
	} elsif ($line =~ /^\S/) { # we're looking at a new section
		if ($line =~ /DEFINITION\s+(.*)$/) { # if line has the definition, keep this.
			$name = $1;
		} elsif ($line =~ /FEATURES/) { # if we hit a line that says FEATURES, we're starting the features.
			print "starting features\n";
			$in_features = 1;
		} elsif ($line =~ /ORIGIN/) { # if we hit a line that says ORIGIN, we have sequence
			# this is the end of the features: process the last feature and finish.
			print "finished features.\n";
			parse_feature_desc ($feat_desc_string);
			$in_features = 0;
			$in_sequence = 1;
			print "starting sequence for $name\n";
		}
	} elsif ($in_sequence == 1) {
		$line =~ s/\d//g;
		$line =~ s/\s//g;
		$sequence .= $line;
		# process lines in the sequence.
	} elsif ($in_features == 1) {
		# three types of lines in features:
		if ($line =~ /^\s{5}(\S.+)$/) {
			# start of a feature descriptor. Parse whatever feature and qualifiers we might've been working on and move on.
			if ($feat_desc_string ne "") {
				my $feat_desc_hash = parse_feature_desc ($feat_desc_string);
			}
			$feat_desc_string = $1;
		} elsif ($line =~ /^\s{21}(\S.+)$/) {
			$feat_desc_string .= "+$1";
		}
	}

	$line = readline FH;
}
close FH;

# print Dumper(@gene_array);

open FASTA_FH, ">", "$gbfile.fasta";
open FEATURES_FH, ">", "$gbfile.features";
my $gene_id = 0;
foreach my $gene (@gene_array) {
	if ($gene->{"type"} eq "gene") {
		my $interval_str = flatten_interval ($gene->{"region"});
		my $geneseq = sequence_for_interval ($interval_str);
		my $genename = $gene->{"features"}->{"gene"};
# 		print FASTA_FH ">$gene_id"."_$genename\n$geneseq\n";
		foreach my $feat (@{$gene->{"contains"}}) {
			my $feat_id = 0;
			foreach my $reg (@{$feat->{"region"}}) {
				my $regseq = sequence_for_interval ($reg);
				my $featname = $feat->{"type"};
				print FASTA_FH ">$gene_id"."_$feat_id"."_$genename"."_$featname\n$regseq\n";
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
# print FASTA_FH ">$name\n$sequence\n";
close FASTA_FH;
close FEATURES_FH;

open FEATURES_FH, ">", "$gbfile.sequin";
print FEATURES_FH ">Features $name\n";
foreach my $gene (@gene_array) {
	my $interval_str = flatten_interval ($gene->{"region"});
	my @gene_interval = ($interval_str);
	print FEATURES_FH sequin_feature(\@gene_interval, $gene);
	foreach my $feat (@{$gene->{"contains"}}) {
		my @regions = ();
		foreach my $r (@{$feat->{"region"}}) { push @regions, $r; }
		print FEATURES_FH sequin_feature (\@regions, $feat);
	}
}
close FEATURES_FH;

sub sequence_for_interval {
	my $interval_str = shift;

	$interval_str =~ /(\d+)\.\.(\d+)/;
	my $start = $1;
	my $end = $2;
	my $geneseq = "";
	if ($start < $end) {
		(undef, $geneseq, undef) = split_seq ($sequence, $start, $end);
	} else {
		(undef, $geneseq, undef) = split_seq ($sequence, $end, $start);
# 		$geneseq = reverse_complement ($geneseq);
	}
	return $geneseq;
}

sub sequin_feature {
	my $regions = shift;
	my $feature = shift;
	my $result = "";

	my $first_int = shift @$regions;
	$first_int =~ /(\d+)\.\.(\d+)/;
	$result = "$1\t$2\t$feature->{type}\n";
	foreach my $int (@$regions) {
		$int =~ /(\d+)\.\.(\d+)/;
		$result .= "$1\t$2\n";
	}
	foreach my $key (keys %{$feature->{"features"}}) {
		$result .= "\t\t\t$key\t$feature->{features}->{$key}\n";
	}
	return $result;
}

sub stringify_feature {
	my $regions = shift;
	my $feature = shift;
	my $result = "";

	my @features = ();
	foreach my $key (keys %{$feature->{"features"}}) {
		my $feat = "$key=$feature->{features}->{$key}";
		push @features, $feat;
	}
	$result = "$feature->{type}\t".join(",",@$regions)."\t".join("#",@features)."\n";
	return $result;
}


sub flatten_interval {
	my $int_array = shift;

	# calculate the largest extent of the main interval.
	my @locs = ();
	foreach my $r (@$int_array) {
		if ($r =~ /(\d+)\.\.(\d+)/) {
			push @locs, $1;
			push @locs, $2;
		}
	}
	# sort all the listed locations, make the start be the smallest possible val and the end be the largest possible val.
	@locs = sort {$a <=> $b} @locs;
	my $main_start = shift @locs;
	my $main_end = pop @locs;

	# was the original interval complemented?
	my $r = @$int_array[0];#join (",", @$int_array);
	my ($loc1, $loc2) = split (/\.\./, $r);
	if ($loc1 < $loc2) {
		return "$main_start..$main_end";
	} else {
		return "$main_end..$main_start";
	}

}

sub parse_feature_desc {
	my $feat_desc_string = shift;
	my $feat_desc_hash = {};

	my ($feature_desc, @feature_quals) = split (/\+\//, $feat_desc_string);
	#parse feature_desc into name and interval.
	$feature_desc =~ s/\+//g;
	$feature_desc =~ /^(.{16})(.+)$/;
	my $type = $1;
	my $region = $2;

	$type =~ s/ *$//; # remove trailing blanks.

	$feat_desc_hash->{"type"} = $type;
	$feat_desc_hash->{"region"} = parse_interval($region);
	$feat_desc_hash->{"features"} = parse_features(\@feature_quals);
	if ($feat_desc_hash->{"type"} eq "gene") {
		$curr_gene = {};
		push @gene_array, $curr_gene;
		$curr_gene->{"contains"} = ();
		$curr_gene->{"region"} = $feat_desc_hash->{"region"};
		$curr_gene->{"features"} = $feat_desc_hash->{"features"};
		$curr_gene->{"type"} = "gene";
	} else {
		if (within_interval($curr_gene->{"region"}, $feat_desc_hash->{"region"})) {
			# this feat_desc_hash belongs to the current gene.
			push @{$curr_gene->{"contains"}}, $feat_desc_hash;
		} else {
			$curr_gene = $feat_desc_hash;
			push @gene_array, $feat_desc_hash;
		}
	}

	return $feat_desc_hash;
}

sub parse_interval {
	my $intervalstr = shift;
	my @regions = ();
	if ($intervalstr =~ /^complement\s*\((.+)\)/) {
		# this is a complementary strand feature.
		my $subregions = parse_interval($1);
		foreach my $subreg (@$subregions) {
			if ($subreg =~ /(\d+)\.\.(\d+)/) {
				push @regions, "$2..$1";
			}
		}
	} elsif ($intervalstr =~ /^join\s*\((.+)\)$/) {
		# this is a series of intervals
		my @subintervals = split(/,/, $1);
		foreach my $subint (@subintervals) {
			my $subregions = parse_interval($subint);
			push @regions, @$subregions;
		}
	} elsif ($intervalstr =~ /(\d+)\.\.(\d+)/) {
		push @regions, "$intervalstr";
	}
	return \@regions;
}

sub parse_features {
	my $features_ref = shift;

	my @features = @$features_ref;
	my $feature_hash = {};
	while (@features > 0) {
		my $f = shift @features;
		$f =~ s/\+/ /g;
		if ($f =~ /(.+)=(.+)/) {
			my $key = $1;
			my $val = $2;
			if ($val =~ /"(.*)"/) {
				$val = $1;
			}
			if ($key eq "translation") {
				next;
				$val =~ s/ //g;
			}
			$feature_hash->{$key} = $val;
		} elsif ($f =~ /(.+)/) {
			my $key = $1;
			$feature_hash->{$key} = "";
		} else {
			print "haven't dealt with this: $f\n";
		}
	}
	return $feature_hash;
}

sub within_interval {
	my $main_interval = shift;
	my $test_interval = shift;

	if (($main_interval eq "") || ($test_interval eq "")) {
		# if the interval we're testing for is blank, return 0.
		return 0;
	}

	if ((ref $test_interval) !~ /ARRAY/) {
		$test_interval = parse_interval($test_interval);
	}

	if ((ref $main_interval) !~ /ARRAY/) {
		$main_interval = parse_interval($main_interval);
	}

	# calculate the largest extent of the main interval.
	my @locs = ();
	foreach my $r (@$main_interval) {
		if ($r =~ /(\d+)\.\.(\d+)/) {
			push @locs, $1;
			push @locs, $2;
		}
	}
	# sort all the listed locations, make the start be the smallest possible val and the end be the largest possible val.
	@locs = sort {$a <=> $b} @locs;
	my $main_start = shift @locs;
	my $main_end = pop @locs;

	# do the same for the tested intervals.
	@locs = ();
	foreach my $r (@$test_interval) {
		if ($r =~ /(\d+)\.\.(\d+)/) {
			push @locs, $1;
			push @locs, $2;
		}
	}
	# sort all the listed locations, make the start be the smallest possible val and the end be the largest possible val.
	@locs = sort {$a <=> $b} @locs;
	my $test_start = shift @locs;
	my $test_end = pop @locs;
	if (($test_start >= $main_start) && ($test_end <= $main_end)) {
		return 1;
	}
	return 0;
}

