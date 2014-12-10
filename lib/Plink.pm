package Plink;
use strict;
use Data::Dumper;
use Subfunctions qw(debug set_debug consensus_str);


BEGIN {
	require Exporter;
	# set the version for version checking
	our $VERSION     = 1.00;
	# Inherit from Exporter to export functions and variables
	our @ISA         = qw(Exporter);
	# Functions and variables which are exported by default
	our @EXPORT      = qw();
	# Functions and variables which can be optionally exported
	our @EXPORT_OK   = qw(parse_plink parse_ped parse_map);
}

=head1

Plink data structure:
These come from the ped file:
	$plink_hash->{"individuals"}: a hash of hashes, keyed by individual_id
		->{"family_id"}: family ID
		->{"individual_id"}: individual ID
		->{"paternal_id"}: paternal ID (0 for no information)
		->{"maternal_id"}: maternal ID (0 for no information)
		->{"sex"}: sex (1 if male, 2 if female, 0 if unknown)
		->{"phenotype"}: phenotype (a coded value for association genetics)
		->{"paternal"}: the paternal haplotype
		->{"maternal"}: the maternal haplotype
		->{"genotype"}: the consensus genotype
	$plink_hash->{"names"}: an array of the names of the individuals
These come from the map file:
	$plink_hash->{"snps"}: a hash of hashes, keyed by snp name
		->{"chromosome"}: chromosome (scaffold)
		->{"name"}: snp name
		->{"genetic_map_position"}: 0 (position in a genetic map)
		->{"base_pair"}: position in the physical map (base pair)
	$plink_hash->{"snp_names"}: an array of the names of the snps
	$plink_hash->{"chromosome"}: a hash of arrays, keyed by chromosome name
		each array has the base pair location of the snps in order
=cut

sub parse_plink {
	my $inputped = shift;
	my $inputmap = shift;
	my $plink_hash = shift;

	unless (defined $plink_hash) {
		$plink_hash = {};
	}
	parse_ped ($inputped, $plink_hash);
	parse_map ($inputmap, $plink_hash);
	return $plink_hash;
}

sub parse_ped {
	my $inputped = shift;
	my $plink_hash = shift;

	unless (defined $plink_hash) {
		$plink_hash = {};
	}

	$plink_hash->{"individuals"} = {};
	$plink_hash->{"names"} = ();

	open PED_FH, "<", $inputped;
	# ped file:
	# 1 ALAA20-2_DNA55 0 0 2 2 C C
	# col 1: family ID
	# col 2: individual ID
	# col 3: paternal ID (0 for no information)
	# col 4: maternal ID (0 for no information)
	# col 5: sex (1 if male, 2 if female, 0 if unknown)
	# col 6: phenotype (a coded value for association genetics)
	# next cols are pairs of allelic values corresponding to the mapped snps, 0 for missing data
	foreach my $line (<PED_FH>) {
		if ($line =~ /^(.+?)\s+(.+?)\s+(.+?)\s+(.+?)\s+(.+?)\s+(.+?)\s+(.*)$/) {
			my $indiv_hash = {};
			$indiv_hash->{"family_id"} = "$1";
			$indiv_hash->{"individual_id"} = "$2";
			$indiv_hash->{"paternal_id"} = "$3";
			$indiv_hash->{"maternal_id"} = "$4";
			$indiv_hash->{"sex"} = $5;
			$indiv_hash->{"phenotype"} = $6;
			$indiv_hash->{"alleles"} = "$7";
			push @{$plink_hash->{"names"}}, $indiv_hash->{"individual_id"};
			$plink_hash->{"individuals"}->{$indiv_hash->{"individual_id"}} = $indiv_hash;
		}
	}
	close PED_FH;

	foreach my $indiv_id (@{$plink_hash->{"names"}}) {
		my $indiv = $plink_hash->{"individuals"}->{$indiv_id};
		my $alleles = delete $indiv->{alleles};
		$alleles =~ s/\s+//g;
		$alleles =~ s/0/-/g; # replace any 0's with -'s.
		$indiv->{"paternal"} = ($alleles =~ s/([A-Za-z\-])[A-Za-z\-]/$1\./gr);
		$indiv->{"paternal"} =~ s/\.//g;
		$indiv->{"maternal"} = ($alleles =~ s/[A-Za-z\-]([A-Za-z\-])/\.$1/gr);
		$indiv->{"maternal"} =~ s/\.//g;
		my @seqarray = ($indiv->{"paternal"}, $indiv->{"maternal"});
		my $genotype = consensus_str(\@seqarray);
		$indiv->{"genotype"} = $genotype;
	}

	return $plink_hash;
}

sub parse_map {
	my $inputmap = shift;
	my $plink_hash = shift;

	unless (defined $plink_hash) {
		$plink_hash = {};
	}

	$plink_hash->{"snps"} = {};
	open MAP_FH, "<", $inputmap;
	# map file:
	# col 1: chromosome (scaffold)
	# col 2: snp name
	# col 3: 0 (position in a genetic map)
	# col 4: position in the physical map (base pair)
	foreach my $line (<MAP_FH>) {
		my $snp_hash = {};
		if ($line =~ /^(.+?)\s+(.+?)\s+(.+?)\s+(.+?)$/) {
			$snp_hash->{"chromosome"} = "$1";
			$snp_hash->{"name"} = "$2";
			$snp_hash->{"genetic_map_position"} = $3;
			$snp_hash->{"base_pair"} = $4;
			$plink_hash->{"snps"}->{$snp_hash->{"name"}} = $snp_hash;
			if (!(exists $plink_hash->{"chromosomes"}->{$snp_hash->{"chromosome"}})) {
				$plink_hash->{"chromosomes"}->{$snp_hash->{"chromosome"}} = ();
			}
			my $snp_pair = {};
# 			$snp_pair->{"name"} = $snp_hash->{"name"};
			$snp_pair->{"base_pair"} = $snp_hash->{"base_pair"};
			push @{$plink_hash->{"chromosomes"}->{$snp_hash->{"chromosome"}}}, $snp_hash->{"base_pair"};
			push @{$plink_hash->{"snp_names"}}, $snp_hash->{"name"};
		}
	}
	close MAP_FH;
	foreach my $chr (keys %{$plink_hash->{"chromosomes"}}) {
		my @list = sort {$a <=> $b} @{$plink_hash->{"chromosomes"}->{$chr}};
		$plink_hash->{"chromosomes"}->{$chr} = \@list;
	}
	return $plink_hash;
}


# must return 1 for the file overall.
1;
