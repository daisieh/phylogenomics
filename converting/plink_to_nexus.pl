#!/usr/bin/env perl

use strict;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use File::Temp qw(tempfile);
use File::Basename qw(fileparse);
use FindBin;
use lib "$FindBin::Bin/..";
use lib "$FindBin::Bin/../parsing";
use Subfunctions qw(debug set_debug get_iupac_code);
use Nexus qw(write_nexus_character_block write_nexus_trees_block write_nexus_taxa_block);

my $help = 0;
my $outfile = "";
my $inputmap = "";
my $inputped = "";

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

GetOptions ('map=s' => \$inputmap,
			'ped=s' => \$inputped,
			'output=s' => \$outfile,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help){
    pod2usage(-verbose => 1);
}

if ($outfile eq "") {
	$inputped =~ /(.+)\.ped/;
	$outfile = "$1.nex";
}

if ($outfile !~ /\.nex$/) {
	$outfile = "$outfile.nex";
}

if (($inputmap eq "") && ($inputped eq "")) {
	pod2usage(-msg => "Both an input .ped and an input .map file are required.", -exitval => 2);
}

if ($inputmap !~ /\.map$/) {
	pod2usage(-msg => "File $inputmap is not a .map file.", -exitval => 2);
}

if ($inputped !~ /\.ped$/) {
	pod2usage(-msg => "File $inputped is not a .ped file.", -exitval => 2);
}

print "processing .map file...\n";
my $snps = ();
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
	}
	push @$snps, $snp_hash;
}
close MAP_FH;

my $total_snp_count = @$snps;
my $individuals = {};
my $indiv_array = ();

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
print "processing .ped file...\n";
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
		push @$indiv_array, $indiv_hash->{"individual_id"};
		$individuals->{$indiv_hash->{"individual_id"}} = $indiv_hash;
	}
}
close PED_FH;
print "mapping genotypes...\n";
foreach my $indiv_id (@$indiv_array) {
	print "looking at $indiv_id\n";
	my $indiv = $individuals->{$indiv_id};
	my $alleles = "$indiv->{alleles}";
	$alleles =~ s/\s+//g;
	my $genotype = "";
	while ($alleles ne "") {
		$alleles =~ /(..)(.*)/;
		my $snppair = $1;
		$alleles = $2;
 		if ($snppair eq "00") {
 			$genotype .= "-";
 		} else {
			$genotype .= get_iupac_code($snppair);
		}
	}
	$indiv->{"genotype"} = $genotype;
	print "$indiv_id\t$genotype\n";
	if ((length $genotype) != $total_snp_count) {
		print "$indiv->{individual_id} has an incorrect number of snps specified in the ped file: ". (length $genotype) . " instead of $total_snp_count.\n";
		exit;
	}
}

# print "individuals " .Dumper ($individuals);

# write out as nexus:
print "writing output to $outfile\n";
my $nexusstring = "#NEXUS\n\n";

$nexusstring .= write_nexus_taxa_block($indiv_array);
my $taxahash = {};
foreach my $indiv_id (@$indiv_array) {
	$taxahash->{$indiv_id} = $individuals->{$indiv_id}->{"genotype"};
}
$nexusstring .= write_nexus_character_block($taxahash, $indiv_array);

open OUT_FH, ">", $outfile;
print OUT_FH $nexusstring;
close OUT_FH;

