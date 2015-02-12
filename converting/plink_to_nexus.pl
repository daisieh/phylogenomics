#!/usr/bin/env perl

use strict;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use File::Temp qw(tempfile);
use File::Basename qw(fileparse);
use FindBin;
use lib "$FindBin::Bin/../lib";
use Subfunctions qw(debug set_debug get_iupac_code consensus_str);
use Nexus qw(write_nexus_character_block write_nexus_trees_block write_nexus_taxa_block write_nexus_sets_block);
use Plink qw(parse_plink);
use GFF qw(parse_gff_file parse_gff_block read_gff_block);

my $help = 0;
my $outfile = "";
my $inputmap = "";
my $inputped = "";
my $inputname = "";
my $gfffile = "";

if (@ARGV == 0) {
    pod2usage(-verbose => 2);
}

GetOptions ('map=s' => \$inputmap,
			'ped=s' => \$inputped,
			'input=s' => \$inputname,
			'output=s' => \$outfile,
			'gff=s' => \$gfffile,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help){
    pod2usage(-verbose => 2);
}

if (($inputmap eq "") && ($inputped eq "")) {
	if ($inputname eq "") {
		pod2usage(-msg => "Both an input .ped and an input .map file are required.", -exitval => 2);
	} else {
		$inputmap = "$inputname.map";
		$inputped = "$inputname.ped";
	}
}

if ($inputmap !~ /\.map$/) {
	pod2usage(-msg => "File $inputmap is not a .map file.", -exitval => 2);
}

if ($inputped !~ /\.ped$/) {
	pod2usage(-msg => "File $inputped is not a .ped file.", -exitval => 2);
}

# unless (-e $inputped) {
# 	pod2usage(-msg => "File $inputped does not exist.", -exitval => 2);
# }
#
unless (-e $inputmap) {
	pod2usage(-msg => "File $inputmap does not exist.", -exitval => 2);
}

if ($outfile eq "") {
	$inputped =~ /(.+)\.ped/;
	$outfile = "$1.nex";
}

if ($outfile !~ /\.nex$/) {
	$outfile = "$outfile.nex";
}

my $plink_hash = parse_plink($inputped, $inputmap);
# print Dumper($plink_hash);
my $gene_locations = {};
if ($gfffile ne "") {
	open my $gff_fh, "<:crlf", $gfffile;
	my $gff_array = parse_gff_file ($gfffile);
	print Dumper($gff_array);
	foreach my $gff_hash (@$gff_array) {
		if (!(exists $gene_locations->{$gff_hash->{"seqid"}})) {
			$gene_locations->{$gff_hash->{"seqid"}} = ();
		}
		my $gene_hash = {};
		$gene_hash->{"name"} = $gff_hash->{"Name"};
		$gene_hash->{"start"} = $gff_hash->{"start"};
		$gene_hash->{"end"} = $gff_hash->{"end"};
		push @{$gene_locations->{ $gff_hash->{"seqid"}}}, $gene_hash;
	}
}
# print Dumper($gene_locations);
my $genes = 0;
foreach my $chr (keys %$gene_locations) {
	$genes += @{$gene_locations->{$chr}};
}
print "found $genes genes\n";
my $snp_by_lg = {};

foreach my $chromosome (keys %{$plink_hash->{"chromosomes"}}) {
	my $chr_name = "Chr" . sprintf("%02d", $chromosome);
	my $chr_snps = $plink_hash->{"chromosomes"}->{$chromosome};
	my $in_gene = 0;
	my $gff_genes = $gene_locations->{$chr_name};
	my $curr_gene = shift @$gff_genes;
	foreach my $snp (@$chr_snps) {
		# if we have no gene to look at, end.
		if ($curr_gene == undef) {
			last;
		}
		# if the snp we're looking at is past the end of the gene we were looking at, advance the gene
		while ($snp > $curr_gene->{"end"}) {
			$in_gene = 0;
			$curr_gene = shift @$gff_genes;
		}

		# okay. Now for the current gene possibility:
		my $gff_start = $curr_gene->{"start"};
		my $gff_end = $curr_gene->{"end"};
		my $gff_name = $curr_gene->{"name"};
		if (!(exists $snp_by_lg->{$chr_name}->{$gff_name})) {
			my @snp_gene = ();
			$snp_by_lg->{$chr_name}->{$gff_name} = \@snp_gene;
		}
		# if the snp is between the gene's start and end, it's in the gene.
		if (($snp >= $gff_start) && ($snp <= $gff_end)) {
			push @{$snp_by_lg->{$chr_name}->{$gff_name}}, $snp;

			# the mapped snp with name $chromosome_$snp will belong to gene $gff_name
			my $snp_name = $chromosome."_".$snp;
			$plink_hash->{"snps"}->{$snp_name}->{"gene"} = $gff_name;
		} else {
			print $chromosome."_$snp is not in $gff_name ($gff_start - $gff_end)\n";
		}
	}
}


# write out as nexus:
print "writing output to $outfile\n";

my $nexushash = {};

# set up taxa block:
$nexushash->{"taxa"} = $plink_hash->{"names"};

# set up characters block:
$nexushash->{"characters"} = {};
foreach my $indiv_id (@{$nexushash->{"taxa"}}) {
	$nexushash->{"characters"}->{$indiv_id} = $plink_hash->{"individuals"}->{$indiv_id}->{"genotype"};
}
$nexushash->{"charlabels"} = $plink_hash->{"snp_names"};

# set up sets block:
my $curr_gene = ${$plink_hash->{"snp_names"}}[0];
my $curr_gene_start = 0;
my @charsets = ();
$nexushash->{"sets"}->{"charset"} = \@charsets;
for (my $i = 0; $i < @{$plink_hash->{"snp_names"}}; $i++) {
	my $snp_name = ${$plink_hash->{"snp_names"}}[$i];
	my $snp = $plink_hash->{"snps"}->{$snp_name};
	if ($snp->{"gene"} ne $curr_gene) {
# 		print "$curr_gene goes from $curr_gene_start to ". ($i-1) . "\n";
		my $charset = {};
		$charset->{"name"} = $curr_gene;
		$charset->{"start"} = $curr_gene_start + 1;
		$charset->{"end"} = $i;
		push @charsets, $charset;
		$curr_gene = $snp->{"gene"};
		$curr_gene_start = $i;
	}
}
my $nexusstring = "#NEXUS\n\n";
$nexusstring .= write_nexus_taxa_block($nexushash);
$nexusstring .= write_nexus_character_block($nexushash);
$nexusstring .= write_nexus_sets_block($nexushash);

open OUT_FH, ">", $outfile;
print OUT_FH $nexusstring;
close OUT_FH;

__END__

=head1 NAME

plink_to_nexus

=head1 SYNOPSIS

plink_to_nexus [-map mapfile -ped pedfile] [-input inputname] [-output outputname]


=head1 OPTIONS
    -input:         filename of ped/map file (if both share a name w/o the file extension)
    -ped:           filename of ped file (must specify -map as well)
    -map:           filename of map file (must specify -ped as well)
	-outputfile:    name of output file (will have extension .nex)

=head1 DESCRIPTION

Takes a pair of plink-formatted .map/.ped files and converts them to a nexus file.

=cut

