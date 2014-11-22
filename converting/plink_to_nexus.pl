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
use Nexus qw(write_nexus_character_block write_nexus_trees_block write_nexus_taxa_block);
use Plink qw(parse_plink);

my $help = 0;
my $outfile = "";
my $inputmap = "";
my $inputped = "";
my $inputname = "";

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

GetOptions ('map=s' => \$inputmap,
			'ped=s' => \$inputped,
			'input=s' => \$inputname,
			'output=s' => \$outfile,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help){
    pod2usage(-verbose => 1);
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

unless (-e $inputped) {
	pod2usage(-msg => "File $inputped does not exist.", -exitval => 2);
}

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

my $nexusstring = "#NEXUS\n\n";
$nexusstring .= write_nexus_taxa_block($nexushash);
$nexusstring .= write_nexus_character_block($nexushash);

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

