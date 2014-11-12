#!/usr/bin/env perl

use strict;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/..";
use Subfunctions qw(parse_phylip pad_seq_ends debug set_debug consensus_str);
use Nexus qw(write_nexus_character_block write_nexus_trees_block write_nexus_taxa_block);

my $help = 0;
my $outfile = "";
my $inputname = "";

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

GetOptions ('input=s' => \$inputname,
			'output=s' => \$outfile,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help){
    pod2usage(-verbose => 1);
}

if ($outfile eq "") {
	$outfile = "$inputname.nex";
}

if ($outfile !~ /\.nex$/) {
	$outfile = "$outfile.nex";
}

my $raxml_data = {};
$raxml_data->{"name"} = "RAxML_%.$inputname";
$raxml_data->{"input"} = $inputname;

if (!(-x "RAxML_info.$inputname")) {
	print "running RAxML for the input file...\n";
	if ($inputname =~ /\.fa.*$/) {
		# it's a fasta file, convert to phylip...
	}
}

print "found the RAxML run\n";

open FH, "<", "RAxML_info.$inputname";
$raxml_data->{"params"} = "";
# $raxml_data->{"input"} = "";
my $line = readline FH;
while ($line) {
	if ($line =~ /^\s+$/) {
		# go on
	} elsif ($line =~ /RAxML was called as follows/) {
		while ($line !~ /-s\s(.+?)\s/) {
			$line = readline FH;
		}
		$line =~ /-s\s(.+?)\s/;
		$raxml_data->{"input"} = $1;
		print "$raxml_data->{input} is the input file\n";
	} elsif ($raxml_data->{"input"} eq "") {
		$raxml_data->{"params"} .= $line;
	}
	$line = readline FH;
}
close FH;

##### write in the input file
my ($taxa, $taxon_names) = parse_phylip($raxml_data->{"input"});
$raxml_data->{"taxa"} = $taxon_names;
$raxml_data->{"characters"} = $taxa;


##### write in the best tree
open FH, "<", "RAxML_bipartitions.$inputname";
$raxml_data->{"besttree"} = "";
foreach my $line (<FH>) {
	$raxml_data->{"besttree"} .= $line;
}
close FH;


##### write out the NEXUS file
open OUT_FH, ">", $outfile;
print OUT_FH "#NEXUS\n\n";

# print OUT_FH  write_nexus_character_block ($raxml_data->{"characters"}, $raxml_data->{"taxa"});


my @lines = split(/\n/,$raxml_data->{"params"});
foreach my $line (@lines) {
	print OUT_FH "[ ".$line." ]\n";
}

print OUT_FH write_nexus_trees_block ($raxml_data->{"besttree"}, $raxml_data->{"taxa"});

close OUT_FH;

__END__

=head1 NAME

lookup

=head1 SYNOPSIS

lookup.pl -lookup lookupfile -subject subjectfile -match matchcol

=head1 OPTIONS

  -lookup:          list of items to find.
  -subject:         tab-delimited file to find items in.
  -match:           column number of the subject file that is to be matched (1-indexed)

=head1 DESCRIPTION

Return only the requested items from a subject tab-delimited table file.

=cut

