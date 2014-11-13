#!/usr/bin/env perl

use strict;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use File::Temp qw(tempfile);
use File::Basename qw(fileparse);
use FindBin;
use lib "$FindBin::Bin/..";
use Subfunctions qw(write_phylip parse_phylip parse_fasta pad_seq_ends debug set_debug consensus_str);
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
# $raxml_data->{"name"} = "RAxML_%.$inputname";
# $raxml_data->{"input"} = $inputname;

if (!(-s "RAxML_info.$inputname")) {
	print "running RAxML for the input file $inputname...\n";
	if ($inputname =~ /\.fa.*$/) {
		# it's a fasta file, convert to phylip...
		my ($taxa, $taxanames) = parse_fasta($inputname);
		$raxml_data->{"fulltaxa"} = $taxanames;
		$raxml_data->{"characters"} = $taxa;
	} elsif ($inputname =~ /\.phy.*$/) {
		# it's a phylip file...
		my ($taxa, $taxanames) = parse_phylip($inputname);
		$raxml_data->{"fulltaxa"} = $taxanames;
		$raxml_data->{"characters"} = $taxa;
	}

	# write out a temporary phylip file for input:
	my ($fh, $raxml_input) = tempfile();
	my $phylip_str = write_phylip ($raxml_data->{"characters"}, $raxml_data->{"fulltaxa"});
	print $fh $phylip_str;
	close $fh;
	# make sure that the phylip taxa names are accounted for:
	my @phynames = split (/\n/, $phylip_str, @{$raxml_data->{"fulltaxa"}} + 2);
	shift @phynames;
	pop @phynames;
	$raxml_data->{"taxa"} = ();
	foreach my $name (@phynames) {
		$name =~ /^(.+?)\s/;
		push @{$raxml_data->{"taxa"}}, $1;
	}
	print Dumper($raxml_data->{"taxa"});
	$inputname = fileparse ($raxml_input);
# 	raxmlHPC-PTHREADS -fa -s all_cps.phy -x 141105 -# 100 -m GTRGAMMA -n 141105 -T 16 -p 141105
	my $cmd = "raxmlHPC-PTHREADS -fa -s $raxml_input -x 141105 -# 100 -m GTRGAMMA -n $inputname -T 16 -p 141105";
	print $cmd . "\n";
	system ($cmd);
}

if (!(-s "RAxML_info.$inputname")) {
	print "Couldn't find the RAxML run RAxML_info.$inputname\n";
	exit;
}
print "found the RAxML run RAxML_info.$inputname\n";

open FH, "<", "RAxML_info.$inputname";
$raxml_data->{"params"} = "";
$raxml_data->{"input"} = "";
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
		print "$raxml_data->{input} is the input file\n$line";
	} elsif ($raxml_data->{"input"} eq "") {
		$raxml_data->{"params"} .= $line;
	}
	$line = readline FH;
}
close FH;

##### write in the input file if it hasn't already been read:
if (!(exists $raxml_data->{"characters"})) {
	my ($taxa, $taxon_names) = parse_phylip($raxml_data->{"input"});
	$raxml_data->{"taxa"} = $taxon_names;
	$raxml_data->{"characters"} = $taxa;
}

##### write in the best tree
open FH, "<", "RAxML_bipartitions.$inputname";
$raxml_data->{"trees"}->{"besttree"} = "";
foreach my $line (<FH>) {
	$raxml_data->{"trees"}->{"besttree"} .= $line;
}
close FH;


##### write out the NEXUS file
open OUT_FH, ">", $outfile;
print OUT_FH "#NEXUS\n\n";

# print OUT_FH  write_nexus_character_block ($raxml_data->{"characters"}, $raxml_data->{"taxa"});
# print OUT_FH write_nexus_taxa_block ($raxml_data->{"taxa"});

my @lines = split(/\n/,$raxml_data->{"params"});
foreach my $line (@lines) {
	print OUT_FH "[ ".$line." ]\n";
}

print OUT_FH write_nexus_trees_block ($raxml_data->{"trees"}, $raxml_data->{"fulltaxa"}, $raxml_data->{"taxa"});

close OUT_FH;

__END__

=head1 NAME

lookup

=head1 SYNOPSIS

raxml_condense.pl -input raxml_input -out nexus_output

=head1 OPTIONS

  -input:          list of items to find.
  -output:         tab-delimited file to find items in.

=head1 DESCRIPTION

Return only the requested items from a subject tab-delimited table file.

=cut

