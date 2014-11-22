#!/usr/bin/env perl

use strict;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use File::Temp qw(tempfile);
use File::Basename qw(fileparse);
use FindBin;
use lib "$FindBin::Bin/../lib";
use Nexus qw(write_nexus_character_block write_nexus_trees_block write_nexus_taxa_block);
use Subfunctions qw(write_phylip parse_phylip parse_fasta pad_seq_ends debug set_debug consensus_str);

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
my $running = 0;
($inputname, my $inputpath, undef) = fileparse ($inputname);
my $raxmlinfofile = File::Spec->catpath( undef, $inputpath, "RAxML_info.$inputname" );
print "looking for existing RAxML run $raxmlinfofile...\n";
if (!(-e $raxmlinfofile)) {
	print "running RAxML for the input file $inputname...\n";
	$running = 1;
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
	$inputname = fileparse ($raxml_input);
# 	raxmlHPC-PTHREADS -fa -s all_cps.phy -x 141105 -# 100 -m GTRGAMMA -n 141105 -T 16 -p 141105
	my $cmd = "raxmlHPC-PTHREADS -fa -s $raxml_input -x 141105 -# 100 -m GTRGAMMA -n $inputname -T 16 -p 141105";
	print $cmd . "\n";
	system ($cmd);
	$raxmlinfofile = File::Spec->catpath( undef, $inputpath, "RAxML_info.$inputname" );
}


if (!(-s $raxmlinfofile)) {
	print "Couldn't find the RAxML run $raxmlinfofile\n";
	exit;
}
print "found the RAxML run $raxmlinfofile\n";

open FH, "<", $raxmlinfofile;
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
	if (-e $raxml_data->{"input"}) {
		print "loading input file $raxml_data->{input}\n";
		my ($taxa, $taxon_names) = parse_phylip($raxml_data->{"input"});
		$raxml_data->{"taxa"} = $taxon_names;
		$raxml_data->{"characters"} = $taxa;
	} else {
		print "Input file $raxml_data->{input} could not be found.\n";
		exit;
	}
}


##### trees:
$raxml_data->{"trees"} = ();

##### write in the best tree
print "loading in best tree from RAxML_bipartitions.$inputname\n";
open FH, "<", "RAxML_bipartitions.$inputname";
my $tree = {};
$tree->{"besttree"} = "";
foreach my $line (<FH>) {
	$tree->{"besttree"} .= $line;
}
push @{$raxml_data->{"trees"}}, $tree;
close FH;

#### write in the bootstrap trees
print "loading in bootstrap trees from RAxML_bootstrap.$inputname\n";
open FH, "<", "RAxML_bootstrap.$inputname";
my $i = 1;
foreach my $line (<FH>) {
	$tree = {};
	$tree->{"bootstrap$i"} = $line;
	push @{$raxml_data->{"trees"}}, $tree;
	$i++;
}
close FH;

##### write out the NEXUS file
print "writing out nexus file $outfile\n";
open OUT_FH, ">", $outfile;
print OUT_FH "#NEXUS\n\n";

# print OUT_FH write_nexus_character_block ($raxml_data);
# print OUT_FH write_nexus_taxa_block ($raxml_data);

my @lines = split(/\n/,$raxml_data->{"params"});
foreach my $line (@lines) {
	print OUT_FH "[ ".$line." ]\n";
}
if (exists $raxml_data->{"fulltaxa"}) {
	$raxml_data->{"alttaxa"} = ();
	push @{$raxml_data->{"alttaxa"}}, $raxml_data->{"taxa"};
	$raxml_data->{"taxa"} = delete $raxml_data->{"fulltaxa"};
}

print OUT_FH write_nexus_trees_block ($raxml_data);

close OUT_FH;

if ($running) {
	system ("rm *.$inputname");
}

__END__

=head1 NAME

lookup

=head1 SYNOPSIS

raxml.pl -input raxml_input -out nexus_output

=head1 OPTIONS

  -input:          list of items to find.
  -output:         tab-delimited file to find items in.

=head1 DESCRIPTION

Return only the requested items from a subject tab-delimited table file.

=cut

