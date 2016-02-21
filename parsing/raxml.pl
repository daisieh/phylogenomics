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
my $inputfile = "";
my $raxml = "raxmlHPC-PTHREADS";
my $x = int(rand(1e8));
my $p = int(rand(1e8));
my $threads = 16;
my $arguments = "-# 100 -m GTRGAMMA";

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

GetOptions ('input=s' => \$inputfile,
			'output=s' => \$outfile,
			'raxml=s' => \$raxml,
			'x=i' => \$x,
			'p=i' => \$p,
			'threads=i' => \$threads,
			'arguments=s' => \$arguments,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help){
    pod2usage(-verbose => 2);
}

# Set defaults
if ($outfile eq "") {
	$outfile = "$inputfile.nex";
}

if ($outfile !~ /\.nex$/) {
	$outfile = "$outfile.nex";
}

if ($arguments =~ /^(.*)\s*-T (\d+)\s*(.*)$/) {
	$threads = $2;
	$arguments = $1.$3;
}

if ($arguments =~ /^(.*)\s*-x (\d+)\s*(.*)$/) {
	$x = $2;
	$arguments = $1.$3;
}

if ($arguments =~ /^(.*)\s*-p (\d+)\s*(.*)$/) {
	$p = $2;
	$arguments = $1.$3;
}

my $raxml_data = {};
my $running = 0;
my ($inputname, $inputpath, undef) = fileparse ($inputfile);
my $raxmlinfofile = File::Spec->catpath( undef, $inputpath, "RAxML_info.$inputname" );
print "looking for existing RAxML run $raxmlinfofile...\n";
if (!(-e $raxmlinfofile)) {
	print "running RAxML for the input file $inputname...\n";
	$running = 1;
	if ($inputname =~ /\.fa.*$/) {
		# it's a fasta file, convert to phylip...
		my ($taxa, $taxanames) = parse_fasta($inputfile);
		$raxml_data->{"fulltaxa"} = $taxanames;
		$raxml_data->{"characters"} = $taxa;
	} elsif ($inputname =~ /\.phy.*$/) {
		# it's a phylip file...
		my ($taxa, $taxanames) = parse_phylip($inputfile);
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
	my $cmd = "$raxml -fa -s $raxml_input -n $inputname -x $x -T $threads -p $p $arguments";
	print "RAxML is called as '$cmd'\n";
	system ($cmd);
	$raxmlinfofile = File::Spec->catpath( undef, $inputpath, "RAxML_info.$inputname" );
}


if (!(-s $raxmlinfofile)) {
	print "Couldn't find the RAxML run $raxmlinfofile\n";
	exit;
}
print "found the RAxML run $raxmlinfofile\n";

# metadata that opentree takes:
# <meta datatype="xsd:string" property="ot:branchLengthDescription" xsi:type="nex:LiteralMeta"/>
# <meta datatype="xsd:string" property="ot:branchLengthMode" xsi:type="nex:LiteralMeta">ot:changesCount</meta>
# <meta datatype="xsd:string" property="ot:branchLengthTimeUnit" xsi:type="nex:LiteralMeta">Myr</meta>
# <meta datatype="xsd:string" property="ot:curatedType" xsi:type="nex:LiteralMeta">Maximum likelihood</meta>
# <meta datatype="xsd:string" property="ot:inGroupClade" xsi:type="nex:LiteralMeta">node5</meta>
# <meta datatype="xsd:string" property="ot:nodeLabelMode" xsi:type="nex:LiteralMeta"/>
# <meta datatype="xsd:string" property="ot:nodeLabelTimeUnit" xsi:type="nex:LiteralMeta"/>
# <meta datatype="xsd:string" property="ot:outGroupEdge" xsi:type="nex:LiteralMeta"/>
# <meta datatype="xsd:string" property="ot:specifiedRoot" xsi:type="nex:LiteralMeta">node1</meta>
# <meta datatype="xsd:boolean" property="ot:unrootedTree" xsi:type="nex:LiteralMeta"/>


open FH, "<:crlf", $raxmlinfofile;
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
open FH, "<:crlf", "RAxML_bipartitions.$inputname";
my $tree = {};
$tree->{"besttree"} = "";
foreach my $line (<FH>) {
	$tree->{"besttree"} .= $line;
}
push @{$raxml_data->{"trees"}}, $tree;
close FH;

#### write in the bootstrap trees
print "loading in bootstrap trees from RAxML_bootstrap.$inputname\n";
open FH, "<:crlf", "RAxML_bootstrap.$inputname";
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

raxml.pl

=head1 SYNOPSIS

raxml.pl -input raxml_input -out nexus_output [-raxml raxml_executable] [-x x-seed] [-p p-seed] [-threads num_threads] [-arguments 'raxml-arguments']

=head1 OPTIONS

  -input:          Fasta or Phylip formatted sequence alignment file
  -output:         Name of outputted Nexus file
  -raxml:          optional: the name of the raxml binary (default is 'raxmlHPC-PTHREADS')
  -x:              optional: value of the -x argument (default is a random number)
  -p:              optional: value of the -p argument (default is a random number)
  -threads:        optional: number of threads to be used by raxml (default is 16)
  -arguments:      optional: arguments to be passed to raxml (overrides -x/-p/-T if these are specified in -arguments) (default is "-# 100 -m GTRGAMMA")

=head1 DESCRIPTION

Takes a nucleotide sequence alignment in Fasta or Phylip format, runs RAxML on the sequences, and converts the output into a Nexus-formatted file. Bootstrap trees are saved as "bootstrap1" as well as the best tree as "besttree".

=cut

