use strict;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
require "subfuncs.pl";

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my $runline = "running " . basename($0) . " " . join (" ", @ARGV) . "\n";

my @inputfiles = ();
my $outname = "melded";
my $outformat = "nex";
my $help = 0;

GetOptions ('files|input=s{2,}' => \@inputfiles,
            'outputfile=s' => \$outname,
            'format=s' => \$outformat,
            'help' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

print $runline;

my ($res1, $res2) = meld_sequence_files(\@inputfiles);
my %mastertaxa = %{$res1};
my %regiontable = %{$res2};

my $ntax = keys %mastertaxa;
my $nchar = delete $mastertaxa{"length"};

$regiontable{ "regions" } = "";
$regiontable{ "exclusion-sets" } = "";

open (fileOUT, ">$outname"."_regions.tab") or die "couldn't make $outname"."_regions.tab\n";
print fileOUT "regions\t$regiontable{'regions'}\n";
print fileOUT "exclusion_sets\t$regiontable{'exclusion-sets'}\n";
delete $regiontable{"regions"};
delete $regiontable{"exclusion-sets"};
foreach my $key (keys %regiontable) {
	print fileOUT "$key\t$regiontable{$key}\n";
}
close fileOUT;

#print melded matrix to output file
if ($outformat =~ /fa/) {
	open (fileOUT, ">$outname".".".$outformat) or die "couldn't make $outname"."."."$outformat\n";
	truncate fileOUT, 0;
	foreach my $key ( keys %mastertaxa ) {
		my $value = $mastertaxa{$key};
		print fileOUT ">$key\n$value\n";
	}
	close fileOUT;
} else {
	open (fileOUT, ">$outname".".nex") or die "couldn't make $outname".".nex\n";
	truncate fileOUT, 0;
	print fileOUT "#NEXUS\nBegin DATA;\n";
	print fileOUT "Dimensions ntax=$ntax nchar=$nchar;\n";
	print fileOUT "Format datatype=NUCLEOTIDE gap=-;\nMatrix\n";
	foreach my $key ( keys %mastertaxa ) {
		my $value = $mastertaxa{$key};
		print fileOUT "$key\t$value\n";
	}
	print fileOUT ";\nEnd;\n";
	close fileOUT;
}

print "Melded matrix is in $outname".".$outformat; summary of regions is in $outname"."_regions.tab\n";



__END__

=head1 NAME

meld

=head1 SYNOPSIS

meld -input [at least two files to meld] [-output outfile] [-format nex|fasta]

=head1 OPTIONS

    -files|input:   list of files to meld
	-output:        optional: name of output file, default is "melded"
	-format:		optional: format of output, default is nexus.

=head1 DESCRIPTION

Melds aligned sequence files of either fasta or nexus format.
Sequences with the same names will be concatenated; sequences missing
from any particular alignment will be concatenated as missing data.
A tab-delimited table of the aligned regions and the sequences present
in each will be output to outfile_regions.tab.

The program defaults to output into melded.nex and melded_regions.tab.

=cut


