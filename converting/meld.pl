#!/usr/bin/perl
use strict;
sub parse_nexus;
sub parse_fasta;
my @inputfiles;
push @inputfiles, @ARGV;

my %matrices = ();
my %mastertaxa = ();
my %regiontable = ();
$regiontable{ "regions" } = "";
$regiontable{ "exclusion-sets" } = "";
my $currlength = 0;
my $outname = "melded";
my $outformat = "nexus";
my @matrixnames = ();

if ((@ARGV == 0) || (@ARGV[1] =~ /-h/) || (@ARGV[1] =~ /-u/)) {
	print "meld.pl [-output outfile] [-format nex|fasta] inputfile1 inputfile2...\n";
	print qq{
Melds aligned sequence files of either fasta or nexus format.
Sequences with the same names will be concatenated; sequences missing
from any particular alignment will be concatenated as missing data.
A tab-delimited table of the aligned regions and the sequences present
in each will be output to outfile_regions.tab.

The program defaults to output into melded.nex and melded_regions.tab.
};
	exit;
}


for (my $i=0; $i< scalar(@inputfiles); $i++) {
	my $inputfile = @inputfiles[$i];
	if ($inputfile eq "-output") {
		$outname = @inputfiles[$i+1];
		$i++;
		next;
	} elsif ($inputfile eq "-format") {
		$outformat = @inputfiles[$i+1];
		$i++;
		next;
	}

	push @matrixnames, $inputfile;
	if ($inputfile =~ /\.nex/) {
		$matrices{ $inputfile } = parse_nexus $inputfile;
	} elsif ($inputfile =~/\.fa/) {
		$matrices{ $inputfile } = parse_fasta $inputfile;
	} else {
		print "Couldn't parse $inputfile: not nexus or fasta format\n";
	}
	while ( my ($key, $value) = each(%{$matrices{ $inputfile }}) ) {
        $mastertaxa{ $key } = "";
    }
}

#for each matrix, make a mastertaxa matrix that has missing data for the taxa with no entries in this matrix.
#also, add another column to the regiontable hash
# foreach my $key ( keys %matrices ) {
for (my $i=0; $i<@matrixnames; $i++) {
	my $key = @matrixnames[$i];
	my $ref = $matrices{$key};
	print "adding $key...\n";
	$regiontable{"regions"} = $regiontable{"regions"} . "$key\t";
	my %expandedmatrix = %mastertaxa;
	foreach my $k (keys %{$ref}) {
		#add entries from this matrix into expandedmatrix
		$expandedmatrix{ $k } = $ref->{ $k };
		$regiontable{$k} = $regiontable{$k} . "x\t";
	}
	my $total = $expandedmatrix{'length'};
	my $starts_at = $currlength + 1;
	$currlength = $currlength + $total;
	$regiontable{"exclusion-sets"} = $regiontable{"exclusion-sets"} . "$starts_at" . "-" . "$currlength\t";
	my $replacement = "-" x $total;
	foreach my $k (keys %expandedmatrix) {
		my $l = length($expandedmatrix{$k});
		if ($l == 0) {
			#if the entry in expandedmatrix is empty, fill it with blanks
			$expandedmatrix{ $k } = "$replacement";
			$regiontable{$k} = $regiontable{$k} . "\t";
		}
		$l = length($expandedmatrix{$k});
	}
	$matrices{$key} = \%expandedmatrix;
}

#now, for each matrix, concatenate them to the corresponding entry in mastertaxa
my $l = 0;
# foreach my $key ( keys %matrices ) {
for (my $i=0; $i<@matrixnames; $i++) {
	my $key = @matrixnames[$i];
	my $ref = $matrices{$key};
	$l = $l + $ref->{"length"};
	foreach my $k (keys %{$ref}) {
		$mastertaxa{$k} = $mastertaxa{$k} . $ref->{$k};
	}
}

delete $mastertaxa{"length"};
delete $regiontable{"length"};
my $ntax = keys %mastertaxa;
my $nchar = $l;

open (fileOUT, ">$outname"."_regions.tab") or die "couldn't make $outname"."_regions.tab\n";
truncate fileOUT, 0;
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

sub parse_fasta {
	my $inputfile = shift(@_);

	my %taxa = ();
	open (fileIN, "$inputfile") or die "no file named $inputfile";

	my $input = readline fileIN;
	my $length = 0;
	my $taxonlabel = "";
	my $sequence = "";
	while ($input ne "") {
		if ($input =~ /^>(.+)\s*$/) {
			$taxonlabel = $1;
			if ($length > 0) {
				# we are at the next taxon; push the last one onto the taxon array.
				$taxa {"length"} = $length;
				$length = 0;
			}
		} else {
			$input =~ /^\s*(.+)\s*$/;
			$taxa {$taxonlabel} .= $1;
			$length += length($1);
		}
		$input = readline fileIN;
	}

	close (fileIN);
	return \%taxa;
}

sub parse_nexus {
	my $inputfile = shift(@_);
	open (fileIN, "$inputfile") or die "no file named $inputfile";
	my @inputs = <fileIN>;

	my $input = "";
	if ($inputs[1] eq "") {
		$input = $inputs[0];
		$input =~ s/\r/\n/gs;
	} else {
		foreach my $line (@inputs) {
			$input .= "$line";
		}
	}
	close (fileIN);

	#remove comment blocks
	$input =~ s/\[.*?\]//sg;

	$input =~ /Format(.*?)\;/ig;
	my $format = "$1";
	$format =~ /gap=(.)/;
	my $gapchar = $1;

	#parse nexus block
	$input =~ /Matrix(.*?)\;/isg;
	my $matrix = "$1";

	$matrix =~ /\s*?(\S+?)\s+/s;
	my $firsttaxon = "$1";
	my %taxa = ();

	my @sections = split /$firsttaxon/, $matrix;
	my @taxonlabels;
	foreach my $section (@sections) {
		$section =~ s/\s+$//s;
		if ($section eq "") { next; }
		$section = "$firsttaxon$section";
		$section =~ s/\s+$/\n/s;
		$section =~ s/\t//sg;
		my $numtaxa = ($section =~ s/\n/\n/sg);

		@taxonlabels = split /\n/, $section;

		foreach my $taxonlabel (@taxonlabels) {
			$taxonlabel =~ s/\s+(.+?)$//;
			my $taxondata = $1;
			$taxondata =~ s/$gapchar/-/g;
			$taxa{ $taxonlabel } = $taxa{ $taxonlabel } . $taxondata;
		}
	}

	my $length = length $taxa{ $taxonlabels[0] };
	$taxa{ "length" } = $length;

	return \%taxa;
}

