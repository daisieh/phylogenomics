#!/usr/bin/perl
use strict;
sub uninterleave;
my @inputfiles;
push @inputfiles, @ARGV;

my %matrices = ();
my %mastertaxa = ();
my %regiontable = ();
$regiontable{ "regions" } = "";
$regiontable{ "exclusion-sets" } = "";
my $currlength = 0;
my $outfolder = "./";
my @matrixnames = ();

#foreach my $inputfile (@inputfiles) {
for (my $i=0; $i< scalar(@inputfiles); $i++) {
	my $inputfile = @inputfiles[$i];
	if ($inputfile eq "-output") {
		$outfolder = @inputfiles[$i+1];
		$i++;
		next;
	}
	push @matrixnames, $inputfile;
	$matrices{ $inputfile } = uninterleave $inputfile;
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
foreach my $key ( keys %matrices ) {
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

open (fileOUT, ">$outfolder". "regions.tab") or die "couldn't make $outfolder/regions.tab\n";
truncate fileOUT, 0;
print fileOUT "regions\t$regiontable{'regions'}\n";
print fileOUT "exclusion_sets\t$regiontable{'exclusion-sets'}\n";
delete $regiontable{"regions"};
delete $regiontable{"exclusion-sets"};
foreach my $key (keys %regiontable) {
	print fileOUT "$key\t$regiontable{$key}\n";
}
close fileOUT;

#print melded matrix to melded.nex
open (fileOUT, ">$outfolder/melded.nex") or die "couldn't make melded.nex\n";
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

print "Melded matrix is in melded.nex; summary of regions is in regions.tab\n";

sub uninterleave {
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

