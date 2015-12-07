# written for the Graham lab on June 6, 2013

# if executed in the directory with a bunch of fasta files that have a
# list of genes per species in the format ">rpl32_Temp_temp", it will
# make a set of gene.fasta files in a subdirectory with each species's gene in it.

use strict;

my $species_file = shift;
my $result_dir = shift;

unless ($species_file) {
	my $filename = "tempfile";
	system ("ls > ".$filename);
	$species_file = $filename;
}

unless ($result_dir) {
	$result_dir = "genes";
}
unless (-d $result_dir) {
	system ("mkdir $result_dir");
}

open FH, "<", $species_file;
my @list = <FH>;
my %species_hash = ();
my @file_list;
foreach my $sp (@list) {
	chomp $sp;
	if (-e $sp) {
		if ($sp =~ /(.*?)\.fa.*/) {
			push @file_list, "$sp";
		}
	} elsif (-e "$sp.fasta") {
		push @file_list, "$sp.fasta";
	}
}
close FH;

foreach my $file (@file_list) {
	open FH, "<", "$file";
	my $seqname = readline FH;
	my $seq = readline FH;

	$seqname =~ />(.*)_(.*_.*)$/;
	my $gene_name = $1;
	my $species_name = $2;

	while ($seq) {
		$seqname =~ />(.*)_$species_name$/;
		$gene_name = $1;
		$species_hash{$gene_name}{$species_name} = $seq;
		$seqname = readline FH;
		$seq = readline FH;
	}
}
close FH;

my ($gene, $genehash);
my $species_list;
while (($gene, $genehash) = each %species_hash) {
	my ($sp_name, $seq);
	while (($sp_name, $seq) = each %{$genehash}) {
 		open FH, ">>", "$result_dir/$gene.fasta";
 		print FH ">$sp_name\n$seq\n";
 		close FH;
	}
}
