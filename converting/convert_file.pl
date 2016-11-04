#!/usr/bin/env perl

use Bio::SeqIO;
use Bio::Seq::RichSeq;
use Bio::Align::Utilities qw(cat);
use Bio::AlignIO;
use FindBin;
use lib "$FindBin::Bin/../lib";
use Bioperl qw(convert_aln_to_nexus);

my $usage  = "convert_file.pl in_file out_file\n";
my $in_file   = shift or die $usage;
my $out_file = shift or die $usage;

$in_file =~ /.*\.(.+)/;
my $in_format = $1;

$out_file =~ /.*\.(.+)/;
my $out_format = $1;

if (($out_format eq "fa") || ($out_format eq "fas") || ($out_format eq "fasta")) {
	$out_format = "fasta";
} elsif (($out_format eq "nex") || ($out_format eq "nexus")) {
	$out_format = "nexus";
} elsif (($out_format eq "gb") || ($out_format eq "genbank")) {
	$out_format = "genbank";
} elsif (($out_format eq "phy") || ($out_format eq "phylip")) {
	$out_format = "phylip";
}

my $aln;

if (($in_format eq "gb") || ($in_format eq "genbank")) {
	my $seqio_object = Bio::SeqIO->new(-file => $in_file);
	my $loc_aln = Bio::SimpleAlign->new();
	# write each entry in the input file to the output file
	while (my $inseq = $seqio_object->next_seq) {
		my $loc_seq = Bio::LocatableSeq->new(-seq => $inseq->seq(), -id => $inseq->id, -start => 1, -end => $inseq->length());
		$loc_aln->add_seq($loc_seq);
	}
	$aln = $loc_aln;
} else {
	my $guesser = Bio::Tools::GuessSeqFormat->new( -file => $in_file );
	$in_format  = $guesser->guess;
	my $loc_aln = Bio::AlignIO->new(-file => $in_file, -format => $in_format);
	$aln = $loc_aln->next_aln();
}

if ($out_format eq "nexus") {
	my $result = convert_aln_to_nexus ($aln);
	open my $gene_file, ">$out_file";
	truncate $gene_file, 0;
	print $gene_file $result;
	close $gene_file;
} elsif ($out_format eq "genbank") {
	my $name = $in_file;
	$name =~ s/\..+//;
	$name =~ s/.*\///;
	my $seq_out = Bio::SeqIO->new(-file => ">$out_file", -format => $out_format);
	foreach my $seq ($aln->each_seq()) {
		my $name = $seq->id() . "_$name";
		my $tempseq = Bio::Seq::RichSeq->new(-seq => $seq->seq(), -id => $name, -accession_number => "");
		$seq_out->write_seq($tempseq);
	}
} else {
	open my $fh, ">$out_file";
	my $aln_out = Bio::AlignIO->new(-fh => $fh, -format => "phylip");
	$aln_out->write_aln($aln);
}
