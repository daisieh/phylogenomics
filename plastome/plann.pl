#!/usr/bin/env perl

use Getopt::Long;
use Pod::Usage;
use File::Temp qw (tempfile tempdir);
use FindBin;
use lib "$FindBin::Bin/../lib";
use Blast qw (parse_xml);
use Genbank qw (parse_genbank parse_regionfile write_features_as_table parse_gene_array_to_features set_sequence get_sequence);
use Subfunctions qw (parse_fasta blast_to_genbank merge_to_featuretable);
use Data::Dumper;

my $help = 0;
my $outfile = "";
my $gbfile = "";
my $fastafile = "";
my $orgname = "";
my $samplename = "";

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

GetOptions ('reference=s' => \$gbfile,
			'fastafile=s' => \$fastafile,
			'outfile=s' => \$outfile,
			'organism=s' => \$orgname,
			'sample=s' => \$samplename,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

if ($gbfile !~ /\.gb$/) {
	print "reference file needs to be a fully annotated Genbank file.\n";
	exit;
}

if ($outfile eq "") {
	$outfile = "output";
}


my $gene_array = parse_genbank($gbfile);
open FEATURES_FH, ">", "$gbfile.features";
print FEATURES_FH write_features_as_table ($gene_array);
close FEATURES_FH;

# writes out $outfile.regions
my ($ref_hash, $ref_array) = blast_to_genbank ($gene_array, $fastafile);
open my $outfh, ">", "$outfile.regions" or die "couldn't create $outfile";
foreach my $subj (@$ref_array) {
	print $outfh "$subj($ref_hash->{$subj}->{'strand'})\t$ref_hash->{$subj}->{'start'}\t$ref_hash->{$subj}->{'end'}\n";
}

close $outfh;


my ($fastahash, $fastaarray) = parse_fasta($fastafile);

# there should be only one key, so just one name.
if ($samplename eq "") {
	$samplename = @$fastaarray[0];
}

my $queryseq = $fastahash->{@$fastaarray[0]};
my $genbank_header = "$samplename [organism=$orgname][moltype=Genomic DNA][location=chloroplast][topology=Circular][gcode=11]";
open FASTA_FH, ">", "$outfile.fsa";
print FASTA_FH ">$genbank_header\n$queryseq\n";
close FASTA_FH;

merge_to_featuretable ("$outfile.regions", $fastafile, "$gbfile.features", $outfile, $genbank_header);

# a regionfile is the output of parse_blast.pl comparing the fastafile to the reference fasta file from genbank.pl
# my ($gene_array, $gene_index_array) = parse_regionfile("$outfile.regions");
#
# my ($destination_gene_array, $destination_gene_index_array) = Genbank::parse_featurefile("$gbfile.features");
#
#
# my $gene_hash = {};
# foreach my $id (@$gene_index_array) {
# 	my $gene = shift $gene_array;
# 	$gene_hash->{$id} = $gene;
# }
#
# my $dest_gene_hash = {};
# foreach my $id (@$destination_gene_index_array) {
# 	my $gene = shift $destination_gene_array;
# 	$dest_gene_hash->{$id} = $gene;
# }
#
# # fill in the genes from the regionfile with the info from the destination gene array
# my @final_gene_array = ();
# foreach my $id (@$gene_index_array) {
# 	my $gene = $gene_hash->{$id};
# 	my $dest_gene = $dest_gene_hash->{$id};
#
# 	foreach my $q (keys %{$dest_gene->{"qualifiers"}}) {
# 		if (!(exists $gene->{"qualifiers"}->{$q})) {
# 			$gene->{"qualifiers"}->{$q} = $dest_gene->{"qualifiers"}->{$q}
# 		}
# 	}
# 	my @new_contains = ();
# 	$gene->{"id"} = $id;
# 	foreach my $destcontains (@{$dest_gene->{"contains"}}) {
# 		my $genecontains = shift $gene->{"contains"};
# 		$destcontains->{"region"} = $genecontains->{"region"};
# 		push @new_contains, $destcontains;
# 	}
# 	$gene->{"contains"} = \@new_contains;
# 	push @final_gene_array, $gene;
# }
#
# # print header
# open FH, ">", "$outfile.tbl";
# print FH ">Features\t$genbank_header\n";
#
# # start printing genes
# foreach my $gene (@final_gene_array) {
# 	# first, print overall gene information
# 	my $genename = $gene->{"qualifiers"}->{"gene"};
# 	my $gene_id = $gene->{"id"};
# 	my $feat_id = 0;
# 	foreach my $r (@{$gene->{'region'}}) {
# 		$r =~ /(\d+)\.\.(\d+)/;
# 		print FH "$1\t$2\tgene\n";
# 	}
# 	foreach my $q (keys %{$gene->{'qualifiers'}}) {
# 		print FH "\t\t\t$q\t$gene->{qualifiers}->{$q}\n";
# 	}
#
# 	# then, print each feature contained.
# 	foreach my $feat (@{$gene->{'contains'}}) {
# 		foreach my $reg (@{$feat->{"region"}}) {
# 			my $strand = "+";
# 			my ($start, $end) = split (/\.\./, $reg);
# 			if ($end < $start) {
# 				$strand = "-";
# 				my $oldend = $start;
# 				$start = $end;
# 				$end = $oldend;
# 			}
# 			my $regseq = Genbank::sequence_for_interval ($reg);
# 			my $featname = $feat->{"type"};
# 			$feat_id++;
# 		}
# 		print FH Genbank::sequin_feature ($feat->{'region'}, $feat);
# 	}
# }
# close FH;

# need to annotate inverted repeats

# my $self_array = parse_xml ("$refblast.xml");
# my @irs = ();
# foreach my $hit (@$self_array) {
# 	my @hsps = sort order_by_query_start @{$hit->{"hsps"}};
# 	foreach my $hsp (@hsps) {
# 		# only look at identical pieces that are smaller than the entire reference
# 		if ((($hsp->{"query-to"} - $hsp->{"query-from"}) < ($reflen - 1)) && (($hsp->{"query-to"} - $hsp->{"query-from"}) > 10000)) {
# 			push @irs, $hsp;
# 		}
# 	}
# }
#
# if (@irs > 2) {
# 	die "Error! There seem to be more than two inverted repeats. Are you sure this is a plastome sequence?";
# }
# sub order_by_query_start {
# 	my $bstart = $b->{"query-from"};
# 	my $astart = $a->{"query-from"};
#
# 	if ($astart < $bstart) { return -1; }
# 	if ($astart > $bstart) { return 1; }
# 	return 0;
# }


__END__

=head1 NAME

plann.pl

=head1 SYNOPSIS

plann.pl -reference gbfile.gb -fasta denovoplastome.fasta -out outfile [-organism "Genus species"] [-sample samplename]

=head1 OPTIONS

  -reference:       a well-annotated plastome reference sequence, in genbank file format
  -fastafile:       the plastome sequence to be annotated, in fasta format
  -outfile:         the output name (default is "output")
  -organism:        [optional: scientific name for Genbank annotation]
  -sample:          the name of the plastome sample (default is the name in the fasta file)

=head1 DESCRIPTION


=cut
