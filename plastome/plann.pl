#!/usr/bin/env perl

use Getopt::Long;
use Pod::Usage;
use File::Temp qw (tempfile tempdir);
use FindBin;
use lib "$FindBin::Bin/../lib";
use Blast qw (parse_xml);
use Genbank qw (parse_genbank write_features_as_fasta write_features_as_table parse_regionfile parse_feature_table set_sequence get_sequence sequence_for_interval sequin_feature parse_gene_array_to_features compare_hsps);
use Subfunctions qw (parse_fasta);
use Data::Dumper;

my $help = 0;
my $outfile = "";
my $gbfile = "";
my $fastafile = "";
my $orgname = "";
my $samplename = "";

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

my $gene_array = parse_genbank($gbfile);

open FAS_FH, ">", "$gbfile.fasta";
print FAS_FH write_features_as_fasta ($gene_array);
close FAS_FH;
my ($ref_hash, $ref_array) = parse_fasta ("$gbfile.fasta", 1);
print Dumper ($gene_array);
system("blastn -query $fastafile -subject $gbfile.fasta -outfmt 5 -out $outfile.xml -word_size 10");

print "parsing results to $outfile.regions\n";

my $hit_array = parse_xml ("$outfile.xml");
my $hits = {};
foreach my $hit (@$hit_array) {
	my $subject = $hit->{"subject"}->{"name"};
	my $query = $hit->{"query"}->{"name"};
	my @hsps = sort compare_hsps @{$hit->{"hsps"}};
	my $best_hit = shift @hsps;
	$hits->{$subject}->{"hsp"} = $best_hit;
	$hits->{$subject}->{"slen"} = $hit->{"subject"}->{"length"};
	if ($best_hit->{"hit-from"} < $best_hit->{"hit-to"}) {
		$hits->{$subject}->{"orientation"} = 1;
	} else {
		$hits->{$subject}->{"orientation"} = -1;
	}
}

open REGIONS_FH, ">", "$outfile.regions" or die "couldn't create $outfile.regions";

foreach my $subj (@$ref_array) {
	$subj =~ s/\t.*$//;
	if (!(exists $hits->{$subj}->{"hsp"})) {
		next;
	}
	my $adjusted_hseq = $hits->{$subj}->{"hsp"}->{"hseq"};
	$adjusted_hseq =~ s/-//g;
	my $hlen = length $adjusted_hseq;
	print REGIONS_FH "$subj\t$hits->{$subj}->{hsp}->{'query-from'}\t$hits->{$subj}->{hsp}->{'query-to'}\n";
	print REGIONS_FH $ref_hash->{$subj} . "\n";
	print REGIONS_FH "$hits->{$subj}->{hsp}->{'hseq'}\n$hits->{$subj}->{hsp}->{'qseq'}\n";
}

close REGIONS_FH;

# a regionfile is the output of parse_blast.pl comparing the fastafile to the reference fasta file from genbank.pl
my ($gene_array2, $gene_index_array) = parse_regionfile("$outfile.regions");

my $featuretable = write_features_as_table ($gene_array);
my ($destination_gene_hash_array, $destination_gene_index_array) = parse_gene_array_to_features ($gene_array);

my ($fastahash, undef) = parse_fasta($fastafile);
my $seqlen = 0;
foreach my $k (keys $fastahash) {
	# there should be only one key, so just one name.
	if ($samplename eq "") {
		$samplename = $k;
	}
	$seqlen = length ($fastahash->{$k});
	set_sequence($fastahash->{$k});
}

my $genbank_header = "$samplename [organism=$orgname][moltype=Genomic DNA][location=chloroplast][topology=Circular][gcode=11]";
open FASTA_FH, ">", "$outfile.fsa";
print FASTA_FH ">$genbank_header\n";
print FASTA_FH get_sequence() . "\n";
close FASTA_FH;

my $gene_hash = {};
foreach my $id (@$gene_index_array) {
	my $gene = shift $gene_array2;
	$gene_hash->{$id} = $gene;
}

my $dest_gene_hash = {};
foreach my $id (@$destination_gene_index_array) {
	my $gene = shift $destination_gene_hash_array;
	$dest_gene_hash->{$id} = $gene;
}

# fill in the genes from the regionfile with the info from the destination gene array
my @final_gene_array = ();
foreach my $id (@$gene_index_array) {
	my $gene = $gene_hash->{$id};
	my $dest_gene = $dest_gene_hash->{$id};

	foreach my $q (keys %{$dest_gene->{"qualifiers"}}) {
		if (!(exists $gene->{"qualifiers"}->{$q})) {
			$gene->{"qualifiers"}->{$q} = $dest_gene->{"qualifiers"}->{$q}
		}
	}
	my @new_contains = ();
	$gene->{"id"} = $id;
	foreach my $destcontains (@{$dest_gene->{"contains"}}) {
		my $genecontains = shift $gene->{"contains"};
		$destcontains->{"region"} = $genecontains->{"region"};
		push @new_contains, $destcontains;
	}
	$gene->{"contains"} = \@new_contains;
	push @final_gene_array, $gene;
}

# print header
open FH, ">", "$outfile.tbl";
print FH ">Features\t$genbank_header\n";

# start printing genes
foreach my $gene (@final_gene_array) {
	# first, print overall gene information
	my $genename = $gene->{"qualifiers"}->{"gene"};
	my $gene_id = $gene->{"id"};
	my $feat_id = 0;
	foreach my $r (@{$gene->{'region'}}) {
		$r =~ /(\d+)\.\.(\d+)/;
		print FH "$1\t$2\tgene\n";
	}
	foreach my $q (keys %{$gene->{'qualifiers'}}) {
		print FH "\t\t\t$q\t$gene->{qualifiers}->{$q}\n";
	}

	# then, print each feature contained.
	foreach my $feat (@{$gene->{'contains'}}) {
		foreach my $reg (@{$feat->{"region"}}) {
			my $strand = "+";
			my ($start, $end) = split (/\.\./, $reg);
			if ($end < $start) {
				$strand = "-";
				my $oldend = $start;
				$start = $end;
				$end = $oldend;
			}
			my $regseq = sequence_for_interval ($reg);
			my $featname = $feat->{"type"};
			$feat_id++;
		}
		print FH sequin_feature ($feat->{'region'}, $feat);
	}
}
close FH;



__END__

=head1 NAME

parse_blast

=head1 SYNOPSIS

parse_blast [-blast blast_file] [-outputfile output_file]

=head1 OPTIONS

  -blast:           "outfmt 5" formatted blastn output
  -outputfile:      name of output file

=head1 DESCRIPTION


=cut
