use Bio::SeqIO;
use Bio::Align::Utilities qw(cat);
use Pod::Usage;
use File::Basename;
use Getopt::Long;

require "subfuncs.pl";

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my ($fa_file, $gb_file, $out_file, $help, $multiple) = 0;
GetOptions ('fasta=s' => \$fa_file,
            'genbank|gb_file=s' => \$gb_file,
            'outputfile=s' => \$out_file,
            'help|?' => \$help,
            'multiple!' => \$multiple) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

#my $inseq = Bio::SeqIO->new(-file => "<$fa_file", -format => "fasta");
my $whole_aln = make_aln_from_fasta_file ($fa_file);
my $gb_seqio = Bio::SeqIO->new(-file => $gb_file);
my @gene_seqs;
my @gene_alns;
my @gene_locs;

while (my $seq_object = $gb_seqio->next_seq) {
	for my $feat_object ($seq_object->get_SeqFeatures) {
		if ($feat_object->primary_tag eq "CDS") {
			my $name = main_name_for_gb_feature($feat_object);
			my @locations = $feat_object->location->each_Location;
			my $cat_aln = 0;
			my $cat_genome_seq = 0;
			my $cat_loc = "";
			my $strand = 0;
			foreach $loc (@locations) {
				$strand = $loc->strand;
				my $start = $loc->start;
				my $end = $loc->end;
				my $curr_slice = $whole_aln->slice($start, $end);
				if ($cat_aln == 0) {
					$cat_aln = $curr_slice;
					$cat_genome_seq = $curr_slice->get_seq_by_pos(1)->seq();
					$cat_loc = "$start..$end,";
				} else {
					$cat_aln = cat($cat_aln, $curr_slice);
					$cat_genome_seq = $cat_aln->get_seq_by_pos(1)->seq();
					$cat_loc .= "$start..$end,"
				}
			}
			if ($strand < 0) {
				# must flip each seq in the curr_slice
				my $flipped_aln = Bio::SimpleAlign->new();
				foreach $seq ( $cat_aln->each_seq() ) {
					$seq = $seq->revcom();
					$flipped_aln->add_seq($seq);
				}
				$cat_aln = $flipped_aln;
			}

# 			$cat_aln = $cat_aln->slice(1, $cat_aln->length()-3);
			$cat_aln->description($name);
			push @gene_alns, $cat_aln;
			push @gene_seqs, $cat_genome_seq;
			push @gene_locs, $cat_loc;
		}
	}
}

open OUT_FH, ">", $out_file;

for (my $j=0;$j<@gene_alns;$j++) {
	my $seq = @gene_alns[$j]->get_seq_by_pos(1);
	my $genomeseq = uc(@gene_seqs[$j]);
 	my $gene_name = @gene_alns[$j]->description();
	my $working_seq = $genomeseq;
	my $i=0;
 	print "$gene_name\t" . @gene_locs[$j] . "\n";
	while (length($working_seq) > 50) {
		$working_seq =~ /(.{10})(.{10})(.{10})(.{10})(.{10})(.*)/;
		print OUT_FH ($i*50)."\t$gene_name\tseq\t" . "$1 $2 $3 $4 $5" . "\n";
		$working_seq = $6;
		$i++;
	}
	print OUT_FH ($i*50)."\t$gene_name\tseq\t$working_seq\n";
	my $pseq = $seq->translate();
	$working_seq = $pseq->seq;
	$i=0;
	while (length($working_seq) > 50) {
		$working_seq =~ /(.{10})(.{10})(.{10})(.{10})(.{10})(.*)/;
		print OUT_FH ($i*50)."\t$gene_name\tprot\t" . "$1 $2 $3 $4 $5". "\n";
		$working_seq = $6;
		$i++;
	}
	print OUT_FH ($i*50)."\t$gene_name\tprot\t$working_seq\n";

# 	print OUT_FH ">$gene_name\tprot\n" . $pseq->seq . "\n";
# 	print OUT_FH ">$gene_name\tgenomeseq\n" . $genomeseq . "\n";
}

# for (my $j=0;$j<@gene_alns; $j++) {
# 	my $gene_name = $gene_alns[$j]->description();
# 	foreach my $seq ($gene_alns[$j]->each_seq()) {
# 		my $name = $seq->id() . "_$gene_name: " . $seq->length();
# 		print OUT_FH ">$name\n";
# 		print OUT_FH $seq->seq() . "\n";
# 		print OUT_FH ">$name prot\n";
# 		print OUT_FH $gene_prots[$j] . "\n";
# 	}
# }
#
#
close OUT_FH;
__END__


=head1 NAME

parse_fasta_to_genes

=head1 SYNOPSIS

parse_fasta_to_genes [-fasta fa_file] [-genbank gb_file] [-outputfile output_file] [-multiple]

=head1 OPTIONS

  -fasta:           fasta file of aligned sequences
  -genbank|gb_file: genbank file with CDS coordinates
  -outputfile:      prefix of output files
  -multiple:		output separate fasta files for each CDS

=head1 DESCRIPTION

Given a fasta file of aligned sequences and a corresponding genbank file
with the CDS coordinates, will create a fasta file with each
CDS corresponding to a separate sequence block.

=cut

