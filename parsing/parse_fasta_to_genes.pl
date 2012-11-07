use Bio::SeqIO;
use Bio::Align::Utilities qw(cat);
use Pod::Usage;
use File::Basename;
require "subfuncs.pl";

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my ($fa_file, $gb_file, $out_file, $help) = 0;
GetOptions ('fasta=s' => \$fa_file,
            'genbank|gb_file=s' => \$gb_file,
            'outputfile=s' => \$out_file,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

my $whole_aln = make_aln_from_fasta_file ($fa_file);
my $gb_seqio = Bio::SeqIO->new(-file => $gb_file);
my @gene_alns;

while (my $seq_object = $gb_seqio->next_seq) {
	for my $feat_object ($seq_object->get_SeqFeatures) {
		if ($feat_object->primary_tag eq "gene") {
			my $name = main_name_for_gb_feature($feat_object);
			my @locations = $feat_object->location->each_Location;
			my $cat_aln = 0;
			my $strand = 0;
			foreach $loc (@locations) {
				$strand = $loc->strand;
				my $start = $loc->start;
				my $end = $loc->end;
				my $curr_slice = $whole_aln->slice($start, $end);
				if ($cat_aln == 0) {
					$cat_aln = $curr_slice;
				} else {
					$cat_aln = cat($cat_aln, $curr_slice);
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

			$cat_aln = $cat_aln->slice(1, $cat_aln->length()-3);
			$cat_aln->description($name);
			push @gene_alns, $cat_aln;
		}
	}
}

open my $gene_file, ">$out_file.fasta";

foreach my $aln (@gene_alns) {
	my $gene_name = $aln->description();
	foreach my $seq ($aln->each_seq()) {
		my $name = $seq->id() . "_$gene_name: " . $seq->length();
		print $gene_file ">$name\n";
		print $gene_file $seq->seq() . "\n";
	}
}
close $gene_file;

__END__

=head1 NAME

parse_fasta_to_genes

=head1 SYNOPSIS

parse_fasta_to_genes [-fasta fa_file] [-genbank gb_file] [-outputfile output_file]

=head1 OPTIONS

  -fasta:           fasta file of aligned sequences
  -genbank|gb_file: genbank file with CDS coordinates
  -outputfile:      prefix of output files

=head1 DESCRIPTION

Given a fasta file of aligned sequences and a corresponding genbank file
with the CDS coordinates, will create a fasta file with each
CDS corresponding to a separate sequence block.

=cut

