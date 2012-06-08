use Bio::SeqIO;
require "subfuncs.pl";

# Bring in the file and format, or die with a nice
# usage statement if one or both arguments are missing.
my $usage  = "test_biopl.pl genbank_file fasta_file result_dir\n";
my $gb_file   = shift or die $usage;
my $fa_file = shift or die $usage;
my $result_dir = shift or die $usage;

my $start_pos = 10;
my $stop_pos = 50;
my $inseq;
$inseq = Bio::SeqIO->new(-file => "<$gb_file", -format => "genbank") or die "not genbank\n";
eval {$seq = $inseq->next_seq;} or die "1 not fasta\n";
my @seq_features = $seq->get_SeqFeatures;
foreach my $feature (@seq_features) {
	if ($feature->has_tag("protein_id")) {
		my @feat_tags;
		my $flag = 1;
		eval {@feat_tags = $feature->get_tag_values("gene")} or {$flag = 0};
		if ($flag) {
			my $gene_name = @feat_tags[0];
			$start_pos = $feature->start;
			$stop_pos = $feature->end;
			print "writing $gene_name, $start_pos-$stop_pos...\n";
			my $outfile = "$result_dir\/$gene_name.nex";
			my $curr_aln = make_aln_from_fasta_file($fa_file);
			my $sub_aln = $curr_aln->slice($start_pos, $stop_pos);
			my $result = convert_aln_to_nexus ($sub_aln);
			open my $gene_file, ">$outfile";
			truncate $gene_file, 0;
			print $gene_file $result;
			close $gene_file;
		}
	}
}

