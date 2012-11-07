use Bio::SeqIO;
use Bio::Align::Utilities qw(cat);

require "subfuncs.pl";

# Bring in the file and format, or die with a nice
# usage statement if one or both arguments are missing.
my $usage  = "test_biopl.pl fasta_file intergenic_file result_file\n";
my $fa_file = shift or die $usage;
my $intergenic_file = shift or die $usage;
my $result_file = shift or die $usage;

my $whole_aln = make_aln_from_fasta_file ($fa_file);
my $cat_aln = 0;
my $result_str = "";
my $partition_str = "begin mrbayes;\n";
my $num_partitions = 0;
my $partition_list = "";
open my $F, "<$intergenic_file" or die "couldn't open fasta file";

my $curr_seg = readline $F;
while ($curr_seg) {
	$curr_seg =~ /(.+?)\t(.+?)\t(.+?)$/;
	my $name = $1;
	my $start = $2;
	my $end = $3;
	my $curr_slice = $whole_aln->slice($start, $end);

	if ($cat_aln == 0) {
		$start = 1;
		$cat_aln = $curr_slice;
	} else {
		$start = $cat_aln->length() + 1;
		$cat_aln = cat($cat_aln, $curr_slice);
	}
	$end = $start + $curr_slice->length() - 1;
	$partition_list .= "$name, ";
	$num_partitions++;
	$partition_str .= "\tCHARSET $name = $start-$end;\n";
	$curr_seg = readline $F;
}
$partition_list =~ s/, $/;/;
$partition_str .= "\tpartition spacers = $num_partitions: $partition_list\nend;\n";

my $result = convert_aln_to_nexus ($cat_aln);
open my $gene_file, ">$result_file.nex";
truncate $gene_file, 0;
print $gene_file $result;
print $gene_file $partition_str;

close $gene_file;
