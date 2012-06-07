use Bio::SeqIO;
use Bio::AlignIO;
require "subfuncs.pl";

# Bring in the file and format, or die with a nice
# usage statement if one or both arguments are missing.
my $usage  = "sliding_window.pl fasta_file window_size\n";
my $fa_file = shift or die $usage;
my $window_size = shift or die $usage;

my $start_pos = 1;
my $stop_pos = $window_size;
my $len;
my $flag = 1;
my $curr_aln;
$curr_aln = make_aln_from_fasta_file($fa_file);

my $sub_aln = $curr_aln->slice($start_pos, $stop_pos);

my $result = return_nexus_string ($sub_aln);
print "$result\n";

while ($flag) {
	my $gene_name = "$start_pos"."_"."$stop_pos";
	my $result = return_nexus_string ($curr_aln);

	$flag = perc_diff_partition ($curr_aln, $start_pos, $stop_pos);
	if ($flag > 0) {
		print "$start_pos\t$stop_pos\t$flag\n";
	}
	$start_pos += $window_size;
	$stop_pos += $window_size;
}
