#parse_fasta_to_genes.pl [-fasta fa_file] [-genbank gb_file] [-outputfile output_file] [-multiple]

#makeblastdb -in ref_gene_set -out ref_db -dbtype 'nucl'
#blastn -query new_file -subject ref_db -outfmt 3

my $blastfile = shift;
my $outfile = shift;

open FH, "<", $blastfile;
open OUT_FH, ">", $outfile;

my $subject = "";
my $subject_length = 0;
my $line = "";
my $curr_pos = 0;
while (defined $line) {
	$line = readline FH;
	if ($line =~ /^\s+$/) { next; }

	# start of new subject
	if ($line =~ /Subject= (.*)$/) {
		my $query_start = 0;
		my $query_end = 0;
		my $subject_start = 0;
		my $subject_end = 0;
		$subject = $1;
		print "start subject $subject\n";
		# reading a subject
		while ($line !~ /Effective/) {
			$line = readline FH;
			if ($line =~ /Length=(\d+)/) {
				$subject_length = $1;
			}
			if ($line =~ /Query_\d+\s+(\d+).+?(\d+)$/) {
				$query_end = $2;
				if ($query_start == 0) {
					$query_start = $1;
				}
				$line = readline FH;
				if ($line =~ /Subject_\d+\s+(\d+)\s+(\S+)\s+(\d+)/) {
					$subject_end = $3;
					my $subj = $2;
					if ($subject_start == 0) {
						$subject_start = $1;
						if (($subject_end - $subject_start) != ($query_end - $query_start)) {
							$query_start = $query_end - length($subj)+1;
						}
					}
				}
			}
		}
		select OUT_FH;
		$| = 1;
		print OUT_FH "$subject\t$query_start\t";
		if (($query_end - $query_start) != ($subject_length - 1)) {
			if (($subject_end - $subject_start) == -($subject_length - 1)) {
				print OUT_FH "$query_end\treverse strand";
			} elsif (($subject_end==0) && ($subject_start == 0)) {
				print OUT_FH "$query_end\tno match";
			} elsif ($subject_end==$subject_start) {
				print OUT_FH "$query_end\tcaught twice";
			} else {
				$query_end = $query_start + $subject_length - 1;
				print OUT_FH "$query_end";
			}
		} else {
			print OUT_FH "$query_end";
		}
		print OUT_FH "\n";
		select STDOUT;
	}
}
close FH;
close OUT_FH;
