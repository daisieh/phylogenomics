
my $blastfile = shift;

open FH, "<", $blastfile;

my $subject = "";
my $subject_length = 0;
my $line = "";
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
		# reading a subject
		while ($line !~ /Effective/) {
			$line = readline FH;
			if ($line =~ /Length=(\d+)/) {
				$subject_length = $1;
			}
			if ($line =~ /Query_\d+\s+(\d+).+?(\d+)$/) {
				if ($query_start == 0) {
					$query_start = $1;
				}
				$query_end = $2;
			}
			if ($line =~ /Subject_\d+\s+(\d+).+?(\d+)/) {
				if ($subject_start == 0) {
					$subject_start = $1;
				}
				$subject_end = $2;
			}
		}
		print "$subject\t$query_start\t$query_end";
		if (($subject_end - $subject_start) != ($subject_length - 1)) {
			print "\t>>>>>problem: ($subject_end - $subject_start) != ($subject_length - 1)<<<<<<";
		}
		print "\n";
	}
}
close FH;
