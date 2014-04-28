use Getopt::Long;
use Pod::Usage;

my $help = 0;
my $blastfile = "";
my $outfile = "";


GetOptions ('blastfile=s' => \$blastfile,
            'outputfile=s' => \$outfile,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

if ($blastfile eq "") {
    pod2usage(-verbose => 1);
}

if ($outfile eq "") {
    pod2usage(-verbose => 1);
}

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

__END__

=head1 NAME

parse_blast

=head1 SYNOPSIS

First run: blastn -query comparison.fasta -subject reference.fasta -outfmt 3 -out blast_file
Then run: parse_blast [-blast blast_file] [-outputfile output_file]

=head1 OPTIONS

  -blast:           "outfmt 3" formatted blastn output
  -outputfile:      name of output file

=head1 DESCRIPTION

Parses an "outfmt 3" formatted blastn file to generate a list of regions to be used in
Genbank annotations.

=cut
