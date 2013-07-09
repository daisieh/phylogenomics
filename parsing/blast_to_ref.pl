use Bio::LocatableSeq;
use File::Temp qw/ tempfile tempdir /;

#take in a reference fasta file and an aligned sequence file

my $ref_file = shift @ARGV;
my $align_file = shift @ARGV;
my $out_file = shift @ARGV;


open refFH, "<", $ref_file;
my @reflines = <refFH>;
close refFH;
my $refid = shift @reflines;
$refid =~ />(.+)/;
$refid = $1;
my $refseq = "";
foreach my $line (@reflines) {
	chomp $line;
	$refseq .= $line;
}

#blast the seq
my ($fh, $blastfile) = tempfile(UNLINK => 1);
system ("blastn -query $align_file -subject $ref_file > $blastfile");
open BLAST_FH, "<", $blastfile;
my @lines = <BLAST_FH>;
close BLAST_FH;

foreach my $line (@lines) {
	if ($line =~ /Strand=Plus\/Minus/) {
		my ($fh, $tempreffile) = tempfile(UNLINK => 1);
		open refFH, ">", $tempreffile;
		my $ref_seq = new Bio::LocatableSeq(-seq => $refseq);
		$refseq = $ref_seq->revcom()->seq;
		print refFH ">$refid\n$refseq\n";
		system ("blastn -query $align_file -subject $tempreffile > $blastfile");
		last;
	}
}

open outFH, ">", $out_file;

print outFH ">".$refid . "\n" . $refseq;

open BLAST_FH, "<", $blastfile;
my @lines = <BLAST_FH>;
close BLAST_FH;

my ($query_start, $query_end, $subject_start, $subject_end) = 0;
my $query_seq = "";
my $query_id = "";
my ($curr_query_start, $curr_query_seq, $curr_query_end) = 0;

my $line = shift @lines;
while ($line) {
	if ($line =~ /Query=\s+(.*)$/) {
		print outFH "$query_seq\n>$1\n";
		$query_id = $1;
		($query_start, $query_end, $subject_start, $subject_end) = 0;
		$query_seq = "";
	} elsif ($line =~ /Query\s+(\d+)\s+(.+?)\s+(\d+).*$/) {
		$curr_query_start = $1;
		$curr_query_seq = $2;
		$curr_query_end = $3;
		if ($curr_query_start < $query_end) {
			while ($line !~ /Query=\s+(.*)$/) {
				$line = shift @lines;
			}
 			unshift @lines, $line;
		}
	} elsif ($line =~ /Sbjct\s+(\d+)\s+(.+?)\s+(\d+).*$/) {
		my $curr_subject_start = $1;
		my $curr_subject_seq = $2;
		my $curr_subject_end = $3;
		while ($curr_subject_seq =~ /(.*)-(.*)/) {
			$curr_subject_seq = "$1$2";
			$curr_query_seq =~ /(.{length($1)}).(.*)/;
			$curr_query_seq = "$1$2";
		}
		if ($query_start == 0) {
			$query_start = $curr_query_start;
			$query_end = $curr_query_end;
			$subject_end = $curr_subject_end;
			$query_seq = "n" x ($curr_subject_start - 1);
			$query_seq .= uc($curr_query_seq);
		} elsif (($query_end + 1) == $curr_query_start) {
			$query_end = $curr_query_end;
			$subject_end = $curr_subject_end;
			$query_seq .= uc($curr_query_seq);
		} else {
			$query_seq .= "n" x ($curr_subject_start - $subject_end - 1);
			$query_seq .= uc($curr_query_seq);
			$query_end = $curr_query_end;
			$subject_end = $curr_subject_end;
		}
	}
	$line = shift @lines;
}
print outFH "$query_seq\n";

close outFH;
