use strict;
use Bio::LocatableSeq;
use File::Temp qw/ tempfile tempdir /;

#take in a reference fasta file and an aligned sequence file

my $ref_file = shift @ARGV;
my $align_file = shift @ARGV;
my $out_file = shift @ARGV;

my @references = ();

open refFH, "<", $ref_file;
my $line = readline refFH;
$line =~ />(.+)/;
my $refid = $1;
my $refseq = "";

while ($line = readline refFH) {
	if ($line !~ />(.+)/) {
		chomp $line;
		$refseq .= $line;
	} else {
		push @references, "$refid\t$refseq";
		$line =~ />(.+)/;
		$refid = $1;
		$refseq = "";
	}
}
push @references, "$refid\t$refseq";

close refFH;

open outFH, ">", $out_file;
foreach my $ref (@references) {
	my ($refid, $refseq) = split (/\t/, $ref);
	print outFH "###\n>$refid\n$refseq";
	print outFH blast_to_ref($align_file, $refseq);
}
close outFH;

sub blast_to_ref {
	my $align_file = shift;
	my $refseq = shift;

	my ($fh, $tempreffile) = tempfile(UNLINK => 1);
	open refFH, ">", $tempreffile;
	my $ref_seq = new Bio::LocatableSeq(-seq => $refseq);
	print refFH ">$refid\n$refseq\n";
	close refFH;

	#blast the seq
	my ($fh, $blastfile) = tempfile(UNLINK => 1);
	system ("blastn -query $align_file -subject $tempreffile > $blastfile");
	open BLAST_FH, "<", $blastfile;
	my @lines = <BLAST_FH>;
	close BLAST_FH;

	foreach my $line (@lines) {
		if ($line =~ /Strand=Plus\/Minus/) {
			$refseq = $ref_seq->revcom()->seq;
			open refFH, ">", $tempreffile;
			print refFH ">$refid\n$refseq\n";
			close refFH;
			system ("blastn -query $align_file -subject $tempreffile > $blastfile");
			last;
		}
	}

	my $result = "";

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
			$result .= "$query_seq\n>$1\n";
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
	$result .= "$query_seq\n";
	return $result;
}
