use strict;
use Bio::LocatableSeq;
use File::Temp qw/ tempfile tempdir /;
use Getopt::Long;
use Pod::Usage;
use File::Basename;

require "subfuncs.pl";

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my $runline = "running " . basename($0) . " " . join (" ", @ARGV) . "\n";

my $ref_file = 0;
my $align_file = 0;
my $out_file = 0;
my $help = 0;

GetOptions ('fasta|input=s' => \$align_file,
            'outputfile=s' => \$out_file,
            'reference=s' => \$ref_file,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

unless($ref_file and $align_file and $out_file) {
	print "hi\n";
    pod2usage(-msg => "Must specify all options.");
}

print $runline;

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

my @result_files = ();
foreach my $ref (@references) {
	my ($fh, $tempreffile) = tempfile(UNLINK => 1);
	my ($refid, $refseq) = split (/\t/, $ref);
	print $fh ">$refid\n$refseq";
	if (length($refseq) < 50) {
		print $fh blast_to_short_ref($align_file, $refseq);
	} else {
		print $fh blast_to_ref($align_file, $refseq);
	}
	close $fh;
	push @result_files, $tempreffile;
}

open outFH, ">", $out_file;
foreach my $resultfile (@result_files) {
	print outFH "###\n";
	open resFH, "<", $resultfile;
	my @lines = <resFH>;
	print outFH join("",@lines);
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

	my ($query_end, $subject_end) = 0;
	my $query_seq = "";
	my ($curr_query_start, $curr_query_seq, $curr_query_end) = 0;

	while (my $line = shift @lines) {
		if ($line =~ /Query=\s+(.*)$/) {
			if ($query_end> 0) {
				$query_seq .= "n" x (length($refseq) - length($query_seq));
			}
			$result .= "$query_seq\n>$1\n";
			($query_end, $subject_end) = 0;
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
			if ($query_end == 0) {
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
	}
	$query_seq .= "n" x (length($refseq) - length($query_seq));
	$result .= "$query_seq\n";
	return $result;
}

sub blast_to_short_ref {
	my $align_file = shift;
	my $refseq = shift;
	my ($fh, $tempreffile) = tempfile(UNLINK => 1);
	open refFH, ">", $tempreffile;
	my $ref_seq = new Bio::LocatableSeq(-seq => $refseq);
	print refFH ">$refid\n$refseq\n";
	close refFH;

	#blast the seq
	my ($fh, $blastfile) = tempfile(UNLINK => 1);
	system ("blastn -task blastn-short -query $align_file -subject $tempreffile > $blastfile");
	open BLAST_FH, "<", $blastfile;
	my @lines = <BLAST_FH>;
	close BLAST_FH;

	foreach my $line (@lines) {
		if ($line =~ /Strand=Plus/) {
			if ($line =~ /Plus\/Minus/) {
				$refseq = $ref_seq->revcom()->seq;
				open refFH, ">", $tempreffile;
				print refFH ">$refid\n$refseq\n";
				close refFH;
				system ("blastn -task blastn-short -query $align_file -subject $tempreffile > $blastfile");
			}
			last;
		}
	}

	my $result = "";

	open BLAST_FH, "<", $blastfile;
	my @lines = <BLAST_FH>;
	close BLAST_FH;

	my $query_end = 0;
	my $query_seq = "";
	my ($curr_query_seq, $curr_query_end) = 0;
	while (my $line = shift @lines) {
		if ($line =~ /Query=\s+(.*)$/) {
			if ($query_end> 0) {
				$query_seq .= "n" x (length($refseq) - length($query_seq));
			}
			$result .= "$query_seq\n>$1\n";
			$query_end = 0;
			$query_seq = "";
		} elsif ($line =~ /Query\s+(\d+)\s+(.+?)\s+(\d+).*$/) {
			$curr_query_seq = $2;
			$curr_query_end = $3;
		} elsif ($line =~ /Sbjct\s+(\d+)\s+(.+?)\s+(\d+).*$/) {
			my $curr_subject_start = $1;
			my $curr_subject_seq = $2;
			my $curr_subject_end = $3;
			while ($curr_subject_seq =~ /(.*)-(.*)/) {
				$curr_subject_seq = "$1$2";
				$curr_query_seq =~ /(.{length($1)}).(.*)/;
				$curr_query_seq = "$1$2";
			}
			if ($query_end == 0) {
				$query_end = $curr_query_end;
 				$query_seq = "n" x ($curr_subject_start - 1);
				$query_seq .= uc($curr_query_seq);
			}
		}
	}
	$query_seq .= "n" x (length($refseq) - length($query_seq));
	$result .= "$query_seq\n";
	return $result;
}


__END__

=head1 NAME

blast_to_ref

=head1 SYNOPSIS

blast_to_ref -fasta fastafile -reference reffile -output outputfile

=head1 OPTIONS

  -fasta|input:     fasta file of aligned sequences.
  -reference:       fasta file with sequences of interest.
  -outputfile:      output file name.

=head1 DESCRIPTION

Takes an aligned fasta file and finds aligned regions that match reference sequence(es).
Uses BLASTN to find regions of similarity.

=cut

