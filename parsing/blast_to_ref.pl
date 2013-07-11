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
my $blast_file = "";
my $evalue = 10;

GetOptions ('fasta|input=s' => \$align_file,
            'outputfile=s' => \$out_file,
            'reference=s' => \$ref_file,
            'blastfile=s' => \$blast_file,
            'evalue=f' => \$evalue,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

unless($ref_file and $align_file and $out_file) {
    pod2usage(-msg => "Must specify reference file, alignment file, and output file.");
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
		push @references, "$refseq";
		$line =~ />(.+)/;
		$refid .= "+$1";
		$refseq = "";
	}
}
push @references, "$refseq";

close refFH;

my @result_files = ();
foreach my $refseq (@references) {
	my ($fh, $ref_result_file) = tempfile(UNLINK => 1, SUFFIX => ".fasta");
	print $fh ">$refid\n$refseq";
	print $fh blast_to_ref($align_file, $refseq, $blast_file, $evalue);
	close $fh;
	push @result_files, $ref_result_file;
}

my ($res1, $res2) = meld_matrices(\@result_files);
my %mastertaxa = %{$res1};
my %regiontable = %{$res2};
open (fileOUT, ">", $out_file);

my $value = $mastertaxa{$refid};

print fileOUT ">$refid\n$value\n";
delete $mastertaxa{$refid};
delete $mastertaxa{"length"};


# currently printing in fasta format: perhaps add a flag to alter this?
foreach my $key ( keys %mastertaxa ) {
	my $value = $mastertaxa{$key};
	print fileOUT ">$key\n$value\n";
}
close fileOUT;



# SUBFUNCTIONS START

sub blast_to_ref {
	my $align_file = shift;
	my $refseq = shift;
	my $blastfile = shift;
	my $evalue = shift;

	my $blast_task = "-evalue $evalue ";
	if (length ($refseq) < 50) {
		$blast_task .= "-task blastn-short";
	}

	my (undef, $tempreffile) = tempfile(UNLINK => 1);
	open refFH, ">", $tempreffile;
	my $ref_seq = new Bio::LocatableSeq(-seq => $refseq);
	print refFH ">reference\n$refseq\n";
	close refFH;

	#blast the seq
	if ($blastfile eq "") {
		(undef, $blastfile) = tempfile(UNLINK => 1);
	}
	system ("blastn $blast_task -query $align_file -subject $tempreffile > $blastfile");
	open BLAST_FH, "<", $blastfile;
	my @lines = <BLAST_FH>;
	close BLAST_FH;

	foreach my $line (@lines) {
		if ($line =~ /Strand=Plus\/Plus/) {
			last;
		} elsif ($line =~ /Strand=Plus\/Minus/) {
			$refseq = $ref_seq->revcom()->seq;
			open refFH, ">", $tempreffile;
			print refFH ">$refid\n$refseq\n";
			close refFH;
			system ("blastn $blast_task -query $align_file -subject $tempreffile > $blastfile");
			last;
		}
	}

	open BLAST_FH, "<", $blastfile;
	my @lines = <BLAST_FH>;
	close BLAST_FH;

	my $result = "";
	my $query_seq = "";
	my ($query_end, $subject_end) = 0;
	my ($curr_query_start, $curr_query_seq, $curr_query_end) = 0;

	if ($blast_task =~ /-task blastn-short/) {
		while (my $line = shift @lines) {
			if ($line =~ /Query=\s+(.*)$/) {
				if ($query_seq ne "") {
					$query_seq .= "n" x (length($refseq) - length($query_seq));
				}
				$result .= "$query_seq\n>$1\n";
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
				if ($query_seq eq "") {
					$query_seq = "n" x ($curr_subject_start - 1);
					$query_seq .= uc($curr_query_seq);
				}
			}
		}
	} else {
		while (my $line = shift @lines) {
			if ($line =~ /Query=\s+(.*)$/) {
				if ($query_seq ne "") {
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
				if ($query_seq eq "") {
					$query_seq = "n" x ($curr_subject_start - 1);
					$query_seq .= uc($curr_query_seq);
				} elsif (($query_end + 1) == $curr_query_start) {
					$query_seq .= uc($curr_query_seq);
				} else {
					$query_seq .= "n" x ($curr_subject_start - $subject_end - 1);
					$query_seq .= uc($curr_query_seq);
				}
				$subject_end = $curr_subject_end;
				$query_end = $curr_query_end;
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
  -blastfile:	    optional: if specified, put BLASTN output in this file.
  -evalue:	        optional: if specified, sets evalue threshold for BLASTN hits.

=head1 DESCRIPTION

Takes a fasta file and finds aligned regions in each sequence in the fasta file that
match the reference sequence(es). Returns a fasta file of aligned regions of similarity.
Uses BLASTN to find regions of similarity.

=cut

