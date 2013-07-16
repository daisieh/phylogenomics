use strict;
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
my @refids = ();
open refFH, "<", $ref_file;
my $refseq = "";
my $refid = "";

while (my $line = readline refFH) {
	if ($line =~ />(.+)/) {
		$refid = $1;
		push @refids, $refid;
		if ($refseq ne "") {
			push @references, "$refseq";
			$refseq = "";
		}
	} else {
 		chomp $line;
		$refseq .= $line;
	}
}

push @references, "$refseq";

close refFH;

my (undef, $taxa_names) = parse_fasta($align_file);

my %result_matrices = ();
my $result_matrix;

for (my $i=0;$i<@refids;$i++) {
	print "blasting $refids[$i]...\n";
	my $this_blast_file = "";
	if ($blast_file ne "") {
		$this_blast_file = "$blast_file.$refids[$i]";
	}
	$result_matrix = blast_to_ref($align_file, $references[$i], $this_blast_file, $evalue);
	$result_matrix->{'reference'} = $references[$i];
	$result_matrices{$refids[$i]} = $result_matrix;
}

$refid = join ("+", @refids);
foreach my $key (keys %result_matrices) {
	$result_matrices{$key}->{$refid} = delete ($result_matrices{$key}->{'reference'});
	print "result_matrix $key has " . keys (%{$result_matrices{$key}}) . " keys\n";
}

my ($res1, $res2) = meld_matrices(\@refids, \%result_matrices);
my %mastertaxa = %{$res1};
my %regiontable = %{$res2};

open (fileOUT, ">", $out_file);

# debug code:
# my $x = 0;
# foreach my $key (@refids) {
# 	foreach my $keyid (keys %{$result_matrices{$key}}) {
# 		print fileOUT ">$keyid\n$result_matrices{$key}->{$keyid}\n";
# 	}
# 	print fileOUT $x++ . "###\n";
# }

my $value = $mastertaxa{$refid};
#
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

	my $revcomp = "";

	foreach my $line (@lines) {
		if ($line =~ /Strand=Plus\/Plus/) {
			last;
		} elsif ($line =~ /Strand=Plus\/Minus/) {
			$revcomp = reverse_complement ($refseq);
			open refFH, ">", $tempreffile;
			print refFH ">$refid\n$revcomp\n";
			close refFH;
			system ("blastn $blast_task -query $align_file -subject $tempreffile > $blastfile");
			last;
		}
	}

	my @seqids = ();
	my @seqs = ();
	open BLAST_FH, "<", $blastfile;
	my @lines = <BLAST_FH>;
	close BLAST_FH;

	my %result_matrix = ();
	my $curr_query_id = "";
	my $query_seq = "";
	my ($query_end, $subject_end) = 0;
	my ($curr_query_start, $curr_query_seq, $curr_query_end) = 0;

	if ($blast_task =~ /-task blastn-short/) {
		while (my $line = shift @lines) {
			if ($line =~ /Query=\s+(.*)$/) {
				if ($query_seq ne "") {
					$query_seq .= "n" x (length($refseq) - length($query_seq));
					$result_matrix{$curr_query_id} = $query_seq;
				}
				$curr_query_id = $1;
				push @seqids, $curr_query_id;
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
					$result_matrix{$curr_query_id} = $query_seq;
				}
				$curr_query_id = $1;
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
	$result_matrix{$curr_query_id} = $query_seq;

	# flip the seqs if we did the blast for the revcomp
	if ($revcomp ne "") {
		foreach my $id (keys %result_matrix) {
			$result_matrix{$id} = reverse_complement ($result_matrix{$id});
		}
	}
	return \%result_matrix;
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

