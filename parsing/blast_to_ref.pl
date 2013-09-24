use strict;
use File::Temp qw/ tempfile tempdir /;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use XML::Simple;

require "subfuncs.pl";

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my $runline = "running " . basename($0) . " " . join (" ", @ARGV) . "\n";

my $ref_file = 0;
my $align_file = 0;
my $out_file = 0;
my $separate = 0;
my $help = 0;
my $blast_file = "";
my $evalue = 10;
my $ref_out = 0;
my $no_blast = 0;
my $debug = 0;

GetOptions ('fasta|input=s' => \$align_file,
            'outputfile=s' => \$out_file,
            'reference=s' => \$ref_file,
            'separate' => \$separate,
            'blastfile=s' => \$blast_file,
            'evalue=f' => \$evalue,
            'include_ref' => \$ref_out,
            'no_blast|noblast' => \$no_blast,
            'debug' => \$debug,
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
	if ($line =~ />(.+)$/) {
		$refid = $1;
		$refid =~ s/\s/_/g;
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

my $result_matrices = ();

my (undef, $tempreffile) = tempfile(UNLINK => 1);
if ($blast_file eq "") {
	(undef, $blast_file) = tempfile(UNLINK => 1);
}

if ($no_blast != 1) {
	system ("makeblastdb -in $ref_file -dbtype nucl -out $tempreffile.db");
	system ("blastn -task blastn -query $align_file -db $tempreffile.db -outfmt 5 -out $blast_file");
} else { debug ("no blast\n"); }
$result_matrices = blast_to_ref("$blast_file");

for (my $i=0;$i<@refids;$i++) {
	$result_matrices->{$refids[$i]}->{'reference'} = $references[$i];
}

if ($separate == 0) {
	$refid = join ("|", @refids);
} else {
	$refid = "reference";
}

for (my $i=0;$i<@refids;$i++) {
	my $key = $refids[$i];
	my $reflen = length($references[$i]);
	$result_matrices->{$key}->{$refid} = delete ($result_matrices->{$key}->{'reference'});
	foreach my $k (keys (%{$result_matrices->{$key}})) {
		# pad out the sequence at the end so that they're aligned for matrixmeld.
		${$result_matrices->{$key}}{$k} .= "-" x ($reflen - length(${$result_matrices->{$key}}{$k}));
	}
}

if ($separate == 0) {
	my ($res1, $res2) = meld_matrices(\@refids, $result_matrices);
	my %mastertaxa = %{$res1};
	my %regiontable = %{$res2};

	my $value = $mastertaxa{$refid};
	open (fileOUT, ">", "$out_file.fasta");
	if ($ref_out == 1) {
		print fileOUT ">$refid\n$value\n";
	}
	delete $mastertaxa{$refid};
	delete $mastertaxa{"length"};

	# currently printing in fasta format: perhaps add a flag to alter this?
	foreach my $key ( keys %mastertaxa ) {
		my $value = $mastertaxa{$key};
		print fileOUT ">$key\n$value\n";
	}
	close fileOUT;
} else {
	foreach my $refname (@refids) {
		open (fileOUT, ">", "$out_file.$refname.fasta");
		my $value = ${$result_matrices->{$refname}}{$refid};
		if ($ref_out == 1) {
			print fileOUT ">$refid\n$value\n";
		}
		delete ${$result_matrices->{$refname}}{$refid};
		delete ${$result_matrices->{$refname}}{"length"};
		foreach my $key ( keys (%{$result_matrices->{$refname}}) ) {
			my $value = ${$result_matrices->{$refname}}{$key};
			print fileOUT ">$key\n$value\n";
		}
		close fileOUT;
	}
}


# SUBFUNCTIONS START
sub blast_to_ref {
	my $blast_xml = shift;
	my %result_matrix = ();

	my $parser = new XML::Simple();

	open XML_FH, "<", $blast_xml;
	my %xml_strings = ();

	my $xml = readline XML_FH;
	my $query_name = "";

	while (my $line = readline XML_FH) {
		if ($line =~ /<\?xml/) {
			# we've reached the next xml object
			$xml_strings{$query_name} = "$xml";
			$xml = "$line";
		} else {
			if ($line =~ /<BlastOutput_query-def>(.*?)<\/BlastOutput_query-def>/) {
				$query_name = $1;
			}
			$xml .= $line;
		}
	}
	$xml_strings{$query_name} = "$xml";

	close XML_FH;

	foreach my $key (keys %xml_strings) {
		my $tree = $parser->XMLin($xml_strings{$key}, ForceArray => 1);
		my $iterations = ($tree->{"BlastOutput_iterations"}[0]->{"Iteration"}); # key Iteration has a value that is an anonymous array of iteration hashes.
		foreach my $iteration (@$iterations) { # for each iteration
			my $iterid = $iteration->{"Iteration_query-ID"}[0];
			my $iterquerydef = $iteration->{"Iteration_query-def"}[0];
			my $hits = $iteration->{"Iteration_hits"}[0]->{"Hit"};
			foreach my $hit (@$hits) {
				my $hsps = $hit->{"Hit_hsps"}[0]->{"Hsp"}; # Key Hsp has a value that is an anonymous array of the hit hashes.
				my $hitdef = $hit->{"Hit_def"}[0];
				debug ("HIT $hitdef QUERY $iterquerydef\n");
				my @sorted_hsps = sort { $b->{"Hsp_score"}[0] - $a->{"Hsp_score"}[0] } @$hsps;
				my @selected_hsps = ();
				my ($query_start, $query_end, $hit_start, $hit_end, $align_strand) = 0;
				while (my $hsp = shift @sorted_hsps) {
					# if $hsp is on the minus strand, reverse-comp before dealing with it.
					if ($hsp->{"Hsp_hit-frame"}[0] < 0) {
						my $hit_to = $hsp->{"Hsp_hit-to"}[0];
						my $hit_from = $hsp->{"Hsp_hit-from"}[0];
						$hsp->{"Hsp_hit-to"}[0] = $hit_from;
						$hsp->{"Hsp_hit-from"}[0] = $hit_to;
						$hsp->{"Hsp_qseq"}[0] = lc(reverse_complement($hsp->{"Hsp_qseq"}[0]));
						$hsp->{"Hsp_hseq"}[0] = lc(reverse_complement($hsp->{"Hsp_hseq"}[0]));
					}
					my $aln_length = $hsp->{"Hsp_align-len"}[0];
					my $aln_percent = sprintf("%.2f",$hsp->{"Hsp_identity"}[0] / $aln_length);
# 					if ($hsp->{"Hsp_evalue"}[0] > 10) {
# 						# this is no good: move on to the next hit.
# 						next;
# 					}
					if (@selected_hsps == 0) {
						# this is the first chunk of sequence.
						$query_start = $hsp->{"Hsp_query-from"}[0];
						$query_end = $hsp->{"Hsp_query-to"}[0];
						$hit_start = $hsp->{"Hsp_hit-from"}[0];
						$hit_end = $hsp->{"Hsp_hit-to"}[0];
						$align_strand = $hsp->{"Hsp_hit-frame"}[0];
						debug ("\tfirst seq with score ".$hsp->{"Hsp_score"}[0].", strand $align_strand: ".$hsp->{"Hsp_hit-from"}[0] ."-".$hsp->{"Hsp_hit-to"}[0]."\t".$hsp->{"Hsp_query-from"}[0]."-".$hsp->{"Hsp_query-to"}[0]."\n");
						push @selected_hsps, $hsp;
					} else {
						# if it's not the first chunk, we need to check to see if it falls in a reasonable part of the sequence: both hit and query have to lie before or after the working region.
						debug ("looking for hsps that work with $align_strand: $hit_start-$hit_end\t$query_start-$query_end\n");
						debug ("\tlooking at next seq with score ".$hsp->{"Hsp_score"}[0].", strand ".$hsp->{"Hsp_hit-frame"}[0].": ".$hsp->{"Hsp_hit-from"}[0] ."-".$hsp->{"Hsp_hit-to"}[0]."\t".$hsp->{"Hsp_query-from"}[0]."-".$hsp->{"Hsp_query-to"}[0]."\n");
						if ($hsp->{"Hsp_hit-frame"}[0] != $align_strand) {
							debug ("\t\tNO: next seq on wrong strand\n");
						} else {
							if ($align_strand > 0) {
								if (($hsp->{"Hsp_hit-to"}[0] < $hit_start) && ($hsp->{"Hsp_query-to"}[0] < $query_start)) {
									# we know $hsp will fall completely on the left, so we're good.
									debug ("\t\tYES: next seq falls to the left of $hit_start and $query_start: ".$hsp->{"Hsp_hit-from"}[0] ."-".$hsp->{"Hsp_hit-to"}[0]."\t".$hsp->{"Hsp_query-from"}[0]."-".$hsp->{"Hsp_query-to"}[0]."\n");
									$query_start = $hsp->{"Hsp_query-from"}[0];
									$hit_start = $hsp->{"Hsp_hit-from"}[0];
									push @selected_hsps, $hsp;
								} elsif (($hsp->{"Hsp_hit-from"}[0] > $hit_end) && ($hsp->{"Hsp_query-from"}[0] > $query_end)) {
									# we know $hsp will fall completely on the left, so we're good.
									debug ("\t\tYES: next seq falls to the right of $hit_end and $query_end: ".$hsp->{"Hsp_hit-from"}[0] ."-".$hsp->{"Hsp_hit-to"}[0]."\t".$hsp->{"Hsp_query-from"}[0]."-".$hsp->{"Hsp_query-to"}[0]."\n");
									$query_end = $hsp->{"Hsp_query-to"}[0];
									$hit_end = $hsp->{"Hsp_hit-to"}[0];
									push @selected_hsps, $hsp;
								} else {
									debug ("\t\tNO: hit and query fell on opposite sides\n");
								}
							} else {
								# we're dealing with a minus strand:
								# hit should be the same, but query is backwards.
								if (($hsp->{"Hsp_hit-to"}[0] < $hit_start) && ($hsp->{"Hsp_query-from"}[0] > $query_end)) {
									# we know $hsp will fall completely on the left, so we're good.
									debug ("\t\tYES: next seq falls to the left of $hit_start and -$query_end: ".$hsp->{"Hsp_hit-from"}[0] ."-".$hsp->{"Hsp_hit-to"}[0]."\t".$hsp->{"Hsp_query-from"}[0]."-".$hsp->{"Hsp_query-to"}[0]."\n");
									$query_end = $hsp->{"Hsp_query-to"}[0];
									$hit_start = $hsp->{"Hsp_hit-from"}[0];
									push @selected_hsps, $hsp;
								} elsif (($hsp->{"Hsp_hit-from"}[0] > $hit_end) && ($hsp->{"Hsp_query-to"}[0] < $query_start)) {
									# we know $hsp will fall completely on the left, so we're good.
									debug ("\t\tYES: next seq falls to the right of $hit_end and -$query_start: ".$hsp->{"Hsp_hit-from"}[0] ."-".$hsp->{"Hsp_hit-to"}[0]."\t".$hsp->{"Hsp_query-from"}[0]."-".$hsp->{"Hsp_query-to"}[0]."\n");
									$query_start = $hsp->{"Hsp_query-from"}[0];
									$hit_end = $hsp->{"Hsp_hit-to"}[0];
									push @selected_hsps, $hsp;
								} else {
									debug ("\t\tNO: hit and query fell on opposite sides\n");
								}
							}
						}
					}

					# remove all hits that would fall within this range, starting from the far end of the array (this ensures we don't mess up the count when we splice things out).
					my $this_hsp = 0;
					debug ("\t\tremoving any hsps that fall in the current range: $hit_start-$hit_end\t$query_start-$query_end\n");
					for (my $j=@sorted_hsps-1; $j>=0; $j--) {
						$this_hsp = $sorted_hsps[$j];
						if (($this_hsp->{"Hsp_hit-from"}[0] >= $hit_start) && ($this_hsp->{"Hsp_hit-to"}[0] <= $hit_end)) {
							debug ("\t\t\tremoving hsp ".$this_hsp->{"Hsp_hit-from"}[0]."-".$this_hsp->{"Hsp_hit-to"}[0]." with score ". $this_hsp->{"Hsp_score"}[0] . " because it's inside the hit\n");
							splice (@sorted_hsps,$j,1);
						} elsif (($this_hsp->{"Hsp_query-from"}[0] >= $query_start) && ($this_hsp->{"Hsp_query-to"}[0] <= $query_end)) {
							debug ("\t\t\tremoving hsp ".$this_hsp->{"Hsp_query-from"}[0]."-".$this_hsp->{"Hsp_query-to"}[0]." with score ". $this_hsp->{"Hsp_score"}[0] . " because it's inside the query\n");
							splice (@sorted_hsps,$j,1);
						} else {
							debug ("\t\t\thsp ".$this_hsp->{"Hsp_hit-from"}[0]."-".$this_hsp->{"Hsp_hit-to"}[0]."\t".$this_hsp->{"Hsp_query-from"}[0]."-".$this_hsp->{"Hsp_query-to"}[0]." is not inside\n");
						}
					}
					debug ("\tthere are still ". @sorted_hsps . " hsps left\n");
				}
				debug ("final sequence:\n");
				my $hit_end = 0;
				my $sequence = "";
				my @selected_hsps = sort { $a->{"Hsp_hit-from"}[0] - $b->{"Hsp_hit-from"}[0] } @selected_hsps;
				foreach my $hsp (@selected_hsps) { # for each hsp for this query
					my $aln_length = $hsp->{"Hsp_align-len"}[0];
 					my $aln_percent = sprintf("%.2f",$hsp->{"Hsp_identity"}[0] / $aln_length);
					debug ("\t".$hsp->{"Hsp_hit-from"}[0]."-".$hsp->{"Hsp_hit-to"}[0]."\t".$hsp->{"Hsp_query-from"}[0]."-".$hsp->{"Hsp_query-to"}[0]."\t$aln_length\t$aln_percent\t".$hsp->{"Hsp_score"}[0]."\n");
					my $last_end = $hit_end;
					$hit_end = $hsp->{"Hsp_hit-to"}[0];
					# remove gaps from query seq that might have been inserted into the ref seq.
					while ($hsp->{"Hsp_hseq"}[0] =~ /^(.*?)(-+)(.*)$/) {
						my $left = length($1);
						my $gap = length($2);
						my $right = length($3);
						$hsp->{"Hsp_hseq"}[0] = $1 . $3;
						$hsp->{"Hsp_qseq"}[0] =~ /^(.{$left})(.{$gap})(.{$right})/;
						$hsp->{"Hsp_qseq"}[0] = $1 . $3;
					}
					$sequence .= "n" x ($hsp->{"Hsp_hit-from"}[0] - $last_end - 1) . $hsp->{"Hsp_qseq"}[0];
				}
				if ($sequence =~ /\w/) {
					$result_matrix{$hitdef}->{$iterquerydef} = $sequence;
				} else {
					debug ("\tno match\n");
				}
			}
		}
	}
	return \%result_matrix;
}

sub debug {
	my $str = shift;

	if ($debug == 1) {
		print STDERR "$str";
	}
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
  -include_ref:     optional: if specified, includes the reference sequence in the fasta output.
  -no_blast:        optional: (use with -blastfile flag) uses the specified blast xml file as input (for debugging)
  -debug:           optional: spews lots of information to STDERR.

=head1 DESCRIPTION

Takes a fasta file and finds aligned regions in each sequence in the fasta file that
match the reference sequence(es). Returns a fasta file of aligned regions of similarity.
Uses BLASTN to find regions of similarity.

=cut

