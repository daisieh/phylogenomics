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

system ("makeblastdb -in $ref_file -dbtype nucl -out $tempreffile.db");
system ("blastn -query $align_file -db $tempreffile.db -outfmt 5 -out $blast_file");
$result_matrices = blast_to_ref_xml("$blast_file");

for (my $i=0;$i<@refids;$i++) {
	$result_matrices->{$refids[$i]}->{'reference'} = $references[$i];
}

$refid = join ("+", @refids);
for (my $i=0;$i<@refids;$i++) {
	my $key = $refids[$i];
	my $reflen = length($references[$i]);
	$result_matrices->{$key}->{$refid} = delete ($result_matrices->{$key}->{'reference'});
	foreach my $k (keys (%{$result_matrices->{$key}})) {
		# pad out the sequence at the end so that they're aligned for matrixmeld.
		${$result_matrices->{$key}}{$k} .= "-" x ($reflen - length(${$result_matrices->{$key}}{$k}));
	}
}

my ($res1, $res2) = meld_matrices(\@refids, $result_matrices);
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
sub blast_to_ref_xml {
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
				# if any of the hits are on the minus strand, reverse-comp before dealing with them.
				foreach my $hsp(@$hsps) {
					my $hit_to = $hsp->{"Hsp_hit-to"}[0];
					my $hit_from = $hsp->{"Hsp_hit-from"}[0];
					my $query_to = $hsp->{"Hsp_query-to"}[0];
					my $query_from = $hsp->{"Hsp_query-from"}[0];
					if ($hit_to < $hit_from) {
						$hsp->{"Hsp_hit-to"}[0] = $hit_from;
						$hsp->{"Hsp_hit-from"}[0] = $hit_to;
						$hsp->{"Hsp_query-to"}[0] = $query_from;
						$hsp->{"Hsp_query-from"}[0] = $query_to;
						$hsp->{"Hsp_qseq"}[0] = lc(reverse_complement($hsp->{"Hsp_qseq"}[0]));
						$hsp->{"Hsp_hseq"}[0] = lc(reverse_complement($hsp->{"Hsp_hseq"}[0]));
					}
				}
				my @sorted_hsps = sort { $a->{"Hsp_hit-from"}[0] - $b->{"Hsp_hit-from"}[0] } @$hsps;
				my $query_end = 0;
				my @selected_hsps = ();
				foreach my $hsp (@sorted_hsps) { # for each hsp for this query
					if ($hsp->{"Hsp_evalue"}[0] > 0.0000001) {
						# this is no good: move on to the next hit.
						next;
					} elsif (@selected_hsps == 0) {
						# this is the first chunk of sequence.
						my @end = sort ({$a <=> $b}($query_end, $hsp->{"Hsp_query-to"}[0], $hsp->{"Hsp_query-from"}[0]));
						$query_end = pop @end;
						push @selected_hsps, $hsp;
					} elsif (($hsp->{"Hsp_query-to"}[0] <= $query_end) && ($hsp->{"Hsp_query-from"}[0] <= $query_end)) {
						# does this seq deal with a part of the query we've already addressed? Then move on.
						next;
					} else {
						# does this seq overlap with the last one in terms of hit coverage?
						my $last_hsp = pop @selected_hsps;
						if ($last_hsp->{"Hsp_hit-to"}[0] > $hsp->{"Hsp_hit-to"}[0]) {
							# if this hsp has a beginning that is smaller than the end of the last one, then we have a partial overlap.
							# compare the two by score:
							if ($last_hsp->{"Hsp_score"}[0] > $hsp->{"Hsp_score"}[0]) {
								push @selected_hsps, $last_hsp;
								next;
							} else {
								my @end = sort ({$a <=> $b}($query_end, $hsp->{"Hsp_query-to"}[0], $hsp->{"Hsp_query-from"}[0]));
								$query_end = pop @end;
								push @selected_hsps, $hsp;
								next;
							}
						} else {
							push @selected_hsps, $last_hsp;
						}
						my @end = sort ({$a <=> $b}($query_end, $hsp->{"Hsp_query-to"}[0], $hsp->{"Hsp_query-from"}[0]));
						$query_end = pop @end;
						push @selected_hsps, $hsp;
					}
				}
				my $hit_end = 1;
				my $sequence = "";
				foreach my $hsp (@selected_hsps) { # for each hsp for this query
					my $last_end = $hit_end;
					$hit_end = $hsp->{"Hsp_hit-to"}[0];
					# remove gaps from query seq that might have been inserted into the ref seq.
					while ($hsp->{"Hsp_hseq"}[0] =~ /^(.*?)(-+)(.*)$/) {
						my $left = $1;
						my $gap = $2;
						my $right = $3;
						$hsp->{"Hsp_hseq"}[0] = $left . $right;
						$hsp->{"Hsp_qseq"}[0] =~ /^(.{length($left)})(.{length($gap)})(.{length($right)})/;
						$hsp->{"Hsp_qseq"}[0] = $1 . $3;
					}
					$sequence .= "-" x ($hsp->{"Hsp_hit-from"}[0] - $last_end) . $hsp->{"Hsp_qseq"}[0];
				}
				$result_matrix{$hitdef}->{$iterquerydef} = $sequence;
			}
		}
	}
	print "result has " . (keys %result_matrix) . "\n";
	return \%result_matrix;
}

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
	system ("makeblastdb -in $align_file -dbtype nucl -out $tempreffile.db");
	system ("blastn $blast_task -query $tempreffile -db $tempreffile.db > $blastfile");
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
			system ("makeblastdb -in $tempreffile -dbtype nucl -out $tempreffile.db");
			system ("blastn $blast_task -query $align_file -db $tempreffile.db > $blastfile");
			last;
		}
	}

	open BLAST_FH, "<", $blastfile;
	my @lines = <BLAST_FH>;
	close BLAST_FH;

	my @seqids = ();
	my @seqs = ();
	my %result_matrix = ();
	my $curr_query_id = "";
	my $query_seq = "";
	my ($query_end, $subject_end) = 0;
	my ($curr_query_start, $curr_query_seq, $curr_query_end) = 0;

	if ($blast_task =~ /-task blastn-short/) {
		my $result_ptr = blast_short_to_alignment ($blastfile);
		%result_matrix = %{$result_ptr};
	} else {
		my $result_ptr = blast_to_alignment ($blastfile);
		%result_matrix = %{$result_ptr};
	}

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

