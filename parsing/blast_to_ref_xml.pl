use strict;
use XML::Simple;
require "subfuncs.pl";

my $blast_xml = shift;

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
			print ">$hitdef | $iterquerydef\n$sequence\n";
		}
	}
}
