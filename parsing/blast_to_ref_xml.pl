use strict;
use XML::Simple;
use Data::Dumper;

my $blast_xml = shift;

my $parser = new XML::Simple();

open XML_FH, "<", $blast_xml;
my %xml_strings = ();

my $curr_xml = readline XML_FH;
my $query_name = "";

while (my $line = readline XML_FH) {
	if ($line =~ /<\?xml/) {
		# we've reached the next xml object
		$xml_strings{$query_name} = "$curr_xml";
		$curr_xml = "$line";
	} else {
		if ($line =~ /<BlastOutput_query-def>(.*?)<\/BlastOutput_query-def>/) {
			$query_name = $1;
		}
		$curr_xml .= $line;
	}
}
$xml_strings{$query_name} = "$curr_xml";

close XML_FH;

foreach my $key (keys %xml_strings) {
	print "key $key\n";
	my $tree = $parser->XMLin($xml_strings{$key});

	my $blastoutput_iterations = $tree->{"BlastOutput_iterations"}; # BlastOutput_iterations is a hash with one key, "Iterations"
	my $iterations = ($tree->{"BlastOutput_iterations"}->{"Iteration"}); # key Iteration has a value that is an anonymous array of iteration hashes.
	foreach my $iteration (@$iterations) { # for each iteration
		my $iterid = $iteration->{"Iteration_query-ID"};
		my $iterquerydef = $iteration->{"Iteration_query-def"};
		my $hsps = $iteration->{"Iteration_hits"}->{"Hit"}->{"Hit_hsps"}->{"Hsp"}; # hsps should be an anonymous array of the actual hit sequences.
		my @sorted_hsps = sort { $a->{"Hsp_hit-from"} - $b->{"Hsp_hit-from"} } @$hsps;
		foreach my $hsp (@sorted_hsps) { # for each hsp for this query
			print "## $iterid ($iterquerydef) from query " . $hsp->{"Hsp_query-from"} . "-" . $hsp->{"Hsp_query-to"} . " has an hsp of score " . $hsp->{"Hsp_score"} ."\n";
			print "##    which matches " . $hsp->{"Hsp_hit-from"} . "-" . $hsp->{"Hsp_hit-to"} ."\n";
			print "\t" . $hsp->{"Hsp_qseq"} . "\n";
			print "\t" . $hsp->{"Hsp_hseq"} . "\n";
			if ($hsp->{"Hsp_evalue"} < 0.0000001) {

			}
		}
	}

	print Dumper ($tree) . "\n";
}

my @element_stack = ();
sub sort_hsps {
	return ( $a->{"Hsp_query-from"} - $b->{"Hsp_query-from"} );
}
