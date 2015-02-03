#!/usr/bin/env perl
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin/../lib";
use Blast qw (parse_xml revcomp_hsp);
use Genbank qw (parse_genbank get_sequence);
use Subfunctions qw (parse_fasta reverse_complement split_seq find_sequences consensus_str);
use File::Temp qw (tempfile);
use Data::Dumper;
use YAML::Tiny;

my $help = 0;
my $outfile = "";
my $reffile = "";
my $contigfile = "";
my $join = 0;

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

GetOptions ('reffile=s' => \$reffile,
			'contigfile=s' => \$contigfile,
			'outfile=s' => \$outfile,
			'join' => \$join,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help){
    pod2usage(-verbose => 1);
}

if ($reffile eq "") {
    pod2usage(-verbose => 1, -msg => "need to specify reference file");
}

if ($contigfile eq "") {
    pod2usage(-verbose => 1, -msg => "need to specify contig file");
}

if ($outfile eq "") {
    pod2usage(-verbose => 1, -msg => "need to specify output path");
}


my $refseq = "";
if ($reffile =~ /\.gb$/) {
	my $gb = parse_genbank($reffile);
	$refseq = get_sequence($gb);
} else {
	my ($ref_hash, $ref_array) = parse_fasta($reffile);
	$refseq = $ref_hash->{@$ref_array[0]};
}
my $reflen = length ($refseq);

my ($reffh, $refseqfile) = tempfile();
print $reffh ">reference\n$refseq\n";
close $reffh;

print "finding inverted repeats\n";
my ($fh, $refblast) = tempfile();
system("blastn -query $refseqfile -subject $refseqfile -outfmt 5 -out $refblast.xml -evalue 1e-90");

my $self_array = parse_xml ("$refblast.xml");
my @irs = ();
foreach my $hit (@$self_array) {
	my @hsps = sort order_by_query_start @{$hit->{"hsps"}};
	foreach my $hsp (@hsps) {
		# only look at identical pieces that are smaller than the entire reference
		if ((($hsp->{"query-to"} - $hsp->{"query-from"}) < ($reflen - 1)) && (($hsp->{"query-to"} - $hsp->{"query-from"}) > 10000)) {
			push @irs, $hsp;
		}
	}
}

if (@irs > 2) {
	die "Error! There seem to be more than two inverted repeats. Are you sure this is a plastome sequence?";
}

my $curr_pos = 1;
my $regions = ();
my $regions_hash = {};

# LSC goes from 1 to the start of @$irs[0] - 1:
my $region = {};
$regions_hash->{"LSC"} = $region;
$region->{"name"} = "LSC";
$region->{"start"} = 1;
$region->{"end"} = $irs[0]->{"query-from"} - 1;
(undef, $region->{"sequence"}, undef) = split_seq ($refseq, $region->{"start"}, $region->{"end"});
push @$regions, $region;

# IRB goes from the start of @$irs[0] to the end of @$irs[0] (inclusive):
$region = {};
$regions_hash->{"IRB"} = $region;
$region->{"name"} = "IRB";
$region->{"sequence"} = $irs[0]{"hseq"};
$region->{"start"} = $irs[0]->{"query-from"};
$region->{"end"} = $irs[0]->{"query-to"};
push @$regions, $region;

# SSC goes from the end of @$irs[0] + 1 to the start of @$irs[1] - 1:
$region = {};
$regions_hash->{"SSC"} = $region;
$region->{"name"} = "SSC";
$region->{"start"} = $irs[0]->{"query-to"} + 1;
$region->{"end"} = $irs[1]->{"query-from"} - 1;
(undef, $region->{"sequence"}, undef) = split_seq ($refseq, $region->{"start"}, $region->{"end"});
push @$regions, $region;

# IRA goes from the start of @$irs[1] to the end of @$irs[1] (inclusive):
$region = {};
$regions_hash->{"IRA"} = $region;
$region->{"name"} = "IRA";
$region->{"sequence"} = $irs[1]{"hseq"};
$region->{"start"} = $irs[1]->{"query-from"};
$region->{"end"} = $irs[1]->{"query-to"};
push @$regions, $region;

my ($fh, $refregions) = tempfile();
foreach $region (@$regions) {
	print $fh ">" . $region->{"name"} . "\n" . $region->{"sequence"}. "\n";

	# clean up the region hash for later use.
	delete $region->{"sequence"};
	# set up a hash value to receive the hits when we get them.
	$region->{"hits"} = ();
	$region->{"length"} = $region->{"end"} - $region->{"start"} + 1;
}
close $fh;

if (-e "$outfile.xml") {
	print "skipping blastn\n";
} else {
	print "running blastn\n";
	system("blastn -query $contigfile -subject $refregions -outfmt 5 -out $outfile.xml -culling_limit 1 -evalue 1e-70");
}

print "parsing results\n";

my $hit_array = parse_xml ("$outfile.xml");
my @hit_list = ();

foreach my $hit (@$hit_array) {
	# each hit represents a contig that we want to assign to a region.
	my $contig = {};
	$contig->{"name"} = $hit->{"query"}->{"name"};
	$contig->{"length"} = $hit->{"query"}->{"length"};
	push @hit_list, $contig->{"name"};

	# push it into the appropriate region's bucket of hits.
	my $region = $hit->{"subject"}->{"name"};
	push @{$regions_hash->{$region}->{"hits"}}, $contig;
	$contig->{"region"} = $region;

	# each hsp represents a matching segment of this contig to this region.
	foreach my $hsp (@{$hit->{"hsps"}}) {
		if ($hsp->{"hit-from"} > $hsp->{"hit-to"}) {
			# tag this contig as being revcomped, so we can fix it when we deal with whole contigs.
			$contig->{"revcomp"} = " (reverse complement)";
		}
	}

	# consolidate all of the matching segments into one large overall match.
	my @query_ends = ();
	my @hit_ends = ();
	foreach my $hsp (@{$hit->{"hsps"}}) {
		push @query_ends, $hsp->{"query-from"};
		push @query_ends, $hsp->{"query-to"};
		push @hit_ends, $hsp->{"hit-from"};
		push @hit_ends, $hsp->{"hit-to"};
	}
	@query_ends = sort {$a <=> $b} @query_ends;
	@hit_ends = sort {$a <=> $b} @hit_ends;

	my $regoffset = $regions_hash->{$region}->{"start"} - 1;
	$contig->{"hit-from"} = $hit_ends[0] + $regoffset;
	$contig->{"hit-to"} = $hit_ends[@hit_ends-1] + $regoffset;
	$contig->{"query-from"} = $query_ends[0];
	$contig->{"query-to"} = $query_ends[@query_ends-1];
}

open OUTFH, ">", "$outfile.raw.yml";
print OUTFH YAML::Tiny->Dump(@$hit_array);
close OUTFH;

# put the sequences for the matching contigs back into the output hash.
my $contig_seqs = find_sequences ($contigfile, \@hit_list);
# write these best seqs out:
open OUTFH, ">", "$outfile.best.fasta";
foreach my $key (keys %$contig_seqs) {
	print "$key\n";
	print OUTFH ">$key\n";
	print OUTFH $contig_seqs->{$key} . "\n";
}
close OUTFH;

my @all_hits = ();
foreach $region (@$regions) {
	foreach my $contig (@{$region->{"hits"}}) {
		$contig->{"sequence"} = $contig_seqs->{$contig->{"name"}};
		if (exists $contig->{"revcomp"}) {
			delete $contig->{"revcomp"};
			$contig->{"sequence"} = reverse_complement ($contig->{"sequence"});
			$contig->{"name"} .= "_rc";
			# flip the query's indices: the q-from is now going to be (length - q-from) and the q-to is (length - q-to)
			my $old_qto = $contig->{"query-to"};
			$contig->{"query-to"} = $contig->{"length"} - $contig->{"query-from"};
			$contig->{"query-from"} = $contig->{"length"} - $old_qto;
		}

		# do some cleanup of the hit and query windows.
		# each contig's putative hit span is from the amount of query extending before the start of the match (hit-from), which is (hit-from - query-from), plus whatever portion of its length is beyond that (length - query-from + hit-from)
		$contig->{"hit-to"} = $contig->{"hit-from"} - $contig->{"query-from"} + $contig->{"length"};
		$contig->{"hit-from"} = $contig->{"hit-from"} - $contig->{"query-from"};
# 		print "cleaning up " . $contig->{"name"} . ", has length " . $contig->{"length"} . " and its " . $contig->{"query-from"} . "-" . $contig->{"query-to"} . " covers the ref " . $contig->{"hit-from"} . "-" . $contig->{"hit-to"} . "\n";
	}
	my @ordered_hits = sort order_by_hit_start @{$region->{"hits"}};
	push @all_hits, @ordered_hits;
	$region->{"hits"} = \@ordered_hits;
}

open OUTFH, ">", "$outfile.yml";
print OUTFH YAML::Tiny->Dump(@$regions);
close OUTFH;

# do the contigs connect to each other?
my @final_contigs = ();
my $first_hit = shift @all_hits;
push @final_contigs, $first_hit;

# compare the end of the last contig in final_contigs to the start of the next contig in all_hits
# while there's still anything left in all_hits
while (@all_hits > 0) {
	# first contig to compare is the last one in final_contigs
	my $contig_seq1 = @final_contigs[@final_contigs - 1];

	# second contig to compare is the next unanalyzed one from all_hits
	my $contig_seq2 = shift @all_hits;

	print "comparing " . $contig_seq1->{"name"} . ", maps " . $contig_seq1->{"hit-from"}. "-" . $contig_seq1->{"hit-to"} . ", to " . $contig_seq2->{"name"} . ", maps " . ($contig_seq2->{"hit-to"} - $contig_seq2->{"length"}) ."-". $contig_seq2->{"hit-to"} . "\n";

	# if the second contig's putative hit range is within the first, drop it.
	if ($contig_seq2->{"hit-to"} <= $contig_seq1->{"hit-to"}) {
		next;
	}

	# compare these two contigs' ends
	print "can we meld these contigs? ";
	(my $fh1, my $contig1) = tempfile();
	my $contig_end = $contig_seq1->{"sequence"};
	if ($contig_seq1->{"sequence"} =~ /^(.*)(.{50})$/) {
		$contig_seq1->{"sequence"} = $1;
		$contig_end = $2;
	}
	print $fh1 ">" . $contig_seq1->{"name"} . "_end\n$contig_end\n";
	my $contig_start = $contig_seq2->{"sequence"};
	if ($contig_seq2->{"sequence"} =~ /^(.{50})/) {
		$contig_start = $1;
	}
	print $fh1 ">" . $contig_seq2->{"name"} . "_start\n$contig_start\n";
	close $fh1;
	(undef, my $temp_out) = tempfile(OPEN => 0);

	system ("mafft --retree 2 --maxiterate 0 --op 10 $contig1 > $temp_out 2>/dev/null");
	# the resulting sequences are a match if the consensus sequence has few ambiguous characters.
	(my $aligned_bits, my $alignarray) = parse_fasta($temp_out);
	my @seqs = ();
	foreach my $k (@$alignarray) {
		push @seqs, $aligned_bits->{$k};
	}
	my $cons_seq = consensus_str(\@seqs);
	my @ambigs = $cons_seq =~ m/[NMRWSYKVHDB]/g;
	if (@ambigs < 5) { # if there are less than 5 ambiguities when we align them...

		# meld the sequence:
		$contig_seq1->{"sequence"} = $contig_seq1->{"sequence"} . $cons_seq . $contig_seq2->{"sequence"};

		# update the parameters for the newly-melded contig.
		$contig_seq1->{"name"} = $contig_seq1->{"name"} . "+" . $contig_seq2->{"name"};
		$contig_seq1->{"region"} = $contig_seq1->{"region"} . "+" . $contig_seq2->{"region"};
		$contig_seq1->{"length"} = length $contig_seq1->{"sequence"};
		$contig_seq1->{"hit-to"} = $contig_seq2->{"hit-to"};
		$contig_seq1->{"query-to"} = "" . ($contig_seq1->{"length"} - $contig_seq1->{"query-from"});
		print "yes, new length is " . $contig_seq1->{"length"} . ", covers " . $contig_seq1->{"hit-from"} . "-" . $contig_seq1->{"hit-to"} . "\n";
	} elsif (($contig_seq2->{"hit-from"} - $contig_seq1->{"hit-to"}) < 300) {
		# if the two contigs' hit ends are within 100 bp of each other, scaffold them together by adding Ns
		$contig_seq1->{"sequence"} .= "N" x ($contig_seq2->{"hit-from"} - $contig_seq1->{"hit-to"}) . $contig_seq2->{"sequence"};

		# update the parameters for the newly-melded contig.
		$contig_seq1->{"name"} = $contig_seq1->{"name"} . "+" . $contig_seq2->{"name"};
		$contig_seq1->{"region"} = $contig_seq1->{"region"} . "+" . $contig_seq2->{"region"};
		$contig_seq1->{"length"} = length $contig_seq1->{"sequence"};
		$contig_seq1->{"hit-to"} = $contig_seq2->{"hit-to"};
		$contig_seq1->{"query-to"} = "" . ($contig_seq1->{"length"} - $contig_seq1->{"query-from"});
		print "maybe? new length is " . $contig_seq1->{"length"} . ", covers " . $contig_seq1->{"hit-from"} . "-" . $contig_seq1->{"hit-to"} . "\n";
	} else {
		# if not meldable, push the second contig onto the final contigs too.
		$contig_seq1->{"sequence"} .= $contig_end;
		$contig_seq2->{"sequence"} = $contig_start . $contig_seq2->{"sequence"};
		push @final_contigs, $contig_seq2;
		print "no, span too large: " . $contig_seq2->{"hit-from"} ."-". $contig_seq1->{"hit-to"} . " = ". ($contig_seq2->{"hit-from"} - $contig_seq1->{"hit-to"}) . "\n";
	}
}

my $final_len = 0;
foreach my $c (@final_contigs) {
	$final_len += $c->{"length"};
}
print "final assembly has " . @final_contigs . " contigs, total length $final_len\n";

open OUTFH, ">", "$outfile.final.yml";
print OUTFH YAML::Tiny->Dump(@final_contigs);
close OUTFH;

open OUTFH, ">", "$outfile.draft.fasta";
if ($join) {
	print OUTFH ">$outfile.draft.fasta\n";
	foreach my $c (@final_contigs) {
		print OUTFH $c->{"sequence"} . "NNNNNNNNNNN";
	}
} else {
	foreach my $c (@final_contigs) {
		print OUTFH ">" . $c->{"name"} . "\n" . $c->{"sequence"} . "\n";
	}
}
close OUTFH;

### Sorting functions

# if $a starts earlier than $b, return -1
sub order_by_hit_start {
	my $bstart = $b->{"hit-from"};
	my $astart = $a->{"hit-from"};

	if ($astart < $bstart) { return -1; }
	if ($astart > $bstart) { return 1; }
	return 0;
}

sub order_by_ref_start {
	my $bstart = $b->{"hit-from"};
	my $astart = $a->{"hit-from"};

	if ($astart < $bstart) { return -1; }
	if ($astart > $bstart) { return 1; }
	return 0;
}

sub order_by_query_start {
	my $bstart = $b->{"query-from"};
	my $astart = $a->{"query-from"};

	if ($astart < $bstart) { return -1; }
	if ($astart > $bstart) { return 1; }
	return 0;
}


__END__

=head1 NAME

contigs_to_cp.pl

=head1 SYNOPSIS

contigs_to_cp.pl [-reffile reffile] [-contigfile contigfile] [-outputfile outputfile]

=head1 OPTIONS

  -reffile:         genbank or fasta file of reference plastome
  -contigfile:      fasta file of putative cp contigs
  -outputfile:      name of output file

=head1 DESCRIPTION

Aligns a list of putative cp contigs along a reference plastome. Outputs a YAML file of the best-matching contigs, in order.

=cut
