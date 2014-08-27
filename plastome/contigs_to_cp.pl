use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin/../parsing/";
use Blast qw (parse_xml revcomp_hsp);
use Genbank qw (parse_genbank get_sequence);
use lib "$FindBin::Bin/../";
use Subfunctions qw (parse_fasta reverse_complement split_seq find_sequences consensus_str);
use File::Temp qw (tempfile);
use Data::Dumper;
use YAML::Tiny;

my $help = 0;
my $outfile = "";
my $reffile = "";
my $contigfile = "";


GetOptions ('reffile=s' => \$reffile,
			'contigfile=s' => \$contigfile,
			'outfile=s' => \$outfile,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

my $refseq = "";
if ($reffile =~ /\.gb$/) {
	my $gb = parse_genbank($reffile);
	$refseq = get_sequence($gb);
	print "$refseq\n";
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
		if (($hsp->{"query-to"} - $hsp->{"query-from"}) < ($reflen - 1)) {
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
	}
	my @ordered_hits = sort order_by_ref_start @{$region->{"hits"}};
	push @all_hits, @ordered_hits;
	$region->{"hits"} = \@ordered_hits;
}

# clean up and trim sequences.
# first, is the first LSC contig extending into the previous IR?
# my $LSC_start = $all_hits[0];
# if ($LSC_start->{"hit-from"} == 1) {
# 	print "trimming start of LSC\n";
# 	(undef, $LSC_start->{"sequence"}, undef) = split_seq ($LSC_start->{"sequence"}, $LSC_start->{"query-from"}, $LSC_start->{"query-to"});
# 	my $offset = $LSC_start->{"query-from"} - 1;
# 	$LSC_start->{"query-from"} = $LSC_start->{"query-from"} - $offset;
# 	$LSC_start->{"query-to"} = $LSC_start->{"query-to"} - $offset;
# }

open OUTFH, ">", "$outfile.yml";
print OUTFH YAML::Tiny->Dump($regions);
close OUTFH;

# do the contigs connect to each other?
my @final_contigs = ();
my $first_hit = shift @all_hits;
push @final_contigs, $first_hit;
# print "final_contigs " . Dumper(@final_contigs) . "\n";

# for (my $i=0; $i < (@all_hits -1); $i++) {

while (@all_hits > 0) { # while there's still anything left in all_hits
	# compare the end of the last contig in final_contigs to the start of the next contig in all_hits
	my $contig_seq1 = pop @final_contigs;
	(my $fh1, my $contig1) = tempfile();
	my $contig_end = $contig_seq1->{"sequence"};
	if ($contig_seq1->{"sequence"} =~ /^(.*)(.{50})$/) {
		$contig_seq1->{"sequence"} = $1;
		$contig_end = $2;
	}
	print $fh1 ">" . $contig_seq1->{"name"} . "_end\n$contig_end\n";
	close $fh1;
#
	my $contig_seq2 = shift @all_hits;
	(my $fh2, my $contig2) = tempfile();
	my $contig_start = $contig_seq2->{"sequence"};
	if ($contig_seq2->{"sequence"} =~ /^(.{50})/) {
		$contig_start = $1;
	}
	print $fh2 ">" . $contig_seq2->{"name"} . "_start\n$contig_start\n";
	close $fh2;
#
	(undef, my $temp_out) = tempfile(OPEN => 0);
	system ("blastn -query $contig1 -subject $contig2 -word_size 20 -outfmt 5 -out $temp_out");
	my $contig_hits = parse_xml ("$temp_out");

	# if there is a hit and the score is high enough, meld these two contigs and push that conjoined contig onto the final contigs.
	if ((@$contig_hits > 0) && (@$contig_hits[0]->{"hsps"}[0]->{"score"} >= 20)) {
		my $this_hit = @$contig_hits[0]->{"hsps"}[0];
		my $offset = ($this_hit->{"query-from"} - 1);
		my $hit_offset = $offset - ($this_hit->{"hit-from"} - 1);
		$contig_seq1->{"name"} = $contig_seq1->{"name"} . "+" . $contig_seq2->{"name"};
		$contig_seq1->{"region"} = $contig_seq1->{"region"} . "+" . $contig_seq2->{"region"};
		$contig_start = "-" x $hit_offset . "$contig_start";
		$contig_seq1->{"sequence"} .= consensus_str([$contig_end, $contig_start]) . $contig_seq2->{"sequence"};
		$contig_seq1->{"length"} = length $contig_seq1->{"sequence"};
		$contig_seq1->{"hit-to"} = $contig_seq2->{"hit-to"};
		$contig_seq1->{"query-to"} = $contig_seq2->{"query-to"};
		push @final_contigs, $contig_seq1;
	} else {
		$contig_seq1->{"sequence"} .= $contig_end;
		$contig_seq2->{"sequence"} = $contig_start . $contig_seq2->{"sequence"};
		push @final_contigs, $contig_seq1;
		push @final_contigs, $contig_seq2;
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
foreach my $c (@final_contigs) {
	print OUTFH ">" . $c->{"name"} . "\n" . $c->{"sequence"} . "\n";
}
close OUTFH;


# confirm ordering of outputs
(undef, my $temp_out) = tempfile(OPEN => 0);
system ("blastn -query $reffile -subject $outfile.draft.fasta -outfmt 5 -out $temp_out");
my $res = parse_xml($temp_out);
foreach my $r (@$res) {
	my @new_hsps = ();
	foreach my $y (@{$r->{"hsps"}}) {
		if ($y->{"hit-frame"} == 1) {
			push @new_hsps, $y;
		}
	}
	my @sorted_hsps = sort order_by_query_start @new_hsps;
	$r->{"hsps"} = \@sorted_hsps;
}
print YAML::Tiny->Dump($res) . "\n";

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

  -reffile:         fasta file of reference plastome
  -contigfile:      fasta file of putative cp contigs
  -outputfile:      name of output file

=head1 DESCRIPTION

Aligns a list of putative cp contigs along a reference plastome. Outputs a YAML file of the best-matching contigs, in order.

=cut
