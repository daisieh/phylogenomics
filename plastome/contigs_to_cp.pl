use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin/../parsing/";
use Blast qw (parse_xml revcomp_hsp);
use lib "$FindBin::Bin/../";
use Subfunctions qw (parse_fasta reverse_complement split_seq find_sequences);
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

my ($ref_hash, $ref_array) = parse_fasta($reffile);
my $refseq = $ref_hash->{@$ref_array[0]};
my $reflen = length ($refseq);

print "finding inverted repeats\n";
my ($fh, $refblast) = tempfile();
system("blastn -query $reffile -subject $reffile -outfmt 5 -out $refblast.xml -evalue 1e-90");

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
	print "$region->{name}\t$region->{start}\t$region->{end}\n";
	print $fh ">$region->{name}\n$region->{sequence}\n";

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
			$hsp = revcomp_hsp($hsp);
			$contig->{"revcomp"} = " (reverse complement)";
		}
	}

	# consolidate all of the matching segments into one large overall match.
	my @hsps = sort order_by_hit_start @{$hit->{"hsps"}};
	my $first_hsp = $hsps[0];
	my $last_hsp = pop @hsps;
	my $regoffset = $regions_hash->{$region}->{"start"} - 1;
	$contig->{"hit-from"} = $first_hsp->{"hit-from"} + $regoffset;
	$contig->{"hit-to"} = $last_hsp->{"hit-to"} + $regoffset;
	$contig->{"query-from"} = $first_hsp->{"query-from"};
	$contig->{"query-to"} = $last_hsp->{"query-to"};
}

# put the sequences for the matching contigs back into the output hash.
my $contig_seqs = find_sequences ($contigfile, \@hit_list);

foreach $region (@$regions) {
	foreach my $contig (@{$region->{"hits"}}) {
		$contig->{"sequence"} = $contig_seqs->{$contig->{"name"}};
	}
	my @ordered_hits = sort order_by_hit_start @{$region->{"hits"}};
	$region->{"hits"} = \@ordered_hits;
}

# open OUTFH, ">", "$outfile.yml" or die "couldn't create $outfile.yml";

YAML::Tiny->DumpFile("$outfile.yml", $regions);

# close OUTFH;

# if $a starts earlier than $b, return -1
sub order_by_hit_start {
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
