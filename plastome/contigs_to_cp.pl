use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin/../parsing/";
use Blast qw (parse_xml revcomp_hsp);
use lib "$FindBin::Bin/../";
use Subfunctions qw (parse_fasta reverse_complement split_seq);
use File::Temp qw (tempfile);
use Data::Dumper;

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
	# set up a hash value to receive the hits when we get them.
	$region->{"hits"} = ();
	print "$region->{name}\t$region->{start}\t$region->{end}\n";
	print $fh ">$region->{name}\n$region->{sequence}\n";
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
my $hits = {};

foreach my $hit (@$hit_array) {
	my $subject = $hit->{"subject"}->{"name"};
	my $query = $hit->{"query"}->{"name"};
	foreach my $hsp (@{$hit->{"hsps"}}) {
		if ($hsp->{"hit-from"} > $hsp->{"hit-to"}) {
			$hsp = revcomp_hsp($hsp);
		}
		$hsp->{"name"} = $query;
	}
	my @hsps = sort order_by_hit_start @{$hit->{"hsps"}};
	push @{$regions_hash->{$subject}->{"hits"}}, \@hsps;
#
}

# open OUTFH, ">", $outfile or die "couldn't create $outfile";

foreach $region (@$regions) {
	print "Comparing to $region->{name}: " . @{$region->{"hits"}} . " hits\n";
	foreach my $hsp (@{$region->{"hits"}}) {
		print $hsp->{"name"} . " going from " . $hsp->{"hit-from"} . " to " . $hsp->{"hit-to"} . ", evalue " . $hsp->{"evalue"} . "\n";
	}
}

# close OUTFH;

sub order_by_hit_start {
	my $score = $b->{"hit-from"} - $a->{"hit-from"};
	my $bstart = 0;
	my $astart = 0;
	if ($b->{"hit-from"} < $b->{"hit-to"}) {
		$bstart = $b->{"hit-from"};
	} else {
		print "reversed\n";
		$bstart = $b->{"hit-to"};
	}
	if ($a->{"hit-from"} < $a->{"hit-to"}) {
		$astart = $a->{"hit-from"};
	} else {
		print "reversed\n";
		$astart = $a->{"hit-to"};
	}
	$score = $astart - $bstart;
	return $score;
}

sub order_by_query_start {
	my $score = $b->{"query-from"} - $a->{"query-from"};
	my $bstart = 0;
	my $astart = 0;
	if ($b->{"query-from"} < $b->{"query-to"}) {
		$bstart = $b->{"query-from"};
	} else {
		$bstart = $b->{"query-to"};
	}
	if ($a->{"query-from"} < $a->{"query-to"}) {
		$astart = $a->{"query-from"};
	} else {
		$astart = $a->{"query-to"};
	}
	$score = $astart - $bstart;
	return $score;
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

Aligns a list of putative cp contigs along a reference plastome.

=cut
